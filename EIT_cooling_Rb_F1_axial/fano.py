import os
import numpy as np
from qutip import basis, steadystate, expect
import config as cfg


def run_fano():
    # --- 1. DIRECTORY SETUP ---
    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.join(script_dir, "results_fano", cfg.RUN_NAME)
    os.makedirs(results_dir, exist_ok=True)

    print(f"Starting Fano steady-state simulation for: {cfg.RUN_NAME}")

    # --- 2. INDEX MAPPING ---
    g1_idxs = [i for i, l in enumerate(cfg.atom_labels) if l[0] == "g1"]
    g2_idxs = [i for i, l in enumerate(cfg.atom_labels) if l[0] == "g2"]
    e_idxs = [i for i, l in enumerate(cfg.atom_labels) if l[0].startswith("e")]

    e3_idxs = [i for i, l in enumerate(cfg.atom_labels) if l[1] == 3]
    e2_idxs = [i for i, l in enumerate(cfg.atom_labels) if l[0] == "e2"]

    # --- 3. PROBE SCAN RANGE ---
    scan_width = 0.15 * cfg.gamma
    det_list = np.linspace(cfg.Delta_p_center - scan_width, cfg.Delta_p_center + scan_width, 300)

    # --- 4. NORMALIZATION TARGETS ---
    # Each laser's declared Rabi frequency (Omega_c_amp, Omega_p_amp, Omega_r_amp)
    # refers to its own *designed* transition (Fig. 1 scheme). Every other
    # dipole-allowed transition driven by the same beam/polarization is scaled
    # relative to that target via the properly weighted hyperfine dipole ratio
    # C_eg^(q) = d_eg^(q) / d_target, matching Eq. (general_laser_interaction)
    # ("normalized transition amplitude").
    _, Fg_c, mg_c = cfg.COUPLING_GROUND
    _, Fe_t, me_t = cfg.EXCITED_TARGET
    d_c_target = cfg.dipole_element(Fg_c, mg_c, Fe_t, me_t, cfg.pol_c)

    _, Fg_p, mg_p = cfg.PROBE_GROUND
    d_p_target = cfg.dipole_element(Fg_p, mg_p, Fe_t, me_t, cfg.pol_p)

    _, Fg_r, mg_r = cfg.REPUMP_GROUND
    d_r_target = cfg.dipole_element(Fg_r, mg_r, Fe_t, me_t, cfg.pol_r)

    # --- 5. BUILDING THE CONSTANT HAMILTONIAN (H_const) ---
    H_const = 0
    for i, (label, f, m) in enumerate(cfg.atom_labels):
        if label.startswith("e"):
            energy = (getattr(cfg, f"E_{label}") - cfg.E_e2) - cfg.Delta_c
        else:
            energy = 0.0

        g_f = -0.5 if "g1" in label else 0.5 if "g2" in label else 0.67
        zeeman_shift = g_f * cfg.mu_B * cfg.B_field * m
        H_const += (energy + zeeman_shift) * (basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag())

    # Coupling laser (sigma-), F=1 -> e. FIX: uses the hyperfine-weighted
    # dipole_element instead of a bare Clebsch-Gordan coefficient.
    for gi in g1_idxs:
        Fg_i, mg_i = cfg.atom_labels[gi][1], cfg.atom_labels[gi][2]
        for ei in e_idxs:
            Fe_i, me_i = cfg.atom_labels[ei][1], cfg.atom_labels[ei][2]
            d = cfg.dipole_element(Fg_i, mg_i, Fe_i, me_i, cfg.pol_c)
            if abs(d) > 1e-8:
                c_rel = d / d_c_target
                V = (cfg.Omega_c_amp * c_rel / 2.0) * basis(cfg.N_atom, ei) * basis(cfg.N_atom, gi).dag()
                H_const += V + V.dag()

    # --- Resonant Repumper Laser ---
    # FIX (scheme swap): now pi-polarized, drives F=2, m=-2 -> F'=2, m=-2
    # (was sigma- on F=2, m=-1 -> F'=2, m=-2), matching Fig. 1 of the thesis.
    try:
        g_rep_idx = cfg.atom_labels.index(cfg.REPUMP_GROUND)
        e_rep_idx = cfg.atom_labels.index(cfg.EXCITED_TARGET)

        d_rep = cfg.dipole_element(Fg_r, mg_r, Fe_t, me_t, cfg.pol_r)

        if abs(d_rep) > 1e-8:
            c_rel = d_rep / d_r_target  # = 1.0 by construction (this IS the target)
            V_rep = (cfg.Omega_r_amp * c_rel / 2.0) * basis(cfg.N_atom, e_rep_idx) * basis(cfg.N_atom, g_rep_idx).dag()
            H_const += V_rep + V_rep.dag()

    except ValueError:
        print("Error: Repumper states not found in atom_labels.")

    # --- 6. COLLAPSE OPERATORS (Spontaneous Emission) ---
    # FIX: uses the normalized branching fraction b_eg (Eq. branching_fraction)
    # instead of a bare Clebsch-Gordan coefficient, so that
    # sum_g L_{e->g}^dagger L_{e->g} = gamma * |e><e| exactly for every
    # excited sublevel (previously F'=1,2 decayed twice as fast as F'=0,3).
    c_ops = []
    for ei in e_idxs:
        Fe_i, me_i = cfg.atom_labels[ei][1], cfg.atom_labels[ei][2]
        for gi in g1_idxs + g2_idxs:
            Fg_i, mg_i = cfg.atom_labels[gi][1], cfg.atom_labels[gi][2]
            b = cfg.branching_fraction(Fg_i, mg_i, Fe_i, me_i)
            if b > 1e-10:
                c_ops.append(np.sqrt(cfg.gamma * b) * (basis(cfg.N_atom, gi) * basis(cfg.N_atom, ei).dag()))

    # --- 7. PROJECTION OPERATORS & DICTIONARY SETUP ---
    P_e_total = sum([basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag() for i in e_idxs])
    P_e2_only = sum([basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag() for i in e2_idxs])
    P_e3_leakage = sum([basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag() for i in e3_idxs])

    results = {
        "det_list": det_list, "abs_tot": [], "abs_e2": [], "leak_e3": []
    }

    for label, f, m in cfg.atom_labels:
        results[f"pop_{label}_m{m}"] = []

    # --- 8. MAIN SIMULATION LOOP ---
    # Probe laser (sigma-). FIX: was pi-polarized; now correctly sigma-
    # (Delta m = -1), matching Fig. 1: |F=2, m=-1> -> |F'=2, m=-2>.
    for dp in det_list:
        H_probe = 0
        two_photon_det = dp - cfg.Delta_c

        for gi in g2_idxs:
            H_probe += two_photon_det * (basis(cfg.N_atom, gi) * basis(cfg.N_atom, gi).dag())
            Fg_i, mg_i = cfg.atom_labels[gi][1], cfg.atom_labels[gi][2]
            for ei in e_idxs:
                Fe_i, me_i = cfg.atom_labels[ei][1], cfg.atom_labels[ei][2]
                d = cfg.dipole_element(Fg_i, mg_i, Fe_i, me_i, cfg.pol_p)
                if abs(d) > 1e-8:
                    c_rel = d / d_p_target
                    V = (cfg.Omega_p_amp * c_rel / 2.0) * basis(cfg.N_atom, ei) * basis(cfg.N_atom, gi).dag()
                    H_probe += V + V.dag()

        rho_ss = steadystate(H_const + H_probe, c_ops)

        results["abs_tot"].append(expect(P_e_total, rho_ss))
        results["abs_e2"].append(expect(P_e2_only, rho_ss))
        results["leak_e3"].append(expect(P_e3_leakage, rho_ss))

        for i, (label, f, m) in enumerate(cfg.atom_labels):
            P_op = basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag()
            results[f"pop_{label}_m{m}"].append(expect(P_op, rho_ss))

    # --- 9. SAVE RESULTS ---
    for key in results:
        np.save(os.path.join(results_dir, f"{key}.npy"), np.array(results[key]))

    print(f"Success. All data saved in results_fano/{cfg.RUN_NAME}")


if __name__ == "__main__":
    run_fano()
