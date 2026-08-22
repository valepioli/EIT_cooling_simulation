import os
import numpy as np
import qutip
from qutip import basis, tensor, qeye, destroy, fock, mcsolve, qload, qsave, expect
import config as cfg

# FIX (QuTiP 4/5 compatibility): QuTiP 5 removed the old `Options` class
# (solver options are now a plain dict) and `qutip.ui.progressbar.TextProgressBar`
# (the progress bar is now selected with a string, e.g. "text"). Depending on
# how QuTiP was installed -- pip (`pip install qutip`, usually >=5.x) vs. a
# Linux distro package (e.g. `/usr/lib/python3/dist-packages/qutip`, often
# still 4.x) -- either API may be present, and each raises a hard error under
# the other's calling convention:
#   - QuTiP 5 code path under QuTiP 4:  ImportError: cannot import name 'Options'
#   - QuTiP 4 code path under QuTiP 5:  AttributeError: 'dict' object has no
#     attribute 'rhs_reuse' (mcsolve tries to read Options-style attributes
#     off the plain dict)
# `requirements.txt` did not pin a version, so both situations were possible
# depending on the environment. The helper below detects the installed
# major version at runtime and builds the right kind of solver options /
# mcsolve call for it, so the script works unmodified on either QuTiP 4.x
# or 5.x.
QUTIP_MAJOR = int(qutip.__version__.split(".")[0])


def make_mcsolve_call(H, psi0, tlist, c_ops, e_ops, ntraj, nsteps=500000, atol=1e-7, rtol=1e-5):
    """Version-adaptive wrapper around qutip.mcsolve (QuTiP 4.x and 5.x)."""
    if QUTIP_MAJOR >= 5:
        solver_options = {
            "store_states": False,
            "nsteps": nsteps,
            "atol": atol,
            "rtol": rtol,
            "progress_bar": "text",
        }
        return mcsolve(H, psi0, tlist, c_ops, e_ops=e_ops, ntraj=ntraj, options=solver_options)
    else:
        from qutip import Options
        solver_options = Options(
            store_states=False,
            nsteps=nsteps,
            atol=atol,
            rtol=rtol,
        )
        return mcsolve(H, psi0, tlist, c_ops, e_ops=e_ops, ntraj=ntraj,
                        options=solver_options, progress_bar=True)

# --- SIMULATION PARAMETERS ---
# Raised ceiling: Prevents wave function truncation (resolves NaN crashes)
N_vib = 4
eta = cfg.eta
dp = cfg.Delta_p_center
n_phon = cfg.n_phon

# Monte Carlo time and trajectory settings
t_total = 700000        # Requested total time
n_points = 5000
n_traj = 150

# --- PATH SETUP ---
script_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(script_dir, "results_time")
os.makedirs(results_dir, exist_ok=True)
save_file = os.path.join(results_dir, f"time_evol_MC_{cfg.RUN_NAME}")


def build_system():
    print(f"Building Tensor Hamiltonian (Atom levels x {N_vib} Fock states)...")

    a = tensor(qeye(cfg.N_atom), destroy(N_vib))
    n_op = a.dag() * a

    g1_idxs = [i for i, l in enumerate(cfg.atom_labels) if l[0] == "g1"]
    g2_idxs = [i for i, l in enumerate(cfg.atom_labels) if l[0] == "g2"]
    e_idxs = [i for i, l in enumerate(cfg.atom_labels) if l[0].startswith("e")]

    # --- normalization targets (see fano.py for the same convention) ---
    _, Fg_c, mg_c = cfg.COUPLING_GROUND
    _, Fe_t, me_t = cfg.EXCITED_TARGET
    d_c_target = cfg.dipole_element(Fg_c, mg_c, Fe_t, me_t, cfg.pol_c)

    _, Fg_p, mg_p = cfg.PROBE_GROUND
    d_p_target = cfg.dipole_element(Fg_p, mg_p, Fe_t, me_t, cfg.pol_p)

    _, Fg_r, mg_r = cfg.REPUMP_GROUND
    d_r_target = cfg.dipole_element(Fg_r, mg_r, Fe_t, me_t, cfg.pol_r)

    H_const_atom = 0 * basis(cfg.N_atom, 0) * basis(cfg.N_atom, 0).dag()

    # Calculate energy levels and Zeeman shifts
    for i, (label, f, m) in enumerate(cfg.atom_labels):
        if label.startswith("e"):
            energy = (getattr(cfg, f"E_{label}") - cfg.E_e2) - cfg.Delta_c
        else:
            energy = 0.0

        g_f = -0.5 if "g1" in label else 0.5 if "g2" in label else 0.67
        zeeman = g_f * cfg.mu_B * cfg.B_field * m
        H_const_atom += (energy + zeeman) * (basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag())

    # Coupling laser (sigma-), F=1 -> e.
    # FIX: uses the hyperfine-weighted dipole_element (Eq. dipole_matrix_multilevel)
    # instead of a bare Clebsch-Gordan coefficient (cfg.safe_clebsch).
    for gi in g1_idxs:
        Fg_i, mg_i = cfg.atom_labels[gi][1], cfg.atom_labels[gi][2]
        for ei in e_idxs:
            Fe_i, me_i = cfg.atom_labels[ei][1], cfg.atom_labels[ei][2]
            d = cfg.dipole_element(Fg_i, mg_i, Fe_i, me_i, cfg.pol_c)
            if abs(d) > 1e-8:
                c_rel = d / d_c_target
                V = (cfg.Omega_c_amp * c_rel / 2.0) * basis(cfg.N_atom, ei) * basis(cfg.N_atom, gi).dag()
                H_const_atom += V + V.dag()

    H_const = tensor(H_const_atom, qeye(N_vib))

    # Probe laser (sigma-), F=2 -> e, carries the Lamb-Dicke motional coupling.
    # FIX (scheme swap): was pi-polarized (pol_p=0) on the whole F=2 manifold;
    # now correctly sigma- (pol_p=-1), matching Fig. 1: |F=2,m=-1> -> |F'=2,m=-2>.
    # The coupling and repump beams keep eta=0 (motional factor set to the
    # identity), exactly as stated in the thesis for the repumper
    # (Eq. repump_no_motion) and for the coupling beam in the three-level
    # benchmark (orthogonal propagation, eta_c=0).
    H_probe = 0
    two_photon_det = dp - cfg.Delta_c

    H_probe_det_atom = sum([two_photon_det * (basis(cfg.N_atom, gi) * basis(cfg.N_atom, gi).dag()) for gi in g2_idxs])
    H_probe += tensor(H_probe_det_atom, qeye(N_vib))

    LD_term_plus = qeye(N_vib) + 1j * eta * (destroy(N_vib) + destroy(N_vib).dag())
    LD_term_minus = qeye(N_vib) - 1j * eta * (destroy(N_vib) + destroy(N_vib).dag())

    for gi in g2_idxs:
        Fg_i, mg_i = cfg.atom_labels[gi][1], cfg.atom_labels[gi][2]
        for ei in e_idxs:
            Fe_i, me_i = cfg.atom_labels[ei][1], cfg.atom_labels[ei][2]
            d = cfg.dipole_element(Fg_i, mg_i, Fe_i, me_i, cfg.pol_p)
            if abs(d) > 1e-8:
                c_rel = d / d_p_target
                V_atom = (cfg.Omega_p_amp * c_rel / 2.0) * basis(cfg.N_atom, ei) * basis(cfg.N_atom, gi).dag()
                H_probe += tensor(V_atom, LD_term_plus) + tensor(V_atom.dag(), LD_term_minus)

    # Trap Hamiltonian
    H_trap = cfg.nu * n_op

    # --- Resonant repumper (pi-polarized) ---
    # FIX (scheme swap): was sigma- on |F=2,m=-1> -> |F'=2,m=-2>; now
    # correctly pi-polarized on |F=2,m=-2> -> |F'=2,m=-2> (Fig. 1). No
    # motional factor (Eq. repump_no_motion), unchanged.
    H_repump_atom = 0 * basis(cfg.N_atom, 0) * basis(cfg.N_atom, 0).dag()
    try:
        g_rep_idx = cfg.atom_labels.index(cfg.REPUMP_GROUND)
        e_rep_idx = cfg.atom_labels.index(cfg.EXCITED_TARGET)

        d_rep = cfg.dipole_element(Fg_r, mg_r, Fe_t, me_t, cfg.pol_r)

        if abs(d_rep) > 1e-8:
            c_rel = d_rep / d_r_target  # = 1.0 by construction
            V_rep = (cfg.Omega_r_amp * c_rel / 2.0) * basis(cfg.N_atom, e_rep_idx) * basis(cfg.N_atom, g_rep_idx).dag()
            H_repump_atom += V_rep + V_rep.dag()

    except ValueError:
        print("Error: Repumper states not found in atom_labels.")

    H_repump = tensor(H_repump_atom, qeye(N_vib))

    # Total Hamiltonian
    H_tot = H_const + H_probe + H_trap + H_repump

    # Collapse operators for spontaneous emission.
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
                decay_op = np.sqrt(cfg.gamma * b) * (basis(cfg.N_atom, gi) * basis(cfg.N_atom, ei).dag())
                c_ops.append(tensor(decay_op, qeye(N_vib)))

    return H_tot, c_ops, n_op


def run_monte_carlo_evolution():
    H, c_ops, n_op = build_system()

    print("\n[+] System built successfully. Initializing Monte Carlo parameters...")

    # --- CREATION OF PROJECTORS FOR EACH ATOMIC STATE ---
    state_projectors = []
    state_labels = []
    for i, (label, f, m) in enumerate(cfg.atom_labels):
        P_atom = basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag()
        P_tot = tensor(P_atom, qeye(N_vib))

        state_projectors.append(P_tot)
        state_labels.append(f"{label}_F{f}_m{m}")

    # Operators to measure: phonon number + all state projectors
    e_ops = [n_op] + state_projectors

    # FIX (scheme swap): the atom now starts in the designed *probe* ground
    # state |F=2, m=-1> (cfg.PROBE_GROUND), not |F=2, m=-2> (which is now
    # the repump leak state, not part of the closed Lambda cycle).
    ideal_g_probe = cfg.atom_labels.index(cfg.PROBE_GROUND)

    # STARTS EXACTLY FROM n PHONONS HERE
    psi0 = tensor(basis(cfg.N_atom, ideal_g_probe), fock(N_vib, n_phon))

    t_list_full = np.linspace(0.0, t_total, n_points)

    print(f"[QuTiP {qutip.__version__}, using the {'>=5.x dict-options' if QUTIP_MAJOR >= 5 else '4.x Options()'} API]")
    print(f"Target time: {t_total} | Calculating {n_traj} parallel trajectories...")

    # nsteps increased to cover long segments between jumps; atol/rtol
    # balanced to be tolerant enough for speed while remaining safe.
    result = make_mcsolve_call(H, psi0, t_list_full, c_ops, e_ops, n_traj,
                                nsteps=500000, atol=1e-7, rtol=1e-5)

    # Extracting results
    n_expect_full = result.expect[0]

    populations = {}
    for idx, lbl in enumerate(state_labels):
        # Index is idx + 1 because result.expect[0] is occupied by n_op
        populations[lbl] = result.expect[idx + 1]

    data_to_save = {
        't_list': t_list_full,
        'n_expect': n_expect_full,
        'populations': populations,
        'params': {'N_vib': N_vib, 'eta': eta, 'dp': dp, 'ntraj': n_traj}
    }

    qsave(data_to_save, save_file)
    print(f"\n[+] Saved results to {save_file}.qu")
    print(f"--- SIMULATION COMPLETE. Final <n> = {n_expect_full[-1]:.3f} ---")


if __name__ == "__main__":
    run_monte_carlo_evolution()