import os
import numpy as np
from qutip import basis, steadystate, expect
import config as cfg  # Importing parameters and helper functions from config.py

def run_fano():
    # --- 1. DIRECTORY SETUP ---
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Files will be saved in: results_fano / <RUN_NAME>
    results_dir = os.path.join(script_dir, "results_fano", cfg.RUN_NAME)
    os.makedirs(results_dir, exist_ok=True)

    print(f"Starting calculation fano simulation for: {cfg.RUN_NAME}")

    # --- 2. INDEX MAPPING ---
    # We group the 24 atomic levels by their physical characteristics for easier calculations
    g1_idxs = [i for i, l in enumerate(cfg.atom_labels) if l[0] == "g1"]  # F=1 manifold
    g2_idxs = [i for i, l in enumerate(cfg.atom_labels) if l[0] == "g2"]  # F=2 manifold
    e_idxs = [i for i, l in enumerate(cfg.atom_labels) if l[0].startswith("e")] # All excited states
    
    # Identify the specific level targeted by the probe (Ideal ground state: F=2, m=-2)
    ideal_g_probe = [i for i, l in enumerate(cfg.atom_labels) if l[1] == 2 and l[2] == -2][0]
    
    # Identify levels that cause "leaks" or noise in the system
    unwanted_g2_idxs = [i for i in g2_idxs if i != ideal_g_probe]
    e3_idxs = [i for i, l in enumerate(cfg.atom_labels) if l[1] == 3] # F'=3 (major source of noise)
    e2_idxs = [i for i, l in enumerate(cfg.atom_labels) if l[0] == "e2"] # F'=2 (target excited state)

    # --- 3. PROBE SCAN RANGE ---
    # Define the frequency range for the probe laser around the resonance center
    scan_width = 0.1 * cfg.gamma
    det_list = np.linspace(cfg.Delta_p_center - scan_width, cfg.Delta_p_center + scan_width, 300)

    # --- 4. BUILDING THE CONSTANT HAMILTONIAN (H_const) ---
    # This includes atomic energies, Zeeman shifts, the Coupling laser, and the Repumper
    H_const = 0
    for i, (label, f, m) in enumerate(cfg.atom_labels):
        # A. Atomic Energy: Set ground states to 0 and excited states relative to Delta_c
        if label.startswith("e"):
            energy = (getattr(cfg, f"E_{label}") - cfg.E_e2) - cfg.Delta_c
        else:
            energy = 0.0
        
        # B. Zeeman Shift: Calculate shift based on g-factor, B-field, and magnetic sub-level m
        g_f = -0.5 if "g1" in label else 0.5 if "g2" in label else 0.67
        zeeman_shift = g_f * cfg.mu_B * cfg.B_field * m
        
        # Add diagonal energy terms to the Hamiltonian
        H_const += (energy + zeeman_shift) * (basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag())

    # C. Coupling Laser (Interaction term)
    pol_c = -1  # sigma- polarization
    for gi in g2_idxs:
        for ei in e_idxs:
            c = cfg.safe_clebsch(cfg.atom_labels[gi][1], 1, cfg.atom_labels[ei][1], cfg.atom_labels[gi][2], pol_c, cfg.atom_labels[ei][2])
            if abs(c) > 1e-5:
                V = (cfg.Omega_c_amp * c / 2.0) * basis(cfg.N_atom, ei) * basis(cfg.N_atom, gi).dag()
                H_const += V + V.dag()

    # D. Repumper Laser (Interaction term: F=1 -> F'=2)
    pol_repump = -1
    Omega_repump = 0.5 * cfg.gamma
    for gi in g1_idxs:
        for ei in e_idxs:
            if cfg.atom_labels[ei][1] == 2: # Only target F'=2
                c = cfg.safe_clebsch(cfg.atom_labels[gi][1], 1, cfg.atom_labels[ei][1], cfg.atom_labels[gi][2], pol_repump, cfg.atom_labels[ei][2])
                if abs(c) > 1e-5:
                    V = (Omega_repump * c / 2.0) * basis(cfg.N_atom, ei) * basis(cfg.N_atom, gi).dag()
                    H_const += V + V.dag()

    # --- 5. COLLAPSE OPERATORS (Spontaneous Emission) ---
    # Defines the decay paths from excited states back to ground states via vacuum fluctuations
    c_ops = []
    for ei in e_idxs:
        for gi in g1_idxs + g2_idxs:
            for q in [-1, 0, 1]: # All polarizations (sigma-, pi, sigma+)
                c_vac = cfg.safe_clebsch(cfg.atom_labels[gi][1], 1, cfg.atom_labels[ei][1], cfg.atom_labels[gi][2], q, cfg.atom_labels[ei][2])
                if abs(c_vac) > 1e-5:
                    # Append decay operator to list for Lindblad solver
                    c_ops.append(np.sqrt(cfg.gamma) * c_vac * (basis(cfg.N_atom, gi) * basis(cfg.N_atom, ei).dag()))

    # --- 6. PROJECTION OPERATORS ---
    # These act as "virtual sensors" to measure population in specific states
    P_e_total = sum([basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag() for i in e_idxs])
    P_e2_only = sum([basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag() for i in e2_idxs])
    P_unwanted_g2 = sum([basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag() for i in unwanted_g2_idxs])
    P_g1_total = sum([basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag() for i in g1_idxs])
    P_e3_leakage = sum([basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag() for i in e3_idxs])

    # Dictionary to store numerical data during the scan
    results = {
        "det_list": det_list, "abs_tot": [], "abs_e2": [], 
        "pop_unw_g2": [], "pop_g1": [], "leak_e3": []
    }

    # --- 7. MAIN SIMULATION LOOP ---
    pol_p = 0 # Probe polarization (pi)
    for dp in det_list:
        H_probe = 0
        two_photon_det = dp - cfg.Delta_c # Relative detuning between probe and coupling
        
        # Add probe laser interaction and the rotating frame detuning
        for gi in g2_idxs:
            H_probe += two_photon_det * (basis(cfg.N_atom, gi) * basis(cfg.N_atom, gi).dag())
            for ei in e_idxs:
                c = cfg.safe_clebsch(cfg.atom_labels[gi][1], 1, cfg.atom_labels[ei][1], cfg.atom_labels[gi][2], pol_p, cfg.atom_labels[ei][2])
                if abs(c) > 1e-5:
                    V = (cfg.Omega_p_amp * c / 2.0) * basis(cfg.N_atom, ei) * basis(cfg.N_atom, gi).dag()
                    H_probe += V + V.dag()

        # Solve the Master Equation for the Steady State (rho_ss)
        rho_ss = steadystate(H_const + H_probe, c_ops)
        
        # Calculate expected values for each diagnostic operator
        results["abs_tot"].append(expect(P_e_total, rho_ss))
        results["abs_e2"].append(expect(P_e2_only, rho_ss))
        results["pop_unw_g2"].append(expect(P_unwanted_g2, rho_ss))
        results["pop_g1"].append(expect(P_g1_total, rho_ss))
        results["leak_e3"].append(expect(P_e3_leakage, rho_ss))

    # --- 8. SAVE ---
    for key in results:
        np.save(os.path.join(results_dir, f"{key}.npy"), np.array(results[key]))
    print(f"Success. Data saved in results_fano/{cfg.RUN_NAME}")

if __name__ == "__main__":
    run_fano()