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

    # --- 4. BUILDING THE CONSTANT HAMILTONIAN (H_const) ---
    H_const = 0
    # Add energies and Zeeman shifts for all states
    for i, (label, f, m) in enumerate(cfg.atom_labels):
        if label.startswith("e"):
            energy = (getattr(cfg, f"E_{label}") - cfg.E_e2) - cfg.Delta_c
        else:
            energy = 0.0
        
        g_f = -0.5 if "g1" in label else 0.5 if "g2" in label else 0.67
        zeeman_shift = g_f * cfg.mu_B * cfg.B_field * m
        H_const += (energy + zeeman_shift) * (basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag())

    # --- UPDATED: Control Laser su F=2 ---
    pol_c = 1  # sigma+
    for gi in g2_idxs:
        for ei in e_idxs:
            c = cfg.safe_clebsch(cfg.atom_labels[gi][1], 1, cfg.atom_labels[ei][1], cfg.atom_labels[gi][2], pol_c, cfg.atom_labels[ei][2])
            if abs(c) > 1e-5:
                V = (cfg.Omega_c_amp * c / 2.0) * basis(cfg.N_atom, ei) * basis(cfg.N_atom, gi).dag()
                H_const += V + V.dag()

    # --- UPDATED: Resonant Repumper Laser ---
    # Drives F=2, m=2 to F'=2, m=1 with sigma- (pol=-1)
    try:
        g_rep_idx = cfg.atom_labels.index(("g2", 2, 2))
        e_rep_idx = cfg.atom_labels.index(("e2", 2, 1))
        
        c_rep = cfg.safe_clebsch(2, 1, 2, 2, cfg.pol_r, 1)
        
        if abs(c_rep) > 1e-5:
            V_rep = (cfg.Omega_r_amp * c_rep / 2.0) * basis(cfg.N_atom, e_rep_idx) * basis(cfg.N_atom, g_rep_idx).dag()
            H_const += V_rep + V_rep.dag()
            
    except ValueError:
        print("Error: Repumper states not found in atom_labels.")

    # --- 5. COLLAPSE OPERATORS (Spontaneous Emission) ---
    c_ops = []
    for ei in e_idxs:
        for gi in g1_idxs + g2_idxs:
            for q in [-1, 0, 1]: 
                c_vac = cfg.safe_clebsch(cfg.atom_labels[gi][1], 1, cfg.atom_labels[ei][1], cfg.atom_labels[gi][2], q, cfg.atom_labels[ei][2])
                if abs(c_vac) > 1e-5:
                    c_ops.append(np.sqrt(cfg.gamma) * c_vac * (basis(cfg.N_atom, gi) * basis(cfg.N_atom, ei).dag()))

    # --- 6. PROJECTION OPERATORS & DICTIONARY SETUP ---
    P_e_total = sum([basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag() for i in e_idxs])
    P_e2_only = sum([basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag() for i in e2_idxs])
    P_e3_leakage = sum([basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag() for i in e3_idxs])

    results = {
        "det_list": det_list, "abs_tot": [], "abs_e2": [], "leak_e3": []
    }
    
    for label, f, m in cfg.atom_labels:
        results[f"pop_{label}_m{m}"] = []

    # --- 7. MAIN SIMULATION LOOP ---
    pol_p = 1 # sigma+ polarization per il probe
    for dp in det_list:
        H_probe = 0
        two_photon_det = dp - cfg.Delta_c 
        
        # --- UPDATED: Probe scan e detuning applicato a F=1 ---
        for gi in g1_idxs:
            H_probe += two_photon_det * (basis(cfg.N_atom, gi) * basis(cfg.N_atom, gi).dag())
            for ei in e_idxs:
                c = cfg.safe_clebsch(cfg.atom_labels[gi][1], 1, cfg.atom_labels[ei][1], cfg.atom_labels[gi][2], pol_p, cfg.atom_labels[ei][2])
                if abs(c) > 1e-5:
                    V = (cfg.Omega_p_amp * c / 2.0) * basis(cfg.N_atom, ei) * basis(cfg.N_atom, gi).dag()
                    H_probe += V + V.dag()

        # Calculate steady state per quest detuning
        rho_ss = steadystate(H_const + H_probe, c_ops)
        
        # Extract global metrics
        results["abs_tot"].append(expect(P_e_total, rho_ss))
        results["abs_e2"].append(expect(P_e2_only, rho_ss))
        results["leak_e3"].append(expect(P_e3_leakage, rho_ss))
        
        # Extract individual populations
        for i, (label, f, m) in enumerate(cfg.atom_labels):
            P_op = basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag()
            results[f"pop_{label}_m{m}"].append(expect(P_op, rho_ss))

    # --- 8. SAVE RESULTS ---
    for key in results:
        np.save(os.path.join(results_dir, f"{key}.npy"), np.array(results[key]))
        
    print(f"Success. All data saved in results_fano/{cfg.RUN_NAME}")

if __name__ == "__main__":
    run_fano()