import os
import numpy as np
from qutip import basis, tensor, qeye, destroy, fock, mcsolve, qload, qsave, expect, Options
from qutip.ui.progressbar import TextProgressBar
import config as cfg

# --- SIMULATION PARAMETERS ---
N_vib = 15                 # Truncation of the Fock space (Must be > 10 if starting at n=10)
eta = 0.25                 # Lamb-Dicke parameter
dp = cfg.Delta_p_center    # Probe detuning (centered on EIT resonance)

# Monte Carlo time and trajectory settings
t_total = 500000.0         # Total simulation time
n_points = 2000            # Number of data points to save across the total time
n_traj = 100               # Number of quantum trajectories to average over (increase for smoother curves)

# --- PATH SETUP ---
script_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(script_dir, "results_time")
os.makedirs(results_dir, exist_ok=True)
# Changed the save file name slightly to distinguish it from the mesolve runs
save_file = os.path.join(results_dir, f"time_evol_500k_MC_{cfg.RUN_NAME}")

def build_system():
    print(f"Building Tensor Hamiltonian (Atom levels x {N_vib} Fock states)...")
    
    # 1. Operators & Indices
    a = tensor(qeye(cfg.N_atom), destroy(N_vib))
    n_op = a.dag() * a
    
    g1_idxs = [i for i, l in enumerate(cfg.atom_labels) if l[0] == "g1"]
    g2_idxs = [i for i, l in enumerate(cfg.atom_labels) if l[0] == "g2"]
    e_idxs = [i for i, l in enumerate(cfg.atom_labels) if l[0].startswith("e")]

    # 2. Atomic Constant Hamiltonian (Energies, Zeeman, Coupling, Repump)
    H_const_atom = 0 * basis(cfg.N_atom, 0) * basis(cfg.N_atom, 0).dag()
    
    for i, (label, f, m) in enumerate(cfg.atom_labels):
        if label.startswith("e"):
            energy = (getattr(cfg, f"E_{label}") - cfg.E_e2) - cfg.Delta_c
        else:
            energy = 0.0
        
        g_f = -0.5 if "g1" in label else 0.5 if "g2" in label else 0.67
        zeeman = g_f * cfg.mu_B * cfg.B_field * m
        H_const_atom += (energy + zeeman) * (basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag())

    # Coupling Laser
    pol_c = -1
    for gi in g2_idxs:
        for ei in e_idxs:
            c = cfg.safe_clebsch(cfg.atom_labels[gi][1], 1, cfg.atom_labels[ei][1], cfg.atom_labels[gi][2], pol_c, cfg.atom_labels[ei][2])
            if abs(c) > 1e-5:
                V = (cfg.Omega_c_amp * c / 2.0) * basis(cfg.N_atom, ei) * basis(cfg.N_atom, gi).dag()
                H_const_atom += V + V.dag()

    # Repumper Laser
    pol_repump = -1
    Omega_repump = 0.5 * cfg.gamma
    for gi in g1_idxs:
        for ei in e_idxs:
            if cfg.atom_labels[ei][1] == 2:
                c = cfg.safe_clebsch(cfg.atom_labels[gi][1], 1, cfg.atom_labels[ei][1], cfg.atom_labels[gi][2], pol_repump, cfg.atom_labels[ei][2])
                if abs(c) > 1e-5:
                    V = (Omega_repump * c / 2.0) * basis(cfg.N_atom, ei) * basis(cfg.N_atom, gi).dag()
                    H_const_atom += V + V.dag()

    # Tensor H_const with identity in vibration space
    H_const = tensor(H_const_atom, qeye(N_vib))

    # 3. Probe Laser (with Lamb-Dicke coupling to phonons)
    H_probe = 0
    two_photon_det = dp - cfg.Delta_c
    
    H_probe_det_atom = sum([two_photon_det * (basis(cfg.N_atom, gi) * basis(cfg.N_atom, gi).dag()) for gi in g2_idxs])
    H_probe += tensor(H_probe_det_atom, qeye(N_vib))

    pol_p = 0
    LD_term_plus = qeye(N_vib) + 1j * eta * (destroy(N_vib) + destroy(N_vib).dag())
    LD_term_minus = qeye(N_vib) - 1j * eta * (destroy(N_vib) + destroy(N_vib).dag())

    for gi in g2_idxs:
        for ei in e_idxs:
            c = cfg.safe_clebsch(cfg.atom_labels[gi][1], 1, cfg.atom_labels[ei][1], cfg.atom_labels[gi][2], pol_p, cfg.atom_labels[ei][2])
            if abs(c) > 1e-5:
                V_atom = (cfg.Omega_p_amp * c / 2.0) * basis(cfg.N_atom, ei) * basis(cfg.N_atom, gi).dag()
                # Apply Lamb-Dicke spatial phase
                H_probe += tensor(V_atom, LD_term_plus) + tensor(V_atom.dag(), LD_term_minus)

    # 4. Trap frequency Hamiltonian
    H_trap = cfg.nu * n_op
    
    # Final Hamiltonian
    H_tot = H_const + H_probe + H_trap

    # 5. Collapse Operators (Spontaneous Emission)
    c_ops = []
    for ei in e_idxs:
        for gi in g1_idxs + g2_idxs:
            for q in [-1, 0, 1]:
                c_vac = cfg.safe_clebsch(cfg.atom_labels[gi][1], 1, cfg.atom_labels[ei][1], cfg.atom_labels[gi][2], q, cfg.atom_labels[ei][2])
                if abs(c_vac) > 1e-5:
                    decay_op = np.sqrt(cfg.gamma) * c_vac * (basis(cfg.N_atom, gi) * basis(cfg.N_atom, ei).dag())
                    c_ops.append(tensor(decay_op, qeye(N_vib)))

    return H_tot, c_ops, n_op


def run_monte_carlo_evolution():
    H, c_ops, n_op = build_system()
    
    print("\n[+] System built successfully. Initializing Monte Carlo parameters...")
    
    # 1. State Definition
    # Note: mcsolve requires a state vector (ket), NOT a density matrix (rho)
    ideal_g_probe = [i for i, l in enumerate(cfg.atom_labels) if l[1] == 2 and l[2] == -2][0]
    
    # Ensure N_vib is strictly greater than the starting phonon number to avoid ValueError
    psi0 = tensor(basis(cfg.N_atom, ideal_g_probe), fock(N_vib, 10)) 
    
    # 2. Time array creation (No chunking needed for mcsolve)
    t_list_full = np.linspace(0.0, t_total, n_points)
    
    # 3. Solver options
    # store_states=False prevents RAM overload by only keeping the expectation values
    opts = Options(store_states=False)
    
    print(f"Target time: {t_total} | Calculating {n_traj} parallel trajectories...")
    
    # 4. Execute Monte Carlo solver
    # This will automatically utilize all available CPU cores
    result = mcsolve(H, psi0, t_list_full, c_ops, [n_op], 
                     ntraj=n_traj, 
                     options=opts, 
                     progress_bar=TextProgressBar())
    
    # 5. Extract data (result.expect[0] corresponds to the expectation value of n_op)
    n_expect_full = result.expect[0]
    
    # 6. Save data to file
    data_to_save = {
        't_list': t_list_full,
        'n_expect': n_expect_full,
        'params': {'N_vib': N_vib, 'eta': eta, 'dp': dp, 'ntraj': n_traj}
    }
    qsave(data_to_save, save_file)
    print(f"\n[+] Saved results to {save_file}.qu")
    print(f"--- SIMULATION COMPLETE. Final <n> = {n_expect_full[-1]:.3f} ---")

if __name__ == "__main__":
    run_monte_carlo_evolution()