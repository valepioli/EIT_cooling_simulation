import os
import numpy as np
from qutip import basis, tensor, qeye, destroy, fock, mesolve, qload, qsave, expect, Options
from qutip.ui.progressbar import TextProgressBar # Import for the text progress bar
import config as cfg

# --- SIMULATION PARAMETERS ---
N_vib = 15                 # Truncation of the Fock space
eta = 0.25                 # Lamb-Dicke parameter
dp = cfg.Delta_p_center    # Probe detuning (centered on EIT resonance)

# Time chunking settings
t_total = 100000.0           # Total time you want to reach
t_chunk_size = 50.0        # How many time units per chunk
steps_per_chunk = 50       # Resolution inside each chunk

# --- PATH SETUP ---
script_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(script_dir, "results_time")
os.makedirs(results_dir, exist_ok=True)
save_file = os.path.join(results_dir, f"time_evol_{cfg.RUN_NAME}")

def build_system():
    print("Building 600x600 Tensor Hamiltonian (24 Atom levels x 15 Fock states)...")
    
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
    Omega_repump =  cfg.Omega_repump
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

def run_time_evolution():
    H, c_ops, n_op = build_system()
    
    # 1. State or Load Checkpoint
    if os.path.exists(save_file + ".qu"):
        print(f"\n[+] Checkpoint found! Loading data from {save_file}.qu...")
        data = qload(save_file)
        t_list_full = list(data['t_list'])
        n_expect_full = list(data['n_expect'])
        rho_current = data['last_rho']
        current_time = t_list_full[-1]
    else:
        print("\n[!] No checkpoint found. Starting fresh simulation...")
        ideal_g_probe = [i for i, l in enumerate(cfg.atom_labels) if l[1] == 2 and l[2] == -2][0]
        psi0 = tensor(basis(cfg.N_atom, ideal_g_probe), fock(N_vib, 10)) # Starting at n=10
        rho_current = psi0 * psi0.dag()
        
        t_list_full = [0.0]
        n_expect_full = [expect(n_op, rho_current)]
        current_time = 0.0

    print(f"Target time: {t_total} | Current time: {current_time:.1f}")
    
    opts = Options(nsteps=20000, store_states=True)

    # 2. Chunking Loop
    while current_time < t_total:
        end_time = min(current_time + t_chunk_size, t_total)
        t_chunk = np.linspace(current_time, end_time, steps_per_chunk)
        
        print(f"\n--> Simulating chunk: t = [{current_time:.1f} to {end_time:.1f}]")
        
        # AGGIUNTO: progress_bar=TextProgressBar() per visualizzare lo stato di avanzamento
        result = mesolve(H, rho_current, t_chunk, c_ops, [n_op], options=opts, progress_bar=TextProgressBar())
        
        # Extract data (ignoring index 0 to prevent duplicating the overlapping boundary point)
        n_chunk = result.expect[0]
        t_list_full.extend(t_chunk[1:])
        n_expect_full.extend(n_chunk[1:])
        
        # Update current state to the last frame of the chunk
        rho_current = result.states[-1]
        current_time = end_time
        
        # 3. Save Data Incrementally
        data_to_save = {
            't_list': np.array(t_list_full),
            'n_expect': np.array(n_expect_full),
            'last_rho': rho_current,
            'params': {'N_vib': N_vib, 'eta': eta, 'dp': dp}
        }
        qsave(data_to_save, save_file)
        print(f"    Saved checkpoint. Current <n> = {n_chunk[-1]:.3f}")

    print("\n--- SIMULATION COMPLETE ---")

if __name__ == "__main__":
    run_time_evolution()