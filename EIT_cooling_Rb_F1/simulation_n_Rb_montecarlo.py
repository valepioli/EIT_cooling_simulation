import os
import numpy as np
from qutip import basis, tensor, qeye, destroy, fock, mcsolve, qload, qsave, expect, Options
from qutip.ui.progressbar import TextProgressBar
import config as cfg

# --- SIMULATION PARAMETERS ---
# SOFFITTO ALZATO A 15: Prevents wave function truncation (resolves NaN crashes)
N_vib = 15                  
eta = 0.25                 
dp = cfg.Delta_p_center    

# Monte Carlo time and trajectory settings
t_total = 500000.0         # Il tuo tempo richiesto
n_points = 5000            
n_traj = 5                 # ABBASSATO A 5 PER LA VELOCITÀ

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

    H_const_atom = 0 * basis(cfg.N_atom, 0) * basis(cfg.N_atom, 0).dag()
    
    for i, (label, f, m) in enumerate(cfg.atom_labels):
        if label.startswith("e"):
            energy = (getattr(cfg, f"E_{label}") - cfg.E_e2) - cfg.Delta_c
        else:
            energy = 0.0
        
        g_f = -0.5 if "g1" in label else 0.5 if "g2" in label else 0.67
        zeeman = g_f * cfg.mu_B * cfg.B_field * m
        H_const_atom += (energy + zeeman) * (basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag())

    pol_c = -1
    for gi in g2_idxs:
        for ei in e_idxs:
            c = cfg.safe_clebsch(cfg.atom_labels[gi][1], 1, cfg.atom_labels[ei][1], cfg.atom_labels[gi][2], pol_c, cfg.atom_labels[ei][2])
            if abs(c) > 1e-5:
                V = (cfg.Omega_c_amp * c / 2.0) * basis(cfg.N_atom, ei) * basis(cfg.N_atom, gi).dag()
                H_const_atom += V + V.dag()

    H_const = tensor(H_const_atom, qeye(N_vib))

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
                H_probe += tensor(V_atom, LD_term_plus) + tensor(V_atom.dag(), LD_term_minus)

    H_trap = cfg.nu * n_op
    H_tot = H_const + H_probe + H_trap

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
    
    # --- CREAZIONE DEI PROIETTORI PER OGNI STATO ATOMICO ---
    state_projectors = []
    state_labels = []
    for i, (label, f, m) in enumerate(cfg.atom_labels):
        # Proiettore per lo stato atomico i: |i><i|
        P_atom = basis(cfg.N_atom, i) * basis(cfg.N_atom, i).dag()
        # Espansione allo spazio tensoriale totale (Atomo x Fononi)
        P_tot = tensor(P_atom, qeye(N_vib))
        
        state_projectors.append(P_tot)
        state_labels.append(f"{label}_F{f}_m{m}")

    # Operatori da misurare: il numero di fononi + tutti i proiettori di stato
    e_ops = [n_op] + state_projectors
    
    ideal_g_probe = [i for i, l in enumerate(cfg.atom_labels) if l[1] == 2 and l[2] == -2][0]
    
    # QUI PARTE DA 5 FONONI ESATTI
    psi0 = tensor(basis(cfg.N_atom, ideal_g_probe), fock(N_vib, 5)) 
    
    t_list_full = np.linspace(0.0, t_total, n_points)
    
    # IMPOSTAZIONI SOLVER: Ottimizzate per velocità e stabilità senza crash ZVODE
    opts = Options(
        store_states=False, 
        nsteps=500000,    # Aumentato per coprire segmenti lunghi tra i salti
        atol=1e-7,        # Bilanciamento: abbastanza tollerante per la velocità, ma sicuro
        rtol=1e-5         # Bilanciamento: abbastanza tollerante per la velocità, ma sicuro
    )
    
    print(f"Target time: {t_total} | Calculating {n_traj} parallel trajectories...")
    
    result = mcsolve(H, psi0, t_list_full, c_ops, e_ops, 
                     ntraj=n_traj, 
                     options=opts, 
                     progress_bar=TextProgressBar())
    
    # Estrazione dei risultati
    n_expect_full = result.expect[0]
    
    populations = {}
    for idx, lbl in enumerate(state_labels):
        # L'indice è idx + 1 perché result.expect[0] è occupato da n_op
        populations[lbl] = result.expect[idx + 1]
    
    data_to_save = {
        't_list': t_list_full,
        'n_expect': n_expect_full,
        'populations': populations,  # <--- Le popolazioni salvate qui
        'params': {'N_vib': N_vib, 'eta': eta, 'dp': dp, 'ntraj': n_traj}
    }
    
    qsave(data_to_save, save_file)
    print(f"\n[+] Saved results to {save_file}.qu")
    print(f"--- SIMULATION COMPLETE. Final <n> = {n_expect_full[-1]:.3f} ---")

if __name__ == "__main__":
    run_monte_carlo_evolution()