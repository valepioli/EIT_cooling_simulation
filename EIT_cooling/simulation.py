import os
import numpy as np
from qutip import *
import config

# --- GESTIONE PERCORSI ASSOLUTI ---
# Trova la cartella dove si trova questo file (simulation.py)
# Se il file è in /root/EIT_cooling_simulation/EIT_cooling/simulation.py, 
# script_dir sarà /root/EIT_cooling_simulation/EIT_cooling/
script_dir = os.path.dirname(os.path.abspath(__file__))
results_path = os.path.join(script_dir, 'results')

# Crea la cartella 'results' se non esiste
if not os.path.exists(results_path):
    os.makedirs(results_path)
    print(f"Cartella creata: {results_path}")

def run_simulation():
    print(f"--- AVVIO SIMULAZIONE ---")
    print(f"Delta_c: {config.Delta_c}, Omega_c: {config.Omega_c:.3f}, nu: {config.nu}")

    # 1. Operatori
    # Definiamo lo spazio di Hilbert: (Atomo 3 livelli) x (Oscillatore N_vib)
    a = tensor(qeye(3), destroy(config.N_vib))
    n_op = a.dag() * a
    
    # Operatori atomici: 0=|g1>, 1=|g2>, 2=|e>
    sig_eg1 = tensor(basis(3,2) * basis(3,0).dag(), qeye(config.N_vib))
    sig_eg2 = tensor(basis(3,2) * basis(3,1).dag(), qeye(config.N_vib))
    proj_g1 = tensor(basis(3,0) * basis(3,0).dag(), qeye(config.N_vib))
    proj_g2 = tensor(basis(3,1) * basis(3,1).dag(), qeye(config.N_vib))

    # 2. Hamiltoniana
    H = (-config.Delta_p * proj_g1 - config.Delta_c * proj_g2 + config.nu * n_op +
         (config.Omega_c/2.0) * (sig_eg2 + sig_eg2.dag()) +
         (config.Omega_p/2.0) * (sig_eg1 * (1 + 1j * config.eta * (a + a.dag())) + 
                          sig_eg1.dag() * (1 - 1j * config.eta * (a + a.dag()))))

    # Operatori di collasso
    c_ops = [np.sqrt(config.gamma/2) * sig_eg1.dag(), 
             np.sqrt(config.gamma/2) * sig_eg2.dag()]

    # 3. Risoluzione Master Equation
    print("Risoluzione dinamica (mesolve)...")
    
    # Utilizzo dell'oggetto Options (necessario per QuTiP 4.x)
    opts = Options(nsteps=10000, store_states=True)
    
    result = mesolve(H, config.psi0, config.t_list, c_ops, [n_op], 
                     options=opts, 
                     progress_bar=True)

    # 4. Spettro EIT (Steady State Scan)
    print("Calcolo spettro EIT...")
    det_scan = np.linspace(config.Delta_c - 1.5, config.Delta_c + 0.5, 400)
    spec = []
    for d in det_scan:
        H_at = -d*basis(3,0)*basis(3,0).dag() - config.Delta_c*basis(3,1)*basis(3,1).dag() + \
               (config.Omega_c/2)*(basis(3,2)*basis(3,1).dag() + basis(3,1)*basis(3,2).dag()) + \
               (0.1/2)*(basis(3,2)*basis(3,0).dag() + basis(3,0)*basis(3,2).dag())
        
        c_ops_at = [np.sqrt(config.gamma/2)*basis(3,0)*basis(3,2).dag(), 
                    np.sqrt(config.gamma/2)*basis(3,1)*basis(3,2).dag()]
        
        rho_s = steadystate(H_at, c_ops_at)
        spec.append(expect(basis(3,2)*basis(3,2).dag(), rho_s))

    # 5. Salvataggio Dati (Gestione Percorso Robusta)
    data_to_save = {
        't_list': config.t_list,
        'expect_n': result.expect[0],
        'states': result.states,
        'det_scan': det_scan,
        'spec': spec,
        'params': {
            'Delta_p': config.Delta_p,
            'Delta_c': config.Delta_c,
            'nu': config.nu,
            'N_vib': config.N_vib,
            'Omega_c': config.Omega_c
        }
    }
    
    # Prende solo il nome del file dal config (es: "eit_simulation_results.qu" -> "eit_simulation_results")
    filename = os.path.basename(config.SAVE_PATH).replace('.qu', '')
    
    # Crea il percorso finale assoluto: /root/.../EIT_cooling/results/eit_simulation_results
    final_path = os.path.join(results_path, filename)
    
    qsave(data_to_save, final_path)
    print(f"\n--- SIMULAZIONE COMPLETATA ---")
    print(f"Dati salvati con successo in: {final_path}.qu")

if __name__ == "__main__":
    run_simulation()