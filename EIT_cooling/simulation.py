import os
import numpy as np
from qutip import *
import config

# Creazione cartella risultati se non esiste
if not os.path.exists('results'):
    os.makedirs('results')

def run_simulation():
    print(f"--- AVVIO SIMULAZIONE ---")
    print(f"Delta_c: {config.Delta_c}, Omega_c: {config.Omega_c:.3f}, nu: {config.nu}")

    # 1. Operatori
    a = tensor(qeye(3), destroy(config.N_vib))
    n_op = a.dag() * a
    sig_eg1 = tensor(basis(3,2) * basis(3,0).dag(), qeye(config.N_vib))
    sig_eg2 = tensor(basis(3,2) * basis(3,1).dag(), qeye(config.N_vib))
    proj_g1 = tensor(basis(3,0) * basis(3,0).dag(), qeye(config.N_vib))
    proj_g2 = tensor(basis(3,1) * basis(3,1).dag(), qeye(config.N_vib))

    # 2. Hamiltoniana
    H = (-config.Delta_p * proj_g1 - config.Delta_c * proj_g2 + config.nu * n_op +
         (config.Omega_c/2.0) * (sig_eg2 + sig_eg2.dag()) +
         (config.Omega_p/2.0) * (sig_eg1 * (1 + 1j * config.eta * (a + a.dag())) + 
                          sig_eg1.dag() * (1 - 1j * config.eta * (a + a.dag()))))

    c_ops = [np.sqrt(config.gamma/2) * sig_eg1.dag(), np.sqrt(config.gamma/2) * sig_eg2.dag()]

    # 3. Risoluzione Master Equation
    print("Risoluzione dinamica (mesolve)...")
    result = mesolve(H, config.psi0, config.t_list, c_ops, [n_op], 
                     options={'nsteps': 10000, 'store_states': True}, 
                     progress_bar=True)

    # 4. Spettro EIT (Steady State Scan)
    print("Calcolo spettro EIT...")
    det_scan = np.linspace(config.Delta_c - 1.5, config.Delta_c + 0.5, 400)
    spec = []
    for d in det_scan:
        H_at = -d*basis(3,0)*basis(3,0).dag() - config.Delta_c*basis(3,1)*basis(3,1).dag() + \
               (config.Omega_c/2)*(basis(3,2)*basis(3,1).dag() + basis(3,1)*basis(3,2).dag()) + \
               (0.1/2)*(basis(3,2)*basis(3,0).dag() + basis(3,0)*basis(3,2).dag())
        rho_s = steadystate(H_at, [np.sqrt(config.gamma/2)*basis(3,0)*basis(3,2).dag(), 
                                    np.sqrt(config.gamma/2)*basis(3,1)*basis(3,2).dag()])
        spec.append(expect(basis(3,2)*basis(3,2).dag(), rho_s))

    # 5. Salvataggio
    data_to_save = {
        't_list': config.t_list,
        'expect_n': result.expect[0],
        'states': result.states,
        'det_scan': det_scan,
        'spec': spec,
        # Salviamo i parametri correnti per riferimento nel plot
        'params': {
            'Delta_p': config.Delta_p,
            'nu': config.nu,
            'N_vib': config.N_vib
        }
    }
    
    qsave(data_to_save, config.SAVE_PATH.replace('.qu', ''))
    print(f"Dati salvati in: {config.SAVE_PATH}")

if __name__ == "__main__":
    run_simulation()