import os
import numpy as np
from qutip import *
import config

# --- ABSOLUTE PATH MANAGEMENT ---
# Finds the folder where this file (simulation.py) is located
# If the file is in /root/EIT_cooling_simulation/EIT_cooling/simulation.py,
# script_dir will be /root/EIT_cooling_simulation/EIT_cooling/
script_dir = os.path.dirname(os.path.abspath(__file__))
results_path = os.path.join(script_dir, 'results')

# Create 'results' folder if it does not exist
if not os.path.exists(results_path):
    os.makedirs(results_path)
    print(f"Folder created: {results_path}")

def run_simulation():
    print(f"--- STARTING SIMULATION ---")
    print(f"Delta_c: {config.Delta_c}, Omega_c: {config.Omega_c:.3f}, nu: {config.nu}")

    # 1. Operators
    # Define Hilbert space: (3-level atom) x (oscillator N_vib)
    a = tensor(qeye(3), destroy(config.N_vib))
    n_op = a.dag() * a
    
    # Atomic operators: 0=|g1>, 1=|g2>, 2=|e>
    sig_eg1 = tensor(basis(3,2) * basis(3,0).dag(), qeye(config.N_vib))
    sig_eg2 = tensor(basis(3,2) * basis(3,1).dag(), qeye(config.N_vib))
    proj_g1 = tensor(basis(3,0) * basis(3,0).dag(), qeye(config.N_vib))
    proj_g2 = tensor(basis(3,1) * basis(3,1).dag(), qeye(config.N_vib))

    # 2. Hamiltonian
    H = (config.Delta_p * proj_g1 + config.Delta_c * proj_g2 + config.nu * n_op +
         (config.Omega_c/2.0) * (sig_eg2 + sig_eg2.dag()) +
         (config.Omega_p/2.0) * (sig_eg1 * (1 + 1j * config.eta * (a + a.dag())) + 
                          sig_eg1.dag() * (1 - 1j * config.eta * (a + a.dag()))))

    # Collapse operators
    c_ops = [np.sqrt(config.gamma/2) * sig_eg1.dag(), 
             np.sqrt(config.gamma/2) * sig_eg2.dag()]

    # 3. Master equation solution
    print("Solving dynamics (mesolve)...")
    
    # Using Options object (required in QuTiP 4.x)
    opts = Options(nsteps=10000, store_states=True)
    
    result = mesolve(H, config.psi0, config.t_list, c_ops, [n_op], 
                     options=opts, 
                     progress_bar=True)

    # 4. EIT spectrum (steady-state scan)
    print("Computing EIT spectrum...")
    det_scan = np.linspace(config.Delta_c - 1.5, config.Delta_c + 1.5, 400)
    spec = []
    for d in det_scan:
        H_at = +d*basis(3,0)*basis(3,0).dag() + config.Delta_c*basis(3,1)*basis(3,1).dag() + \
               (config.Omega_c/2)*(basis(3,2)*basis(3,1).dag() + basis(3,1)*basis(3,2).dag()) + \
               (config.Omega_p/2.0)*(basis(3,2)*basis(3,0).dag() + basis(3,0)*basis(3,2).dag())
        
        c_ops_at = [np.sqrt(config.gamma/2)*basis(3,0)*basis(3,2).dag(), 
                    np.sqrt(config.gamma/2)*basis(3,1)*basis(3,2).dag()]
        
        rho_s = steadystate(H_at, c_ops_at)
        spec.append(expect(basis(3,2)*basis(3,2).dag(), rho_s))

    # 5. Data saving (robust path handling)
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
    
    # Extracts only filename from config path
    filename = os.path.basename(config.SAVE_PATH).replace('.qu', '')
    
    # Builds full absolute path
    final_path = os.path.join(results_path, filename)
    
    qsave(data_to_save, final_path)
    print(f"\n--- SIMULATION COMPLETED ---")
    print(f"Data successfully saved to: {final_path}.qu")

if __name__ == "__main__":
    run_simulation()
