import numpy as np
from qutip import basis, tensor, fock, qeye

# --- PHYSICAL PARAMETERS ---
gamma = 1.0         
nu = 0.5 * gamma    # Trap frequency
Delta_c = 15.0 * gamma 
Delta_p = Delta_c   

# Stark shift condition for EIT cooling
Omega_c = np.sqrt(4 * np.abs(Delta_c) * nu) 
eta = 0.35          # Lamb-Dicke parameter
Omega_p = 0.1 * gamma 
N_vib = 25          

# --- SIMULATION PARAMETERS ---
t_stop = 3500
t_points = 200
t_list = np.linspace(0, t_stop, t_points)

# Initial state: atom in g1, oscillator in Fock state n=15
psi0 = tensor(basis(3,0), fock(N_vib, 15))

# --- PATH ---
SAVE_PATH = "results/eit_simulation_results.qu"
