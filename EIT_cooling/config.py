import numpy as np
from qutip import basis, tensor, fock, qeye

# --- PARAMETRI FISICI ---
gamma = 1.0         
nu = 0.5 * gamma    # Frequenza trappola
Delta_c = -15.0 * gamma 
Delta_p = Delta_c   

# Condizione Stark shift per EIT cooling
Omega_c = np.sqrt(4 * np.abs(Delta_c) * nu) 
eta = 0.35          # Lamb-Dicke
Omega_p = 0.3 * gamma 
N_vib = 25          

# --- PARAMETRI SIMULAZIONE ---
t_stop = 1500
t_points = 200
t_list = np.linspace(0, t_stop, t_points)

# Stato iniziale: atomo in g1, oscillatore in stato di Fock n=15
psi0 = tensor(basis(3,0), fock(N_vib, 15))

# --- PATH ---
SAVE_PATH = "results/eit_simulation_results.qu"