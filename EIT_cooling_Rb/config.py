import numpy as np
from qutip import clebsch

# --- FILE NAMING ---
RUN_NAME = "fano_profile_0" 

# --- Physical Constants & Scaling ---
gamma = 1.0  # Normalized decay rate
MHz = 1 / 6.067  # Scaling factor

# Hyperfine Energies for 87Rb (5P3/2) relative to F'=2
E_e3 = +266.65 * MHz
E_e2 = 0.0
E_e1 = -156.94 * MHz
E_e0 = -229.16 * MHz

# Trap frequency
nu = 0.016 * gamma

# --- Magnetic Field Parameters ---
B_field = 0.0             # Gauss
mu_B = 1.399 * MHz        # Bohr magneton
g_g1 = -0.5               
g_g2 = 0.5                
m_g1 = -1                 
m_g2 = -2                 

# Laser Parameters
Delta_c = +12.0 * gamma 
zeeman_offset = (g_g2 * m_g2 - g_g1 * m_g1) * mu_B * B_field
Delta_p_center = Delta_c - zeeman_offset 

Omega_c_amp = np.sqrt(4 * np.abs(Delta_c) * nu)
Omega_p_amp = 0.15 * gamma

# --- State Mapping (24 levels) ---
atom_labels = []
for f, m in [(1, m) for m in range(-1, 2)]: atom_labels.append(("g1", f, m))
for f, m in [(2, m) for m in range(-2, 3)]: atom_labels.append(("g2", f, m))
for f, m in [(0, m) for m in range(0, 1)]: atom_labels.append(("e0", f, m))
for f, m in [(1, m) for m in range(-1, 2)]: atom_labels.append(("e1", f, m))
for f, m in [(2, m) for m in range(-2, 3)]: atom_labels.append(("e2", f, m))
for f, m in [(3, m) for m in range(-3, 4)]: atom_labels.append(("e3", f, m))

N_atom = len(atom_labels)

def safe_clebsch(j1, j2, j3, m1, m2, m3):
    if not (abs(j1 - j2) <= j3 <= (j1 + j2)): return 0.0
    if m1 + m2 != m3: return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3: return 0.0
    return float(clebsch(j1, j2, j3, m1, m2, m3))