import numpy as np
from qutip import clebsch

RUN_NAME = "eit_axial" 

# ==========================================
# --- PHYSICAL CONSTANTS & SCALING ---
# ==========================================
gamma = 1.0                # Normalized decay rate
MHz = 1 / 6.067            # Scaling factor (1 MHz in units of Gamma)

E_e3 = +266.65 * MHz
E_e2 = 0.0
E_e1 = -156.94 * MHz
E_e0 = -229.16 * MHz

# ==========================================
# --- AXIAL SIMULATION PARAMETERS (Conveyor Belt) ---
# ==========================================
nu = 0.07 * gamma          # Axial trap frequency: 450 kHz
eta = 0.09                  # Axial Lamb-Dicke parameter (strong confinement)
n_phon = 2                  # Initial phonon population derived from Tz = 52.8 uK
N_vib = 10                  # Manageable Fock space ceiling for n=2

# ==========================================
# --- MAGNETIC FIELD PARAMETERS ---
# ==========================================
B_field = 4.0              # Gauss
mu_B = 1.399 * MHz         # Bohr magneton in MHz/Gauss
g_g1, g_g2 = -0.5, 0.5                 
m_g1, m_g2 = -1, -2                  

# ==========================================
# --- LASER PARAMETERS ---
# ==========================================
Delta_c = +13 * gamma 
zeeman_offset = (g_g2 * m_g2 - g_g1 * m_g1) * mu_B * B_field
Delta_p_center = Delta_c - zeeman_offset 

Omega_c_amp = np.sqrt(4 * np.abs(Delta_c) * nu)
Omega_p_amp = 0.09 * gamma
Omega_r_amp = 0.3 * gamma  
pol_r = -1                  

# ==========================================
# --- STATE MAPPING ---
# ==========================================
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