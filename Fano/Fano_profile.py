import numpy as np
import matplotlib.pyplot as plt
from qutip import *

# Parameters taken from FIG. 2 caption of the PDF:
# Delta = 2.5 * gamma
# Omega_control = 1.0 * gamma
# Omega_probe = 0.05 * gamma

gamma = 1.0                 # Unit of measurement
Delta = 2.5 * gamma         # Control laser detuning
Omega_c = 1.0 * gamma       # Control Rabi frequency
Omega_p = 0.05 * gamma      # Probe Rabi frequency (very weak)
gamma_deph = 0.001 * gamma  # Small dephasing for numerical stability

# 1. Basis and operators definition
g1 = basis(3, 0)  # Ground state for probe
g2 = basis(3, 1)  # Ground state for control
e  = basis(3, 2)  # Excited state

sigma_eg1 = e * g1.dag()
sigma_eg2 = e * g2.dag()

# 2. Dissipation (Kρ in the PDF)
# Spontaneous decay into the two ground states
c_ops = [np.sqrt(gamma/2) * g1 * e.dag(), 
         np.sqrt(gamma/2) * g2 * e.dag(),
         np.sqrt(gamma_deph) * (g1*g1.dag() - g2*g2.dag())]

# 3. Probe detuning scan (Delta_P)
# Scan around atomic resonance
detuning_probe_list = np.linspace(-2 * gamma, 5 * gamma, 500)
absorption = []

for dp in detuning_probe_list:
    # Hamiltonian in rotating frame (Eq. 3 and 4 of the PDF)
    # H = -Delta_p|g1><g1| - Delta|g2><g2| + interaction
    # We set energy zero at excited state 'e'
    H = -dp * (g1 * g1.dag()) - Delta * (g2 * g2.dag()) + \
        (Omega_p/2.0) * (sigma_eg1 + sigma_eg1.dag()) + \
        (Omega_c/2.0) * (sigma_eg2 + sigma_eg2.dag())
    
    rho_ss = steadystate(H, c_ops)
    
    # Absorption I(Delta_P) is proportional to excited state population
    absorption.append(expect(e * e.dag(), rho_ss))

# 4. Plot for comparison with FIG. 2
plt.figure(figsize=(8, 5))
plt.plot(detuning_probe_list, absorption, color='black', lw=1.5)
plt.axvline(Delta, color='red', linestyle='--', label="Two-photon resonance (Dark state)")
plt.title("Excitation Spectrum")
plt.xlabel("Probe Detuning $\Delta_P / \gamma$")
plt.ylabel("Excitation $I(\Delta_P)$")
plt.grid(True, alpha=0.2)
plt.legend()

# --- SAVE SECTION ---
# Save image in current directory
filename = "fano_profile.png"
plt.savefig(filename, dpi=300, bbox_inches='tight')
print(f"Plot saved as {filename}")
# ---------------------

plt.show()