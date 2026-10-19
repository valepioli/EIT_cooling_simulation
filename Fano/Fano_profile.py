import os
import numpy as np
import matplotlib.pyplot as plt
from qutip import *

# --- IMPOSTAZIONI GRAFICHE (Coerenti con gli altri script) ---
base_fs = 18
title_fs = 20
plt.rcParams.update({
    'font.size': base_fs,
    'axes.titlesize': title_fs,
    'axes.labelsize': base_fs,
    'xtick.labelsize': base_fs - 2,
    'ytick.labelsize': base_fs - 2,
    'legend.fontsize': base_fs - 2,
    'lines.linewidth': 2.5
})

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
excited_pop = []

for dp in detuning_probe_list:
    # Hamiltonian in rotating frame (Eq. 3 and 4 of the PDF)
    # H = -Delta_p|g1><g1| - Delta|g2><g2| + interaction
    # We set energy zero at excited state 'e'
    H = -dp * (g1 * g1.dag()) - Delta * (g2 * g2.dag()) + \
        (Omega_p/2.0) * (sigma_eg1 + sigma_eg1.dag()) + \
        (Omega_c/2.0) * (sigma_eg2 + sigma_eg2.dag())
    
    rho_ss = steadystate(H, c_ops)
    
    # Excited state population replaces the generic 'absorption'
    excited_pop.append(expect(e * e.dag(), rho_ss))

# 4. Plot for comparison with FIG. 2
fig, ax = plt.subplots(figsize=(12, 6))

ax.plot(detuning_probe_list, excited_pop, color='black', alpha=0.95, label='Excited State Population')
ax.axvline(Delta, color='red', linestyle='--', alpha=0.9, label="Two-Photon Resonance (Dark State)")

ax.set_title("Excited State Population Profile (3-Level System)")
ax.set_xlabel(r"Probe Detuning $\Delta_p$ [$\gamma$]")
ax.set_ylabel("Excited State Population [a.u.]")
ax.grid(True, alpha=0.3)
ax.legend(loc='best')

# --- SAVE SECTION ---
# Save image in current directory
filename = "fano_profile.png"
fig.tight_layout()
fig.savefig(filename, dpi=300, bbox_inches='tight')
print(f"Plot saved as {filename}")
# ---------------------

plt.show()