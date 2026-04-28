import numpy as np
import matplotlib.pyplot as plt
from qutip import *

# Parametri presi dalla didascalia della FIG. 2 del PDF:
# Delta = 2.5 * gamma
# Omega_control = 1.0 * gamma
# Omega_probe = 0.05 * gamma

gamma = 1.0             # Unità di misura
Delta = 2.5 * gamma     # Detuning del laser di controllo
Omega_c = 1.0 * gamma   # Rabi del controllo
Omega_p = 0.05 * gamma  # Rabi del probe (molto debole)
gamma_deph = 0.001 * gamma # Piccola decoerenza per stabilità numerica

# 1. Definizione Base e Operatori
g1 = basis(3, 0) # Ground per il probe
g2 = basis(3, 1) # Ground per il controllo
e  = basis(3, 2) # Stato eccitato

sigma_eg1 = e * g1.dag()
sigma_eg2 = e * g2.dag()

# 2. Dissipazione (Kρ nel PDF)
# Decadimento spontaneo verso i due stati di ground
c_ops = [np.sqrt(gamma/2) * g1 * e.dag(), 
         np.sqrt(gamma/2) * g2 * e.dag(),
         np.sqrt(gamma_deph) * (g1*g1.dag() - g2*g2.dag())]

# 3. Scan del detuning del probe (Delta_P)
# Scansioniamo intorno alla risonanza atomica
detuning_probe_list = np.linspace(-2 * gamma, 5 * gamma, 500)
absorption = []

for dp in detuning_probe_list:
    # Hamiltoniano in rotating frame (Eq. 3 e 4 del PDF)
    # H = -Delta_p|g1><g1| - Delta|g2><g2| + Interazione
    # Usiamo lo zero di energia sullo stato eccitato 'e'
    H = -dp * (g1 * g1.dag()) - Delta * (g2 * g2.dag()) + \
        (Omega_p/2.0) * (sigma_eg1 + sigma_eg1.dag()) + \
        (Omega_c/2.0) * (sigma_eg2 + sigma_eg2.dag())
    
    rho_ss = steadystate(H, c_ops)
    
    # L'assorbimento I(Delta_P) è proporzionale alla popolazione di 'e' 
    # o alla parte immaginaria della coerenza del probe.
    absorption.append(expect(e * e.dag(), rho_ss))

# 4. Plot per confronto con FIG. 2
plt.figure(figsize=(8, 5))
plt.plot(detuning_probe_list, absorption, color='black', lw=1.5)
plt.axvline(Delta, color='red', linestyle='--', label="Risonanza a 2 fotoni (Dark)")
plt.title("Replica FIG. 2 - Morigi Paper (Excitation Spectrum)")
plt.xlabel("Detuning Probe $\Delta_P / \gamma$")
plt.ylabel("Excitation $I(\Delta_P)$")
plt.grid(True, alpha=0.2)
plt.legend()

# --- MODIFICA QUI PER SALVARE ---
# Salva l'immagine nella cartella corrente
filename = "fano_profile.png"
plt.savefig(filename, dpi=300, bbox_inches='tight')
print(f"Grafico salvato come {filename}")
# --------------------------------

plt.show()