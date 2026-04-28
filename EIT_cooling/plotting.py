import os
from qutip import qload
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# --- 1. CONFIGURAZIONE PERCORSI ---
# Trova la cartella dove si trova questo script di plotting
script_dir = os.path.dirname(os.path.abspath(__file__))
# Punta alla cartella results (che è nello stesso posto di simulation.py)
results_path = os.path.join(script_dir, 'results')

# Nome del file (senza estensione .qu, qload la gestisce da solo)
filename = 'eit_simulation_results'
full_load_path = os.path.join(results_path, filename)

# --- 2. CARICAMENTO DATI ---
print(f"Cerco i dati in: {full_load_path}.qu")

if os.path.exists(f"{full_load_path}.qu"):
    data = qload(full_load_path)
    print("Caricamento completato con successo.")
else:
    raise FileNotFoundError(f"Errore: Il file {full_load_path}.qu non esiste. "
                            "Assicurati di aver completato la simulazione.")

t_list = data['t_list']
expect_n = data['expect_n']
states = data['states']
det_scan = data['det_scan']
spec = data['spec']
params = data['params']

# Estrazione parametri per comodità
Delta_p = params['Delta_p']
nu = params['nu']
N_vib = params['N_vib']

# --- 2. SETUP FIGURA ---
fig = plt.figure(figsize=(18, 5))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

ax3.plot(t_list, expect_n, 'k-', alpha=0.3, label="Decadimento")
line_n, = ax3.plot([], [], 'ro', markersize=8)
ax3.set_title(r"Decadimento $\langle n \rangle$")
ax3.set_xlabel("Tempo")
ax3.set_ylim(0, max(expect_n) + 1)

def update(frame):
    step = frame * 3 
    if step >= len(t_list): step = len(t_list)-1
    
    ax1.clear()
    ax2.clear()
    
    # Recupero stato e calcolo popolazioni
    rho = states[step]
    # Se rho è un ket, lo trasformiamo in operatore densità per ptrace
    if rho.type == 'ket':
        rho = rho * rho.dag()
    
    n_pop = rho.ptrace(1).diag()
    n_avg = expect_n[step]
    
    # 1. Istogramma stati di Fock
    ax1.bar(range(N_vib), n_pop.real, color='cyan', edgecolor='black')
    ax1.set_title(f"Fock States | <n> = {n_avg:.2f}")
    ax1.set_ylim(0, 1)
    ax1.set_xlabel("n")
    
    # 2. Spettro Blue Detuned
    ax2.plot(det_scan, spec, 'k', lw=2)
    ax2.fill_between(det_scan, spec, color='blue', alpha=0.1)
    ax2.axvline(Delta_p, color='green', ls='--', label="Carrier (Dark)")
    ax2.axvline(Delta_p - nu, color='red', lw=2, label="Red Sideband (Cooling)")
    ax2.axvline(Delta_p + nu, color='blue', lw=2, label="Blue Sideband (Heating)")
    ax2.set_title("Meccanismo EIT (Fano su RSB)")
    ax2.legend(loc='upper right')
    ax2.set_xlabel("Detuning")
    
    # 3. Pallino su grafico decadimento
    line_n.set_data([t_list[step]], [n_avg])
    return line_n,

print("Generazione animazione...")
ani = FuncAnimation(fig, update, frames=len(t_list)//3, interval=50, repeat=False)
plt.tight_layout()
plt.show()