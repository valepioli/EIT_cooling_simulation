import os
import numpy as np
from qutip import qload
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# --- IMPOSTAZIONI GRAFICHE ---
plt.rcParams.update({
    'font.size': 16,
    'axes.titlesize': 20,
    'axes.labelsize': 18,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 14,
    'lines.linewidth': 2.5
})

# --- 1. PATH CONFIGURATION ---
script_dir = os.path.dirname(os.path.abspath(__file__))

results_path = os.path.join(script_dir, 'results')
filename = 'eit_simulation_results'
full_load_path = os.path.join(results_path, filename)

plots_path = os.path.join(script_dir, 'plots')
if not os.path.exists(plots_path):
    os.makedirs(plots_path)
    print(f"Folder created: {plots_path}")

# --- 2. DATA LOADING ---
if os.path.exists(f"{full_load_path}.qu"):
    data = qload(full_load_path)
    print(f"Data loaded from: {full_load_path}.qu")
else:
    raise FileNotFoundError(f"File not found: {full_load_path}.qu")

# Estrazione dati e conversione esplicita ad array NumPy
t_list_sim = np.array(data['t_list'])
expect_n = data['expect_n']
states = data['states']
det_scan_sim = np.array(data['det_scan'])
spec = data['spec']
params = data['params']

N_vib = params['N_vib']

# --- CONVERSIONI UNITÀ ---
# Tempo: 1 unità = 26.2 ns = 26.2 * 10^-6 ms
TIME_CONV_FACTOR = 26.2 * 1e-6
# For the Fano/EIT spectrum we use a dimensionless detuning axis in units of gamma.
gamma = params.get('gamma', 1.0)

t_list_ms = t_list_sim * TIME_CONV_FACTOR
det_scan_gamma = det_scan_sim / gamma
Delta_p_gamma = params['Delta_p'] / gamma
nu_gamma = params['nu'] / gamma

# --- 3. FIGURE SETUP (Separate) ---
# Fig 1: Istogramma Stati di Fock (Animato)
fig1, ax1 = plt.subplots(figsize=(8, 6))

# Fig 2: Spettro EIT (Statico, non varia nel tempo)
fig2, ax2 = plt.subplots(figsize=(10, 6))

# Fig 3: Evoluzione Temporale (Animata)
fig3, ax3 = plt.subplots(figsize=(10, 6))

# Fig 4: Profilo completo di assorbimento
fig4, ax4 = plt.subplots(figsize=(10, 6))

# --- PLOT GRAFICO STATICO (EIT Profile) ---
ax2.plot(det_scan_gamma, spec, 'k', lw=2)
ax2.fill_between(det_scan_gamma, spec, color='blue', alpha=0.1)
ax2.axvline(Delta_p_gamma, color='green', ls='--', label=r"Carrier")
ax2.axvline(Delta_p_gamma + nu_gamma, color='red', lw=2, label=r"Red Sideband")
ax2.axvline(Delta_p_gamma - nu_gamma, color='blue', lw=2, label=r"Blue Sideband")
ax2.set_title("Excited State Population Profile(3-Level System)")
ax2.legend(loc='upper left')
ax2.set_xlabel(r"Probe detuning $\Delta_p / \gamma$")
ax2.set_ylabel("Excited State Population [a.u.]")
ax2.grid(True, linestyle='--', alpha=0.6)
ax2.set_xlim(Delta_p_gamma - 2, Delta_p_gamma + 2)

fig2.tight_layout()
eit_plot_path = os.path.join(plots_path, 'eit_spectrum_profile.png')
fig2.savefig(eit_plot_path, dpi=300)
print(f"[+] Static EIT profile saved at: {eit_plot_path}")

# --- PLOT PROFILO COMPLETO DI ASSORBIMENTO ---
ax4.plot(det_scan_gamma, spec, 'k', lw=2.5, label='EIT profile')
ax4.fill_between(det_scan_gamma, spec, color='darkorange', alpha=0.25)
ax4.axvline(0, color='black', ls=':', lw=1.5, label=r"$\Delta_p = 0$")
ax4.axvline(Delta_p_gamma, color='green', ls='--', lw=2.0, label=r"Carrier ($\Delta_p / \gamma$)")
ax4.set_title("Excited State Population Profile (3-Level System)")
ax4.set_xlabel(r"Probe detuning $\Delta_p [\gamma$] ")
ax4.set_ylabel("Excited State Population [a.u.]")
ax4.grid(True, linestyle='--', alpha=0.6)
ax4.legend(loc='best')
ax4.set_xlim(-5, 20)

fig4.tight_layout()
full_absorption_path = os.path.join(plots_path, 'full_absorption_profile.png')
fig4.savefig(full_absorption_path, dpi=300)
print(f"[+] Full absorption profile saved at: {full_absorption_path}")

# --- FUNZIONI DI UPDATE PER LE ANIMAZIONI ---
def update_fig1(frame):
    step = frame * 3 
    if step >= len(t_list_sim): 
        step = len(t_list_sim) - 1
    
    ax1.clear()
    rho = states[step]
    if rho.type == 'ket': 
        rho = rho * rho.dag()
    n_pop = rho.ptrace(1).diag()
    n_avg = expect_n[step]
    current_time_ms = t_list_ms[step]
    
    ax1.bar(range(N_vib), n_pop.real, color='cyan', edgecolor='black')
    ax1.set_title(f"Time: {current_time_ms:.4f} ms | $\\langle n \\rangle$ = {n_avg:.2f}")
    ax1.set_ylim(0, 1)
    ax1.set_xlabel("n (vibrational level)")
    ax1.set_ylabel("Population")

def update_fig3(frame):
    step = frame * 3 
    if step >= len(t_list_sim): 
        step = len(t_list_sim) - 1
    
    ax1.clear() # Fix interno plt per la gestione della memoria in animazioni multiple
    ax3.clear()
    n_avg = expect_n[step]
    
    ax3.plot(t_list_ms, expect_n, 'k-', alpha=0.3, label="Evolution")
    ax3.set_title(r"Evolution of $\langle n \rangle$")
    ax3.set_xlabel("Time [ms]")
    ax3.set_ylabel(r"Average Phonon Number $\langle n \rangle$")
    ax3.set_ylim(0, max(expect_n) + 1)
    ax3.grid(True, linestyle='--', alpha=0.6)
    ax3.legend(loc='upper right')

# --- 4. GENERATION AND SAVING ---
print("Generating animations...")
tot_frames = len(t_list_sim) // 3

ani1 = FuncAnimation(fig1, update_fig1, frames=tot_frames, interval=50, repeat=False)
ani3 = FuncAnimation(fig3, update_fig3, frames=tot_frames, interval=50, repeat=False)

gif1_path = os.path.join(plots_path, 'fock_state_evolution.gif')
gif3_path = os.path.join(plots_path, 'phonon_time_evolution.gif')

try:
    ani1.save(gif1_path, writer='pillow', fps=20)
    ani3.save(gif3_path, writer='pillow', fps=20)
    print(f"[+] Animations saved:\n    -> {gif1_path}\n    -> {gif3_path}")
except Exception as e:
    print(f"Error saving GIFs: {e}")

# Salva l'ultimo stato (frame finale) come immagini statiche PNG
update_fig1(tot_frames - 1)
update_fig3(tot_frames - 1)

final_fock_path = os.path.join(plots_path, 'final_fock_state.png')
final_time_path = os.path.join(plots_path, 'final_phonon_evolution.png')

fig1.tight_layout()
fig3.tight_layout()
fig1.savefig(final_fock_path, dpi=300)
fig3.savefig(final_time_path, dpi=300)
print(f"[+] Final static images saved:\n    -> {final_fock_path}\n    -> {final_time_path}")

# Mostra i grafici interattivi
print("Opening interactive windows...")
plt.show()