import os
import numpy as np
import matplotlib.pyplot as plt
from qutip import qload
import config as cfg

# --- IMPOSTAZIONI GRAFICHE (Font e dimensioni più grandi) ---
plt.rcParams.update({
    'font.size': 16,          # Dimensione generale del testo
    'axes.titlesize': 20,     # Dimensione dei titoli dei grafici
    'axes.labelsize': 18,     # Dimensione delle etichette degli assi (X, Y)
    'xtick.labelsize': 16,    # Dimensione dei numeri sull'asse X
    'ytick.labelsize': 16,    # Dimensione dei numeri sull'asse Y
    'legend.fontsize': 14,    # Dimensione del testo nella legenda
    'lines.linewidth': 2.5    # Spessore delle linee
})

# --- PATH SETUP ---
script_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(script_dir, "results_time")

# Aggiunto il suffisso per caricare i dati corretti e non sovrascrivere le immagini
run_name_no_repumper = f"{cfg.RUN_NAME}_no_repumper"
save_file = os.path.join(results_dir, f"time_evol_MC_{run_name_no_repumper}")

# --- PLOT LIMITS (in unità di simulazione) ---
# Cambia questi valori per lo zoom. Verranno convertiti automaticamente in ms nel grafico.
X_LIMIT_PHONONS = (0, 700000) 
X_LIMIT_POPULATIONS = (0, 700000)

# Fattore di conversione: 1 unità = 26.2 ns = 26.2 * 10^-6 ms
TIME_CONV_FACTOR = 26.2 * 1e-6 

def plot_monte_carlo_results():
    print(f"Loading data from {save_file}.qu ...")
    
    # 1. Carica i dati della simulazione
    try:
        data = qload(save_file)
    except FileNotFoundError:
        print("\n[!] Error: Checkpoint file not found.")
        print(f"Make sure you have run the simulation and the file exists at: {save_file}.qu")
        return

    # Estrai array e parametri
    t_list_sim = np.array(data['t_list'])
    
    # Converti il tempo da unità di simulazione a millisecondi
    t_list_ms = t_list_sim * TIME_CONV_FACTOR
    
    n_expect = data['n_expect']
    populations = data.get('populations', {}) 
    params = data.get('params', {})
    
    n_traj = params.get('ntraj', 'Unknown')
    eta = params.get('eta', 'Unknown')
    
    print("\n[+] Data loaded successfully. Generating plots...")

    save_dir = os.path.join(script_dir, "images", "time_evolution")
    os.makedirs(save_dir, exist_ok=True)  

    # ==========================================
    # GRAFICO 1: EVOLUZIONE DEI FONONI
    # ==========================================
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    
    ax1.plot(t_list_ms, n_expect, color='#1f77b4', label=f'Average $\\langle n \\rangle$ ({n_traj} traj)')
    
    ax1.set_title(f"Time Evolution of Phonon Number\nRun: {run_name_no_repumper} | $\\eta$ = {eta}")
    ax1.set_xlabel("Time [ms]")
    ax1.set_ylabel("Average Phonon Number $\\langle n \\rangle$")
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend(loc='upper right')
    
    # Applica i limiti convertiti in ms
    ax1.set_ylim(bottom=0)
    ax1.set_xlim([x * TIME_CONV_FACTOR for x in X_LIMIT_PHONONS]) 

    plt.tight_layout()
    plot1_filename = os.path.join(save_dir, f"plot_MC_phonons_{run_name_no_repumper}.png")
    fig1.savefig(plot1_filename, dpi=300, bbox_inches='tight')
    print(f"[+] Phonon plot saved to: {plot1_filename}")

    # ==========================================
    # GRAFICO 2: POPOLAZIONI ATOMICHE
    # ==========================================
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    
    # Plotta tutte le popolazioni che superano la soglia dell'1%
    for state_label, pop_array in populations.items():
        if max(pop_array) > 0.01:
            ax2.plot(t_list_ms, pop_array, label=state_label)

    # Nota: se nella simulazione senza repumper tracci anche gli stati fondamentali (ground states), 
    # potresti voler cambiare "Excited State Population" in un più generico "State Population" 
    # dato che la popolazione eccitata andrà a 0 e vedrai solo gli stati fondamentali andare a 1.
    ax2.set_title(f"Time Evolution of Atomic Populations\nRun: {run_name_no_repumper} | $\\eta$ = {eta}")
    ax2.set_xlabel("Time [ms]")
    ax2.set_ylabel("State Population")
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    # Sposta la legenda fuori dal grafico
    ax2.legend(loc='center left', bbox_to_anchor=(1.05, 0.5)) 
    
    # Applica i limiti convertiti in ms (Qui lasciamo da -0.05 a 1.05 per vedere gli stati oscuri che si riempiono)
    ax2.set_ylim(-0.05, 1.05) 
    ax2.set_xlim([x * TIME_CONV_FACTOR for x in X_LIMIT_POPULATIONS]) 

    plt.tight_layout()
    plot2_filename = os.path.join(save_dir, f"plot_MC_populations_{run_name_no_repumper}.png")
    fig2.savefig(plot2_filename, dpi=300, bbox_inches='tight')
    print(f"[+] Populations plot saved to: {plot2_filename}")

    # 4. Mostra i grafici a schermo
    plt.show()

if __name__ == "__main__":
    plot_monte_carlo_results()