import os
import matplotlib.pyplot as plt
from qutip import qload
import config as cfg

# --- PATH SETUP ---
script_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(script_dir, "results_time")
save_file = os.path.join(results_dir, f"time_evol_MC_{cfg.RUN_NAME}")

# --- PLOT LIMITS ---
# Change these values to zoom on the two plots independently
X_LIMIT_PHONONS = (0, 500000) 
X_LIMIT_POPULATIONS = (0, 500000)

def plot_monte_carlo_results():
    print(f"Loading data from {save_file}.qu ...")
    
    # 1. Load the data saved by the simulation
    try:
        data = qload(save_file)
    except FileNotFoundError:
        print("\n[!] Error: Checkpoint file not found.")
        print(f"Make sure you have run the simulation and the file exists at: {save_file}.qu")
        return

    # Extract arrays and parameters
    t_list = data['t_list']
    n_expect = data['n_expect']
    populations = data.get('populations', {}) 
    params = data.get('params', {})
    
    n_traj = params.get('ntraj', 'Unknown')
    eta = params.get('eta', 'Unknown')
    
    print("\n[+] Data loaded successfully. Generating plots...")

    # 2. Create a figure with 2 vertically stacked subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    
    # --- NEW: Global Title with Spatial/Stark Context ---
    fig.suptitle(f"EIT Cooling Time Evolution: {cfg.RUN_NAME}\n"
                 f"Spatial Intensity Fraction: {cfg.intensity_frac:.2f} "
                 f"(Stark Shift Compensated)", fontsize=16, fontweight='bold', y=0.98)

    # --- SUBPLOT 1: Phonon Number ---
    ax1.plot(t_list, n_expect, color='#1f77b4', linewidth=1.5, label=f'Average $\\langle n \\rangle$ ({n_traj} traj)')
    
    ax1.set_title(f"Phonon Cooling Dynamics | $\\eta$ = {eta}", fontsize=14)
    ax1.set_xlabel("Time [Arbitrary Units]", fontsize=12) 
    ax1.set_ylabel("Average Phonon Number $\\langle n \\rangle$", fontsize=12)
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend(loc='upper right', fontsize=11)
    
    # Subplot 1 Limits
    ax1.set_ylim(bottom=0)
    ax1.set_xlim(X_LIMIT_PHONONS) 

    # --- SUBPLOT 2: Atomic Populations ---
    for state_label, pop_array in populations.items():
        # Only plot states that actually get populated to avoid visual clutter
        if max(pop_array) > 0.01:
            ax2.plot(t_list, pop_array, linewidth=1.5, label=state_label)

    ax2.set_title("Atomic Level Occupancy Dynamics", fontsize=14)
    ax2.set_xlabel("Time [Arbitrary Units]", fontsize=12)
    ax2.set_ylabel("State Population", fontsize=12)
    ax2.grid(True, linestyle='--', alpha=0.7)
    # Move legend slightly outside to avoid covering data lines
    ax2.legend(loc='center right', fontsize=10, bbox_to_anchor=(1.20, 0.5)) 
    
    # Subplot 2 Limits
    ax2.set_ylim(-0.05, 1.05) 
    ax2.set_xlim(X_LIMIT_POPULATIONS) 

    plt.tight_layout()
    # Adjust top margin so the tight_layout doesn't overwrite the suptitle
    fig.subplots_adjust(top=0.90)

    # 3. Save the plot as a PNG image
    save_dir = os.path.join(script_dir, "images", "time_evolution")
    os.makedirs(save_dir, exist_ok=True)  
    plot_filename = os.path.join(save_dir, f"plot_MC_{cfg.RUN_NAME}.png")
    
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"[+] Plot successfully saved to: {plot_filename}")

    # 4. Display the plot on screen
    plt.show()

if __name__ == "__main__":
    plot_monte_carlo_results()