import os
import matplotlib.pyplot as plt
from qutip import qload
import config as cfg

# --- PATH SETUP ---
script_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(script_dir, "results_time")
# Make sure this matches the exact name you used in the simulation script
save_file = os.path.join(results_dir, f"time_evol_500k_MC_{cfg.RUN_NAME}")

def plot_monte_carlo_results():
    print(f"Loading data from {save_file}.qu ...")
    
    # 1. Load the data saved by the simulation
    try:
        data = qload(save_file)
    except FileNotFoundError:
        print("\n[!] Error: Checkpoint file not found.")
        print("Please run the Monte Carlo simulation first to generate the data.")
        return

    # Extract arrays and parameters
    t_list = data['t_list']
    n_expect = data['n_expect']
    params = data.get('params', {})
    
    n_traj = params.get('ntraj', 'Unknown')
    eta = params.get('eta', 'Unknown')
    
    print("\n[+] Data loaded successfully. Generating plot...")

    # 2. Create the plot
    plt.figure(figsize=(10, 6))
    
    # Plotting the expectation value of the phonon number
    plt.plot(t_list, n_expect, color='#1f77b4', linewidth=1.5, label=f'Monte Carlo Average ({n_traj} traj)')
    
    # 3. Formatting the plot
    plt.title(f"Time Evolution of Average Phonon Number $\\langle n \\rangle$\n"
              f"Run: {cfg.RUN_NAME} | $\\eta$ = {eta}", fontsize=14)
    plt.xlabel("Time [Arbitrary Units]", fontsize=12)
    plt.ylabel("Average Phonon Number $\\langle n \\rangle$", fontsize=12)
    
    # Add grid, legend, and limits
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc='upper right', fontsize=11)
    
    # Optional: set y-axis lower limit to 0 since phonons cannot be negative
    plt.ylim(bottom=0)
    plt.tight_layout()

    # 4. Save the plot as a PNG image
    save_dir = os.path.join(script_dir, "images/time_evolution")
    plot_filename = os.path.join(save_dir, f"plot_MC_{cfg.RUN_NAME}.png")
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"[+] Plot successfully saved to: {plot_filename}")

    # 5. Display the plot on screen
    plt.show()

if __name__ == "__main__":
    plot_monte_carlo_results()