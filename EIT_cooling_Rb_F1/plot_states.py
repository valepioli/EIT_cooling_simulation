import os
import numpy as np
import matplotlib.pyplot as plt
import config as cfg

def plot_fano_results():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.join(script_dir, "results_fano", cfg.RUN_NAME)
    
    # --- CREAZIONE CARTELLA IMMAGINI ---
    # Creiamo 'fano_images' dentro a 'images'
    images_dir = os.path.join(script_dir, "images", "fano_images")
    
    if not os.path.exists(results_dir):
        print(f"Error: No results found for RUN_NAME '{cfg.RUN_NAME}'.")
        return
        
    # Assicurati che l'intero percorso (images/fano_images) esista
    os.makedirs(images_dir, exist_ok=True)

    # Load the common x-axis (detuning)
    det_list = np.load(os.path.join(results_dir, "det_list.npy"))
    
    # --- PLOT 1: Global Absorptions ---
    plt.figure(figsize=(8, 5))
    abs_tot = np.load(os.path.join(results_dir, "abs_tot.npy"))
    abs_e2 = np.load(os.path.join(results_dir, "abs_e2.npy"))
    leak_e3 = np.load(os.path.join(results_dir, "leak_e3.npy"))
    
    plt.plot(det_list, abs_tot, label="Total Absorption (All excited)", color='black', linewidth=2)
    plt.plot(det_list, abs_e2, label="Absorption in e2", linestyle='--')
    plt.plot(det_list, leak_e3, label="Leakage to e3", linestyle=':')
    
    plt.title("Global Steady-State Absorptions")
    plt.xlabel("Probe Detuning")
    plt.ylabel("Population")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # SALVATAGGIO PLOT 1 
    plot1_path = os.path.join(images_dir, f"{cfg.RUN_NAME}_plot_1_global_absorptions.png")
    plt.savefig(plot1_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {plot1_path}")
    plt.show()

    # --- PLOT 2: Ground States (g1 and g2) ---
    plt.figure(figsize=(10, 6))
    for label, f, m in cfg.atom_labels:
        if label in ["g1", "g2"]:
            filename = os.path.join(results_dir, f"pop_{label}_m{m}.npy")
            if os.path.exists(filename):
                data = np.load(filename)
                # Plot only if the state is actually populated
                if np.max(data) > 1e-4: 
                    plt.plot(det_list, data, label=f"{label} (F={f}, m={m})")

    plt.title("Ground State Populations (g1 & g2)")
    plt.xlabel("Probe Detuning")
    plt.ylabel("Population")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # SALVATAGGIO PLOT 2
    plot2_path = os.path.join(images_dir, f"{cfg.RUN_NAME}_plot_2_ground_states.png")
    plt.savefig(plot2_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {plot2_path}")
    plt.show()

    # --- PLOT 3: Excited States (e) ---
    plt.figure(figsize=(10, 6))
    for label, f, m in cfg.atom_labels:
        if label.startswith("e"):
            filename = os.path.join(results_dir, f"pop_{label}_m{m}.npy")
            if os.path.exists(filename):
                data = np.load(filename)
                # Plot only if the state is actually populated
                if np.max(data) > 1e-6: 
                    plt.plot(det_list, data, label=f"{label} (F={f}, m={m})")

    plt.title("Excited State Populations")
    plt.xlabel("Probe Detuning")
    plt.ylabel("Population")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # SALVATAGGIO PLOT 3
    plot3_path = os.path.join(images_dir, f"{cfg.RUN_NAME}_plot_3_excited_states.png")
    plt.savefig(plot3_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {plot3_path}")
    plt.show()

if __name__ == "__main__":
    print("Loading data and generating Fano plots...")
    plot_fano_results()