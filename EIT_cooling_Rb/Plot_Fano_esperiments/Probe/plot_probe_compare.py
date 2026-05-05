import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import re

# --- PATH SETUP ---
# Current directory of this script 
script_dir = os.path.dirname(os.path.abspath(__file__))

# Go back two levels to reach the base project directory
base_dir = os.path.abspath(os.path.join(script_dir, "..", ".."))

# Add the base directory to the Python path so we can import config.py
sys.path.append(base_dir)
import config as cfg

def plot_all_results():
    # Define the target directory where the result folders are stored
    results_fano_dir = os.path.join(base_dir, "results_fano")

    # 1. Find all folders containing "fano_profile_probe" inside the results_fano directory
    try:
        run_folders = [f for f in os.listdir(results_fano_dir) 
                       if os.path.isdir(os.path.join(results_fano_dir, f)) and "fano_profile_probe" in f]
    except FileNotFoundError:
        print(f"Error: The target directory {results_fano_dir} does not exist.")
        return

    # Sort the folders so 0.2 comes before 0.5, etc.
    run_folders.sort()

    if not run_folders:
        print(f"No folders named 'fano_profile_probe...' found in {results_fano_dir}.")
        return

    # Initialize the subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 16), sharex=True)
    ax1_twin = ax1.twinx()

    # Generate a color palette based on the number of folders found
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(run_folders)))

    lines_ax1 = []
    labels_ax1 = []

    print(f"Found {len(run_folders)} simulations. Generating plot...")

    # --- DATA LOOP ---
    for idx, run_name in enumerate(run_folders):
        results_dir = os.path.join(results_fano_dir, run_name)
        color = colors[idx]
        
        # Extract the numeric value from the folder name for the legend
        match = re.search(r'probe([0-9.]+)', run_name)
        probe_label = f"Probe = {match.group(1)} $\gamma$" if match else run_name

        # --- LOAD DATA ---
        try:
            det_list = np.load(os.path.join(results_dir, "det_list.npy"))
            absorption_total = np.load(os.path.join(results_dir, "abs_tot.npy"))
            absorption_e2_only = np.load(os.path.join(results_dir, "abs_e2.npy"))
            pop_unwanted_g2 = np.load(os.path.join(results_dir, "pop_unw_g2.npy"))
            pop_g1_total = np.load(os.path.join(results_dir, "pop_g1.npy"))
            pop_leakage_e3 = np.load(os.path.join(results_dir, "leak_e3.npy"))
        except FileNotFoundError:
            print(f"Missing numpy files in {run_name}, skipping.")
            continue

        # --- Plot 1: Fano Profile & Sidebands ---
        ax1.plot(det_list, absorption_total, color=color, ls=":", lw=1.5, alpha=0.15)
        l, = ax1_twin.plot(det_list, absorption_e2_only, color=color, ls="-", lw=2)
        lines_ax1.append(l)
        labels_ax1.append(probe_label)

        # --- Plot 2: Ground State Populations ---
        ax2.plot(det_list, pop_unwanted_g2, color=color, ls="--", lw=1.5, alpha=0.8, label=f"{probe_label} (F=2 Unw.)")
        ax2.plot(det_list, pop_g1_total, color=color, ls="-", lw=2, label=f"{probe_label} (F=1 Res.)")

        # --- Plot 3: Relative Leakage ---
        relative_leakage = np.divide(
            pop_leakage_e3, 
            absorption_total, 
            out=np.zeros_like(pop_leakage_e3), 
            where=np.array(absorption_total) > 1e-8
        )
        ax3.plot(det_list, relative_leakage * 100, color=color, ls="-", lw=2, label=probe_label)

    # ==========================================
    # PLOT FORMATTING
    # ==========================================

    ax1.set_ylabel("Total Excited Pop (Background)", color="gray")
    ax1_twin.set_ylabel("F'=2 Population (Solid Lines)", color="black")
    
    # Reference lines (plotted only once)
    ax1.axvline(cfg.Delta_p_center, color="purple", ls="-", alpha=0.3, label="EIT Resonance")
    ax1.axvline(cfg.Delta_p_center + cfg.nu, color="red", lw=2, ls="--", alpha=0.3, label="Red Sideband")
    ax1.axvline(cfg.Delta_p_center - cfg.nu, color="blue", lw=2, ls="--", alpha=0.3, label="Blue Sideband")

    lines_bg, labels_bg = ax1.get_legend_handles_labels()
    ax1.legend(
       lines_bg + lines_ax1,
       labels_bg + labels_ax1,
       loc="upper left",
       fontsize='small',
       ncol=2
    ).set_zorder(10)
    ax1.patch.set_visible(False) 
    ax1.set_xlim(13.385, 13.415) 
    ax1.set_title("EIT / Fano Profile Comparison across Probe Powers")
    ax1.grid(True, alpha=0.3)

    ax2.set_ylabel("Population Probability")
    ax2.set_title("Ground State Dynamics: Optical Pumping Efficiency")
    ax2.set_ylim(-0.01, 0.5) 
    ax2.set_xlim(13.385, 13.415) 
    ax2.grid(True, linestyle='--', alpha=0.4)
    ax2.legend(loc="upper left", fontsize='x-small', ncol=2) 

    ax3.set_xlabel(r"Probe Detuning $\Delta_p$ ($\gamma$ units)")
    ax3.set_ylabel("Leakage Ratio (%)")
    ax3.set_title("Efficiency Loss: Relative contribution of off-resonant F'=3")
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(13.385, 13.415) 
    ax3.legend(loc="lower left", fontsize='small', ncol=2)

    plt.tight_layout()

    # SAVE: Save directly in the current directory (where the script is located)
    output_path = os.path.join(script_dir, "plot_fano_comparison_all.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')

    print(f"Comparative plot saved to: {output_path}")
    plt.show()

if __name__ == "__main__":
    plot_all_results()