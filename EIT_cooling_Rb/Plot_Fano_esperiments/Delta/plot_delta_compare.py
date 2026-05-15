import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import re
import math

# --- PATH SETUP ---
script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.abspath(os.path.join(script_dir, "..", ".."))

sys.path.append(base_dir)
import config as cfg

def extract_delta(folder_name):
    """Helper function to extract the numeric delta value for proper sorting."""
    match = re.search(r'delta([-+]?[0-9.]+)', folder_name)
    return float(match.group(1)) if match else 0.0

def plot_grid_results():
    results_fano_dir = os.path.join(base_dir, "results_fano")

    try:
        run_folders = [f for f in os.listdir(results_fano_dir) 
                       if os.path.isdir(os.path.join(results_fano_dir, f)) and "fano_profile_delta" in f]
    except FileNotFoundError:
        print(f"Error: The target directory {results_fano_dir} does not exist.")
        return

    # Sort the folders NUMERICALLY based on the extracted Delta_c
    run_folders.sort(key=extract_delta)
    num_runs = len(run_folders)
    
    if num_runs == 0:
        print(f"No folders named 'fano_profile_delta...' found in {results_fano_dir}.")
        return

    print(f"Found {num_runs} simulations. Generating 3 separate Grid plots...")

    # --- DYNAMIC GRID SETUP ---
    cols = 3
    rows = math.ceil(num_runs / cols)

    # Generate colors
    colors = plt.cm.viridis(np.linspace(0, 0.9, num_runs))

    # Initialize the three figures with a grid layout
    fig1, axes1 = plt.subplots(rows, cols, figsize=(cols * 4.5, rows * 3.5))
    fig2, axes2 = plt.subplots(rows, cols, figsize=(cols * 4.5, rows * 3.5))
    fig3, axes3 = plt.subplots(rows, cols, figsize=(cols * 4.5, rows * 3.5))

    # Flatten axes arrays for easy iteration
    axes1_flat = np.atleast_1d(axes1).flatten()
    axes2_flat = np.atleast_1d(axes2).flatten()
    axes3_flat = np.atleast_1d(axes3).flatten()

    # --- DATA LOOP ---
    for idx, run_name in enumerate(run_folders):
        results_dir = os.path.join(results_fano_dir, run_name)
        color = colors[idx]
        
        # Extract the specific delta value from the folder name
        local_delta_c = extract_delta(run_name)
        
        # Format the label nicely (e.g., adding explicit + for positive values)
        delta_label = f"$\Delta_c$ = {local_delta_c:+.1f} $\gamma$"

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

        # --- DYNAMIC AUTO-CENTERING ---
        peak_idx = np.argmax(absorption_e2_only)
        peak_x = det_list[peak_idx]
        
        zoom_span = 0.04
        x_min, x_max = peak_x - zoom_span, peak_x + zoom_span

        # --- Plot 1: Fano Profile (Figure 1) ---
        ax1 = axes1_flat[idx]
        ax1_twin = ax1.twinx()
        
        ax1.plot(det_list, absorption_total, color="gray", ls=":", lw=1.5, alpha=0.3)
        ax1_twin.plot(det_list, absorption_e2_only, color=color, ls="-", lw=2)
        
        # RECALCULATE PROBE CENTER based on local Delta_c and config zeeman offset
        # The X-axis is Probe Detuning, so we must shift it exactly as done in config.py
        local_delta_p_center = local_delta_c - cfg.zeeman_offset
        
        # Add Reference Lines using the correct local_delta_p_center and cfg.nu
        ax1.axvline(local_delta_p_center, color="purple", ls="-", alpha=0.2)
        ax1.axvline(local_delta_p_center + cfg.nu, color="red", lw=1.5, ls="--", alpha=0.4)
        ax1.axvline(local_delta_p_center - cfg.nu, color="blue", lw=1.5, ls="--", alpha=0.4)
        
        # Smart Y-Axis limits based on the local peak inside the zoom window
        mask = (det_list >= x_min) & (det_list <= x_max)
        if np.any(mask):
            local_max_e2 = np.max(absorption_e2_only[mask])
            local_max_tot = np.max(absorption_total[mask])
            # Set limits with a 30% margin at the top so peaks don't hit the ceiling
            ax1_twin.set_ylim(-0.02 * local_max_e2, local_max_e2 * 1.3)
            ax1.set_ylim(0, local_max_tot * 1.3)

        ax1.set_xlim(x_min, x_max) 
        ax1.set_title(delta_label, fontsize=11, fontweight='bold', color=color)
        ax1.grid(True, alpha=0.3)
        
        if idx % cols == 0:
            ax1.set_ylabel("Tot. Abs", color="gray", fontsize=9)
        if idx % cols == cols - 1 or idx == num_runs - 1:
            ax1_twin.set_ylabel("F'=2 Pop", color="black", fontsize=9)

        # --- Plot 2: Ground State Populations (Figure 2) ---
        ax2 = axes2_flat[idx]
        ax2.plot(det_list, pop_unwanted_g2, color=color, ls="--", lw=1.5, alpha=0.8, label="F=2 Unw.")
        ax2.plot(det_list, pop_g1_total, color=color, ls="-", lw=2, label="F=1 Res.")
        
        ax2.set_xlim(x_min, x_max)
        ax2.set_ylim(-0.05, 1.05)
        ax2.set_title(delta_label, fontsize=11, fontweight='bold', color=color)
        ax2.grid(True, linestyle='--', alpha=0.4)
        if idx == 0: 
            ax2.legend(loc="upper left", fontsize='x-small', ncol=2)

        # --- Plot 3: Relative Leakage (Figure 3) ---
        ax3 = axes3_flat[idx]
        relative_leakage = np.divide(
            pop_leakage_e3, 
            absorption_total, 
            out=np.zeros_like(pop_leakage_e3), 
            where=np.array(absorption_total) > 1e-8
        )
        ax3.plot(det_list, relative_leakage * 100, color=color, ls="-", lw=2)
        
        ax3.set_xlim(x_min, x_max)
        
        # Smart Y-axis for leakage to prevent clipping
        if np.any(mask):
            local_max_leak = np.max((relative_leakage * 100)[mask])
            ax3.set_ylim(0, max(10, local_max_leak * 1.2)) # Ensure at least 10% on Y-axis
            
        ax3.set_title(delta_label, fontsize=11, fontweight='bold', color=color)
        ax3.grid(True, alpha=0.3)

    # --- HIDE UNUSED SUBPLOTS ---
    for i in range(num_runs, rows * cols):
        axes1_flat[i].set_visible(False)
        axes2_flat[i].set_visible(False)
        axes3_flat[i].set_visible(False)

    # ==========================================
    # FINAL FORMATTING & SAVING
    # ==========================================

    # Add a global legend to Figure 1 for the sidebands
    from matplotlib.lines import Line2D
    custom_lines = [
        Line2D([0], [0], color="purple", lw=1.5, ls="-", alpha=0.5),
        Line2D([0], [0], color="red", lw=2, ls="--", alpha=0.5),
        Line2D([0], [0], color="blue", lw=2, ls="--", alpha=0.5)
    ]
    fig1.legend(custom_lines, ['Default EIT Resonance', 'Red Sideband (Target)', 'Blue Sideband'], 
                loc='upper center', bbox_to_anchor=(0.5, 0.98), ncol=3, fontsize=10)

    fig1.suptitle("Fano Profile Evolution over Coupling Detunings ($\Delta_c$)", fontsize=16, y=1.04)
    fig1.supxlabel(r"Probe Detuning $\Delta_p$ ($\gamma$ units) - Auto-Centered", fontsize=12)
    fig1.tight_layout()

    fig2.suptitle("Ground State Pumping Efficiency over Coupling Detunings", fontsize=16, y=1.02)
    fig2.supxlabel(r"Probe Detuning $\Delta_p$ ($\gamma$ units) - Auto-Centered", fontsize=12)
    fig2.tight_layout()

    fig3.suptitle("Off-Resonant F'=3 Leakage over Coupling Detunings", fontsize=16, y=1.02)
    fig3.supxlabel(r"Probe Detuning $\Delta_p$ ($\gamma$ units) - Auto-Centered", fontsize=12)
    fig3.tight_layout()

    # SAVE ALL THREE FIGURES
    path_fig1 = os.path.join(script_dir, "plot_1_fano_grid.png")
    path_fig2 = os.path.join(script_dir, "plot_2_ground_grid.png")
    path_fig3 = os.path.join(script_dir, "plot_3_leakage_grid.png")

    fig1.savefig(path_fig1, dpi=300, bbox_inches='tight')
    fig2.savefig(path_fig2, dpi=300, bbox_inches='tight')
    fig3.savefig(path_fig3, dpi=300, bbox_inches='tight')

    print(f"Saved: {path_fig1}")
    print(f"Saved: {path_fig2}")
    print(f"Saved: {path_fig3}")
    
    plt.close('all')
    print("All grid plots generated successfully.")

if __name__ == "__main__":
    plot_grid_results()