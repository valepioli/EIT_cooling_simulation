import os
import numpy as np
import matplotlib.pyplot as plt
import config as cfg

def plot_results():
    # --- 1. DIRECTORY SETUP ---
    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.join(script_dir, "results_fano", cfg.RUN_NAME)
    images_dir = os.path.join(script_dir, "images", "fano_images")
    os.makedirs(images_dir, exist_ok=True)

    # --- 2. LOAD DATA ---
    try:
        det_list = np.load(os.path.join(results_dir, "det_list.npy"))
        absorption_total = np.load(os.path.join(results_dir, "abs_tot.npy"))
        absorption_e2_only = np.load(os.path.join(results_dir, "abs_e2.npy"))
        pop_leakage_e3 = np.load(os.path.join(results_dir, "leak_e3.npy"))
        
        # Load individual m-states
        pop_g1 = {m: np.load(os.path.join(results_dir, f"pop_g1_m{m}.npy")) for m in [-1, 0, 1]}
        pop_g2 = {m: np.load(os.path.join(results_dir, f"pop_g2_m{m}.npy")) for m in [-2, -1, 0, 1, 2]}
    except FileNotFoundError:
        print(f"Error: Could not find results for '{cfg.RUN_NAME}'. Run fano.py first.")
        return

    # Create figure and subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 14), sharex=True)
    
    # Set a main title that includes the spatial intensity fraction context
    fig.suptitle(f"Steady-State EIT Dynamics: {cfg.RUN_NAME}\n"
                 f"Spatial Intensity Fraction: {cfg.intensity_frac:.2f} "
                 f"(Stark Shift Compensated)", fontsize=16, fontweight='bold', y=0.98)

    # --- 3. Plot 1: Fano Profile (Absorption) ---
    line1 = ax1.plot(det_list, absorption_total, "k-", lw=2, alpha=0.2, label="Total Absorption (Raw)")
    ax1.set_ylabel("Total Excited Pop", color="gray")
    
    ax1_twin = ax1.twinx()
    line2 = ax1_twin.plot(det_list, absorption_e2_only, "k-", lw=1.5, label="F'=2 Absorption")
    ax1_twin.set_ylabel("F'=2 Population", color="black")

    # Center Carrier and Motional Sidebands
    vline_carrier = ax1.axvline(cfg.Delta_p_center, color="purple", ls="--", alpha=0.7, label=r"Carrier")
    vline_rsb = ax1.axvline(cfg.Delta_p_center + cfg.nu, color="red", ls=":", alpha=0.7, label=r"Red Sideband")
    vline_bsb = ax1.axvline(cfg.Delta_p_center - cfg.nu, color="blue", ls=":", alpha=0.7, label=r"Blue Sideband")

    # Merge all labels into the legend of the first plot
    lines = line1 + line2 + [vline_carrier, vline_rsb, vline_bsb]
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc="upper left")
    
    ax1.set_title("Carrier and Sidebands Profile")
    ax1.grid(True, alpha=0.3)

    # --- 4. Plot 2: INDIVIDUAL Ground State Populations ---
    # Colors for F=2 (Blues/Greens) and F=1 (Reds/Oranges)
    colors_g2 = ['#1f77b4', '#00ced1', '#2ca02c', '#9467bd', '#e377c2']
    colors_g1 = ['#d62728', '#ff7f0e', '#8c564b']

    # Plot F=2 manifold (Solid lines)
    for i, m in enumerate([-2, -1, 0, 1, 2]):
        ax2.plot(det_list, pop_g2[m], label=f"F=2, m={m}", color=colors_g2[i], ls='-', lw=2)

    # Plot F=1 manifold (Dashed lines)
    for i, m in enumerate([-1, 0, 1]):
        ax2.plot(det_list, pop_g1[m], label=f"F=1, m={m}", color=colors_g1[i], ls='--', lw=2)

    ax2.set_ylabel("Population Probability")
    ax2.set_title("Ground State Dynamics: Optical Pumping into Dark States (t → ∞)")
    ax2.set_ylim(-0.05, 1.05) 
    ax2.grid(True, linestyle='--', alpha=0.4)
    # Put legend outside the plot so it doesn't cover the lines
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small', title="Magnetic Sublevels")

    # --- 5. Plot 3: Relative Leakage ---
    # Safely divide to avoid zero-division errors in regions with no absorption
    relative_leakage = np.divide(
        pop_leakage_e3, absorption_total, 
        out=np.zeros_like(pop_leakage_e3), where=np.array(absorption_total) > 1e-8
    )

    ax3.plot(det_list, relative_leakage * 100, "m-", lw=2, label="Fraction of excitation lost to $F'=3$")
    ax3.set_xlabel(r"Probe Detuning $\Delta_p$ ($\gamma$ units)")
    ax3.set_ylabel("Leakage Ratio (%)")
    ax3.set_title("Efficiency Loss: Relative contribution of off-resonant $F'=3$")
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    plt.tight_layout()
    # Adjust layout to accommodate the global suptitle
    fig.subplots_adjust(top=0.92) 
    
    output_path = os.path.join(images_dir, f"plot_{cfg.RUN_NAME}.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_path}")
    plt.show()

if __name__ == "__main__":
    plot_results()