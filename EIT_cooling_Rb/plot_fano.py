import os
import numpy as np
import matplotlib.pyplot as plt
import config as cfg

def plot_results():
    # --- PATH SETUP ---
    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.join(script_dir, "results_fano", cfg.RUN_NAME)
    images_dir = os.path.join(script_dir, "images", "fano_images")
    os.makedirs(images_dir, exist_ok=True)

    # --- LOAD DATA ---
    try:
        det_list = np.load(os.path.join(results_dir, "det_list.npy"))
        absorption_total = np.load(os.path.join(results_dir, "abs_tot.npy"))
        absorption_e2_only = np.load(os.path.join(results_dir, "abs_e2.npy"))
        pop_unwanted_g2 = np.load(os.path.join(results_dir, "pop_unw_g2.npy"))
        pop_g1_total = np.load(os.path.join(results_dir, "pop_g1.npy"))
        pop_leakage_e3 = np.load(os.path.join(results_dir, "leak_e3.npy"))
    except FileNotFoundError:
        print(f"Error: Could not find results for '{cfg.RUN_NAME}'. Run fano.py first.")
        return

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

    # --- Plot 1: Fano Profile & Sidebands ---
    
    # 1. Plot the total absorption (the background "noise" dominated by F'=3)
    line1 = ax1.plot(
        det_list,
        absorption_total,
        "k-",
        lw=2,
        alpha=0.2, # Very transparent to avoid distraction
        label="Total Absorption (Raw)",
    )
    ax1.set_ylabel("Total Excited Pop", color="gray")
    ax1.tick_params(axis='y', labelcolor="gray")

    # 2. Create a twin Y axis to highlight the tiny F'=2 signal
    ax1_twin = ax1.twinx()
    line2 = ax1_twin.plot(
        det_list,
        absorption_e2_only,
        "k-",
        lw=1.5,
        label="F'=2 Absorption (Filtered Fano)",
    )
    ax1_twin.set_ylabel("F'=2 Population", color="black")
    ax1_twin.tick_params(axis='y', labelcolor="black")

    # Combine the legends from both axes
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc="upper right")

    # Highlight EIT Resonance
    ax1.axvline(
        cfg.Delta_p_center,
        color="purple",
        ls="--",
        alpha=0.7,
        label=r"EIT Resonance",
    )

    # Explicit Sidebands for Blue Detuning
    ax1.axvline(
        cfg.Delta_p_center + cfg.nu,
        color="red",
        lw=2,
        ls="-.",
        label=r"Red Sideband ",
    )
    ax1.axvline(
        cfg.Delta_p_center - cfg.nu,
        color="blue",
        lw=2,
        ls="-.",
        label=r"Blue Sideband",
    )

    ax1.set_title(f"EIT / Fano Profile: {cfg.RUN_NAME}")
    ax1.grid(True, alpha=0.3)


    # --- Plot 2: Ground State Populations & Optical Pumping Efficiency ---
    
    # Plot the population of other mF levels in the F=2 manifold
    # High values here indicate poor optical pumping (atoms stuck in unwanted dark states).
    ax2.plot(
        det_list,
        pop_unwanted_g2,
        "g--",
        alpha=0.8,
        label="Spurious F=2 levels (Optical Pumping Leaks)",
    )

    # Plot the total population remaining in the F=1 manifold
    # This shows if the Repumper laser is working effectively.
    ax2.plot(
        det_list,
        pop_g1_total,
        "r:",
        lw=1.5,
        label="Residual F=1 Pop. (Repump inefficient)",
    )

    ax2.set_ylabel("Population Probability")
    ax2.set_title("Ground State Dynamics: Optical Pumping and State Preparation")
    ax2.set_ylim(-0.05, 1.05) 
    ax2.grid(True, linestyle='--', alpha=0.4)
    ax2.legend(loc="center right", fontsize='small')


    # --- Plot 3: Relative Leakage (Efficiency Loss) ---
    # Calculation of ratio to see how F'=3 influences the total absorption
    relative_leakage = np.divide(
        pop_leakage_e3, 
        absorption_total, 
        out=np.zeros_like(pop_leakage_e3), 
        where=np.array(absorption_total) > 1e-8
    )

    ax3.plot(
        det_list,
        relative_leakage * 100, 
        "m-",
        lw=2,
        label="Fraction of excitation lost to $F'=3$",
    )
    ax3.set_xlabel(r"Probe Detuning $\Delta_p$ ($\gamma$ units)")
    ax3.set_ylabel("Leakage Ratio (%)")
    ax3.set_title("Efficiency Loss: Relative contribution of off-resonant $F'=3$")
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    plt.tight_layout()

    # Save image with specific run name
    output_path = os.path.join(images_dir, f"plot_{cfg.RUN_NAME}.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')

    print(f"Plot saved to: {output_path}")
    plt.show()

if __name__ == "__main__":
    plot_results()