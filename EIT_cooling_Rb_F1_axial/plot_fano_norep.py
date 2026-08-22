import os
import numpy as np
import matplotlib.pyplot as plt
import config as cfg

def plot_results():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Aggiunto il suffisso per caricare i dati corretti e non sovrascrivere le immagini
    run_name_no_repumper = f"{cfg.RUN_NAME}_no_repumper"
    results_dir = os.path.join(script_dir, "results_fano", run_name_no_repumper)
    images_dir = os.path.join(script_dir, "images", "fano_images")
    os.makedirs(images_dir, exist_ok=True)

    # --- LOAD DATA ---
    try:
        det_list = np.load(os.path.join(results_dir, "det_list.npy"))
        excited_pop_total = np.load(os.path.join(results_dir, "abs_tot.npy"))
    except FileNotFoundError:
        print(f"Error: Could not find results for '{run_name_no_repumper}'. Run fano.py first.")
        return

    # Increase global font sizes for all figures
    base_fs = 18
    title_fs = 20
    plt.rcParams.update({
        'font.size': base_fs,
        'axes.titlesize': title_fs,
        'axes.labelsize': base_fs,
        'xtick.labelsize': base_fs - 2,
        'ytick.labelsize': base_fs - 2,
        'legend.fontsize': base_fs - 2
    })

    # --- FIG A: Zoomed profile around carrier and sidebands ---
    fig_zoom, ax_zoom = plt.subplots(1, 1, figsize=(12, 6))
    ax_zoom.plot(det_list, excited_pop_total, color='black', lw=2.5, alpha=0.95, label='Excited State Population')
    ax_zoom.axvline(cfg.Delta_p_center, color='purple', ls='--', alpha=0.9, label='Carrier')
    ax_zoom.axvline(cfg.Delta_p_center + cfg.nu, color='red', ls=':', alpha=0.9, label='Red Sideband')
    ax_zoom.axvline(cfg.Delta_p_center - cfg.nu, color='blue', ls=':', alpha=0.9, label='Blue Sideband')
    ax_zoom.set_xlabel(r"Probe Detuning $\Delta_p$ [$\gamma$]")
    ax_zoom.set_ylabel('Excited State Population [a.u.]')
    ax_zoom.set_title(f"Excited State Population (Zoom): {run_name_no_repumper}")
    ax_zoom.grid(True, alpha=0.3)
    
    # Zoom limits: shift a bit to the right of center
    pad = 0.03 * cfg.gamma
    x_min = cfg.Delta_p_center - cfg.nu - 0.02 * cfg.gamma
    x_max = cfg.Delta_p_center + cfg.nu + pad
    ax_zoom.set_xlim(x_min, x_max)
    ax_zoom.tick_params(axis='both', which='major', labelsize=base_fs - 2)
    ax_zoom.legend(loc='best')
    out_zoom = os.path.join(images_dir, f"plot_{run_name_no_repumper}.png")
    fig_zoom.tight_layout()
    fig_zoom.savefig(out_zoom, dpi=300, bbox_inches='tight')
    print(f"Zoomed excited state population plot saved to: {out_zoom}")

    # --- FIG B: Full profile with no x-limit ---
    fig_full_unlimited, ax_full_unlimited = plt.subplots(1, 1, figsize=(14, 6))
    ax_full_unlimited.plot(det_list, excited_pop_total, color='black', lw=2.5, alpha=0.95, label='Excited State Population')
    ax_full_unlimited.axvline(cfg.Delta_p_center, color='purple', ls='--', alpha=0.9, label='Carrier')
    ax_full_unlimited.axvline(cfg.Delta_p_center + cfg.nu, color='red', ls=':', alpha=0.9, label='Red Sideband')
    ax_full_unlimited.axvline(cfg.Delta_p_center - cfg.nu, color='blue', ls=':', alpha=0.9, label='Blue Sideband')
    ax_full_unlimited.set_xlabel(r"Probe Detuning $\Delta_p$ [$\gamma$]")
    ax_full_unlimited.set_ylabel('Excited State Population [a.u.]')
    ax_full_unlimited.set_title(f"Excited State Population Profile: {run_name_no_repumper}")
    ax_full_unlimited.grid(True, alpha=0.3)
    ax_full_unlimited.tick_params(axis='both', which='major', labelsize=base_fs - 2)
    ax_full_unlimited.legend(loc='best')
    out_full_unl = os.path.join(images_dir, f"fano_full_unlimited_{run_name_no_repumper}.png")
    fig_full_unlimited.tight_layout()
    fig_full_unlimited.savefig(out_full_unl, dpi=300, bbox_inches='tight')
    print(f"Full (unlimited) profile plot saved to: {out_full_unl}")

    plt.show()

if __name__ == "__main__":
    plot_results()