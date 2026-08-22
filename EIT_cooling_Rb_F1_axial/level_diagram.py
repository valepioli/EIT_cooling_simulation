import matplotlib.pyplot as plt
import numpy as np
import os

# Path Setup
script_dir = os.path.dirname(os.path.abspath(__file__))
images_dir = os.path.join(script_dir, "images")
os.makedirs(images_dir, exist_ok=True)


def plot_consistent_eit_diagram():
    # 1. PHYSICAL PARAMETERS (consistent with config.py)
    gamma_val = 6.067  # MHz
    detuning_mhz = 13.0 * gamma_val  # matches Delta_c = +13.0 * gamma in config.py

    e_offsets = {0: -229.16, 1: -156.94, 2: 0.0, 3: 266.65}

    # Coordinates
    y_f1 = 0
    y_f2 = 140
    y_excited_base = 550  # F'=2
    y_e = {i: y_excited_base + e_offsets[i] for i in range(4)}
    y_virtual = y_e[2] + detuning_mhz

    # Colors (Modern Palette)
    c_probe = "#FF2D55"    # Pink/Red (sigma-)   -- FIX: probe is now sigma-
    c_coupler = "#FF2D55"  # same polarization family as probe (sigma-)
    c_coupler = "#AF52DE"  # Purple (sigma-, coupling), kept distinct from probe for clarity
    c_repump = "#34C759"   # Green (pi)           -- FIX: repump is now pi
    c_levels = "#2C3E50"   # Dark Blue Gray

    fig, ax = plt.subplots(figsize=(12, 9), facecolor='white')

    # 2. DRAW MANIFOLDS
    def draw_manifold(y, f_val, label, active_m_list=None):
        if active_m_list is None:
            active_m_list = []
        for m in range(-f_val, f_val + 1):
            is_target = (m in active_m_list)
            alpha = 1.0 if is_target else 0.2
            lw = 3 if is_target else 1.5
            ax.hlines(y, m - 0.3, m + 0.3, colors=c_levels, lw=lw, alpha=alpha)
        ax.text(f_val + 0.5, y, label, va='center', fontweight='bold', alpha=0.7)

    draw_manifold(y_f1, 1, "F=1 (Coupling Ground)", active_m_list=[-1])
    draw_manifold(y_f2, 2, "F=2 (Probe/Repump Ground)", active_m_list=[-2, -1])

    # Highlight specific ground states from simulation (FIX: probe is now m=-1, repump leak is m=-2)
    ax.hlines(y_f2, -1.3, -0.7, colors=c_probe, lw=4)   # |F=2, m=-1>  (probe)
    ax.hlines(y_f2, -2.3, -1.7, colors=c_repump, lw=4)  # |F=2, m=-2>  (repump leak state)
    ax.hlines(y_f1, -1.3, -0.7, colors=c_coupler, lw=4)  # |F=1, m=-1> (coupling)

    for i in range(4):
        draw_manifold(y_e[i], i, f"F'={i}", active_m_list=[-2] if i == 2 else None)

    # 3. DRAW LASERS (STRAIGHT LINES)
    ax.hlines(y_virtual, -2.5, -0.5, colors='gray', ls='--', lw=1, alpha=0.5)
    ax.text(-2.6, y_virtual, r"$\Delta_c$", color='gray', ha='right', va='center', fontweight='bold')

    # PROBE: F=2, m=-1 -> F'=2, m=-2 (sigma-, pol=-1). FIX: was pi on m=-2.
    ax.annotate("", xy=(-1.95, y_virtual), xytext=(-1.05, y_f2 + 5),
                arrowprops=dict(arrowstyle="->", color=c_probe, lw=3))

    # COUPLING: F=1, m=-1 -> F'=2, m=-2 (sigma-, pol=-1). Unchanged.
    ax.annotate("", xy=(-2.05, y_virtual), xytext=(-1.05, y_f1 + 5),
                arrowprops=dict(arrowstyle="->", color=c_coupler, lw=3))

    # REPUMPER: F=2, m=-2 -> F'=2, m=-2 (pi, pol=0), resonant. FIX: was sigma- on m=-1.
    ax.annotate("", xy=(-2.1, y_e[2] - 2), xytext=(-2.05, y_f2 + 5),
                arrowprops=dict(arrowstyle="->", color=c_repump, lw=3))

    # 4. ANNOTATIONS
    ax.annotate('', xy=(3.2, y_e[3]), xytext=(3.2, y_e[2]), arrowprops=dict(arrowstyle='<->', color='gray'))
    ax.text(3.3, (y_e[3] + y_e[2]) / 2, "266.65 MHz", color='gray', fontsize=9, va='center')

    ax.annotate('', xy=(-3.5, y_f2), xytext=(-3.5, y_f1), arrowprops=dict(arrowstyle='<->', color='black'))
    ax.text(-3.6, (y_f1 + y_f2) / 2, "6.834 GHz", ha='right', va='center', fontweight='bold')

    # State Labels
    ax.text(-1.05, y_f2 - 35, r"$|m_F=-1\rangle$", color=c_probe, fontsize=11, ha='center', fontweight='bold')
    ax.text(-1.05, y_f1 - 35, r"$|m_F=-1\rangle$", color=c_coupler, fontsize=11, ha='center', fontweight='bold')
    ax.text(-2.05, y_f2 - 35, r"$|m_F=-2\rangle$", color=c_repump, fontsize=11, ha='center', fontweight='bold')
    ax.text(-2.0, y_e[2] + 20, r"$|e\rangle$", color='black', fontsize=12, ha='center', fontweight='bold')

    # Laser Labels
    ax.text(-1.1, (y_f2 + y_virtual) / 2, r"Probe ($\sigma^-$)", color=c_probe, fontweight='bold', ha='left')
    ax.text(-2.15, (y_f1 + y_virtual) / 2 - 50, r"Coupling ($\sigma^-$)", color=c_coupler, fontweight='bold', ha='right')
    ax.text(-2.15, (y_f2 + y_e[2]) / 2, r"Repump ($\pi$)", color=c_repump, fontweight='bold', ha='right')

    # Formatting
    ax.set_title("87Rb EIT Level Diagram (probe/coupling $\\sigma^-$, repump $\\pi$)", fontsize=15, pad=20)
    ax.set_xlabel("Magnetic Sublevel ($m_F$)", fontsize=12)
    ax.set_ylabel("Energy (Scaled)", fontsize=12)
    ax.set_xlim(-4, 4)
    ax.set_ylim(-80, 900)
    ax.grid(True, axis='y', alpha=0.2)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=c_probe, lw=3, label=r'Probe ($\sigma^-$, $\Delta m = -1$, scanned)'),
        Line2D([0], [0], color=c_coupler, lw=3, label=r'Coupling ($\sigma^-$, $\Delta m = -1$, detuned)'),
        Line2D([0], [0], color=c_repump, lw=3, label=r'Repump ($\pi$, $\Delta m = 0$, resonant)'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=True)

    plt.tight_layout()
    output_path = os.path.join(images_dir, "eit_diagram.png")
    plt.savefig(output_path, dpi=300)
    print(f"Diagram saved to: {output_path}")


if __name__ == "__main__":
    plot_consistent_eit_diagram()
