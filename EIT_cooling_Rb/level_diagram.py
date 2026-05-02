import matplotlib.pyplot as plt
import numpy as np
import os

# Path Setup
script_dir = os.path.dirname(os.path.abspath(__file__))
images_dir = os.path.join(script_dir, "images")
os.makedirs(images_dir, exist_ok=True)

def plot_consistent_eit_diagram():
    # 1. PHYSICAL PARAMETERS (consistent with your simulation)
    gamma_val = 6.067 # MHz
    detuning_mhz = 12 * gamma_val # Delta_c = 12.0 * gamma in your script
    
    e_offsets = {0: -229.16, 1: -156.94, 2: 0.0, 3: 266.65}
    
    # Coordinates
    y_f1 = 0
    y_f2 = 140
    y_excited_base = 550 # F'=2
    y_e = {i: y_excited_base + e_offsets[i] for i in range(4)}
    y_virtual = y_e[2] + detuning_mhz
    
    # Colors (Modern Palette)
    c_probe = "#007AFF"    # Blue (Pi)
    c_coupler = "#FF2D55"  # Pink/Red (Sigma-)
    c_repump = "#FF9500"   # Orange (Sigma-)
    c_levels = "#2C3E50"   # Dark Blue Gray
    
    fig, ax = plt.subplots(figsize=(12, 9), facecolor='white')

    # 2. DRAW MANIFOLDS
    def draw_manifold(y, f_val, label, active_m=None):
        for m in range(-f_val, f_val + 1):
            is_target = (active_m is not None and m == active_m)
            alpha = 1.0 if is_target else 0.2
            lw = 3 if is_target else 1.5
            ax.hlines(y, m-0.3, m+0.3, colors=c_levels, lw=lw, alpha=alpha)
        ax.text(f_val + 0.5, y, label, va='center', fontweight='bold', alpha=0.7)

    draw_manifold(y_f1, 1, "F=1", active_m=-1)
    draw_manifold(y_f2, 2, "F=2 (Ground)", active_m=None) # We highlight specific m later
    
    # Highlight specific ground states from simulation
    ax.hlines(y_f2, -2.3, -1.7, colors=c_probe, lw=4) # |g1> in simulation logic
    ax.hlines(y_f2, -1.3, -0.7, colors=c_coupler, lw=4) # |g2> in simulation logic

    for i in range(4):
        draw_manifold(y_e[i], i, f"F'={i}", active_m=-2 if i==2 else None)

    # 3. DRAW LASERS (STRAIGHT LINES)
    # Virtual State Line
    ax.hlines(y_virtual, -2.5, -0.5, colors='purple', ls='--', lw=1, alpha=0.5)
    ax.text(-2.6, y_virtual, r"$\Delta_c$", color='purple', ha='right', va='center', fontweight='bold')

    # PROBE: F=2, m=-2 -> F'=2, m=-2 (Pi, pol=0)
    ax.annotate("", xy=(-2, y_virtual), xytext=(-2, y_f2+5),
                arrowprops=dict(arrowstyle="->", color=c_probe, lw=3))
    
    # COUPLER: F=2, m=-1 -> F'=2, m=-2 (Sigma-, pol=-1)
    ax.annotate("", xy=(-1.95, y_virtual), xytext=(-1, y_f2+5),
                arrowprops=dict(arrowstyle="->", color=c_coupler, lw=3))
    
    # REPUMPER: F=1, m=-1 -> F'=2, m=-2 (Sigma-, pol=-1)
    ax.annotate("", xy=(-1.9, y_e[2]-5), xytext=(-1, y_f1+5),
                arrowprops=dict(arrowstyle="->", color=c_repump, lw=2.5, ls=(0,(5,2))))

    # 4. ANNOTATIONS
    # Energy differences
    ax.annotate('', xy=(3.2, y_e[3]), xytext=(3.2, y_e[2]), arrowprops=dict(arrowstyle='<->', color='gray'))
    ax.text(3.3, (y_e[3]+y_e[2])/2, "266.65 MHz", color='gray', fontsize=9, va='center')
    
    ax.annotate('', xy=(-3.5, y_f2), xytext=(-3.5, y_f1), arrowprops=dict(arrowstyle='<->', color='black'))
    ax.text(-3, (y_f1+y_f2)/2, "6.834 GHz", ha='right', va='center', fontweight='bold')

    # State Labels
    ax.text(-2, y_f2-35, r"$|m_F=-2\rangle$", color=c_probe, fontsize=11, ha='center', fontweight='bold')
    ax.text(-1, y_f2-35, r"$|m_F=-1\rangle$", color=c_coupler, fontsize=11, ha='center', fontweight='bold')
    ax.text(-2.1, y_e[2]+20, r"$|e\rangle$", color='black', fontsize=12, ha='center', fontweight='bold')

    # Laser Labels
    ax.text(-2.1, (y_f2+y_virtual)/2, "Probe ($\pi$)", color=c_probe, fontweight='bold', ha='right')
    ax.text(-1.1, (y_f2+y_virtual)/2 + 20, "Coupling ($\sigma^-$)", color=c_coupler, fontweight='bold', ha='left')
    ax.text(-0.9, (y_f1+y_e[2])/2, "Repumper\n($\sigma^-$)", color=c_repump, fontweight='bold', ha='left')

    # Formatting
    ax.set_title("Rb87 EIT Cooling Level Diagram", fontsize=16, pad=20)
    ax.set_xlabel("Magnetic Sublevel ($m_F$)", fontsize=12)
    ax.set_ylabel("Energy (Scaled)", fontsize=12)
    ax.set_xlim(-4, 4)
    ax.set_ylim(-80, 900)
    ax.grid(True, axis='y', alpha=0.2)
    
    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=c_probe, lw=3, label='Probe ($\Delta m = 0$)'),
        Line2D([0], [0], color=c_coupler, lw=3, label='Coupler ($\Delta m = -1$)'),
        Line2D([0], [0], color=c_repump, lw=2.5, ls='--', label='Repumper ($\Delta m = -1$)'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=True)

    plt.tight_layout()
    output_path = os.path.join(images_dir, "eit_diagram.png")
    plt.savefig(output_path, dpi=300)
    print(f"Diagram saved to: {output_path}")
    plt.show()

if __name__ == "__main__":
    plot_consistent_eit_diagram()