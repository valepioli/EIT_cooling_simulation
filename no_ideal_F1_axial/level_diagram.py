import matplotlib.pyplot as plt
import numpy as np
import os

# Path Setup
script_dir = os.path.dirname(os.path.abspath(__file__))
images_dir = os.path.join(script_dir, "images")
os.makedirs(images_dir, exist_ok=True)

def plot_consistent_eit_diagram():
    # 1. PHYSICAL PARAMETERS (consistent with config.py & Conveyor_belt.pdf)
    gamma_val = 6.067 # MHz
    
    # AC Stark Shifts at the trap center (1.2 mK)
    stark_shift_g = -25.0 # MHz (Ground state attractive potential)
    stark_shift_e = 45.0  # MHz (Excited state repulsive potential)
    
    # Effective Detuning from the SHIFTED state
    effective_detuning_mhz = 13.0 * gamma_val 
    
    e_offsets = {0: -229.16, 1: -156.94, 2: 0.0, 3: 266.65}
    
    # Coordinates (Unperturbed)
    y_f1_unpert = 0
    y_f2_unpert = 140
    y_excited_base_unpert = 550 # F'=2 unperturbed
    
    # Coordinates (Stark Shifted)
    y_f1 = y_f1_unpert + stark_shift_g
    y_f2 = y_f2_unpert + stark_shift_g
    y_excited_base = y_excited_base_unpert + stark_shift_e
    y_e = {i: y_excited_base + e_offsets[i] for i in range(4)}
    
    # The virtual state is defined relative to the SHIFTED excited state
    y_virtual = y_e[2] + effective_detuning_mhz
    
    # Colors (Modern Palette)
    c_probe = "#007AFF"    # Blue (Pi)
    c_coupler = "#FF2D55"  # Pink/Red (Sigma-)
    c_repump = "#34C759"   # Green (Sigma-)
    c_levels = "#2C3E50"   # Dark Blue Gray
    c_unpert = "#A0A0A0"   # Light Gray for unperturbed levels
    
    fig, ax = plt.subplots(figsize=(13, 10), facecolor='white')

    # 2. DRAW MANIFOLDS
    def draw_manifold(y, f_val, label, active_m_list=None, is_unperturbed=False):
        if active_m_list is None:
            active_m_list = []
            
        for m in range(-f_val, f_val + 1):
            if is_unperturbed:
                ax.hlines(y, m-0.3, m+0.3, colors=c_unpert, lw=1.5, linestyles='dashed', alpha=0.6)
            else:
                is_target = (m in active_m_list)
                alpha = 1.0 if is_target else 0.2
                lw = 3 if is_target else 1.5
                ax.hlines(y, m-0.3, m+0.3, colors=c_levels, lw=lw, alpha=alpha)
                
        # Only label the solid lines
        if not is_unperturbed:
            ax.text(f_val + 0.5, y, label, va='center', fontweight='bold', alpha=0.7)

    # Draw Unperturbed Levels (Visual Reference for Stark Shift)
    draw_manifold(y_f1_unpert, 1, "", is_unperturbed=True)
    draw_manifold(y_f2_unpert, 2, "", is_unperturbed=True)
    draw_manifold(y_excited_base_unpert, 2, "", is_unperturbed=True)

    # Draw Shifted Levels
    draw_manifold(y_f1, 1, "F=1 (Coupler Ground)", active_m_list=[-1])
    draw_manifold(y_f2, 2, "F=2 (Probe/Repump Ground)", active_m_list=[-2, -1]) 
    
    # Highlight specific ground states from simulation
    ax.hlines(y_f2, -2.3, -1.7, colors=c_probe, lw=4)  # |F=2, m=-2> 
    ax.hlines(y_f2, -1.3, -0.7, colors=c_repump, lw=4) # |F=2, m=-1>
    ax.hlines(y_f1, -1.3, -0.7, colors=c_coupler, lw=4) # |F=1, m=-1>

    for i in range(4):
        draw_manifold(y_e[i], i, f"F'={i}", active_m_list=[-2] if i==2 else None)

    # 3. DRAW LASERS (STRAIGHT LINES)
    # Virtual State Line
    ax.hlines(y_virtual, -2.5, -0.5, colors='purple', ls='--', lw=1.5, alpha=0.6)
    ax.text(-2.6, y_virtual, r"$\Delta_c$", color='purple', ha='right', va='center', fontweight='bold')

    # PROBE: F=2, m=-2 -> F'=2, m=-2 (Pi, pol=0)
    ax.annotate("", xy=(-2.05, y_virtual), xytext=(-2.05, y_f2+5),
                arrowprops=dict(arrowstyle="->", color=c_probe, lw=3))
    
    # COUPLER: F=1, m=-1 -> F'=2, m=-2 (Sigma-, pol=-1)
    ax.annotate("", xy=(-1.95, y_virtual), xytext=(-1.05, y_f1+5),
                arrowprops=dict(arrowstyle="->", color=c_coupler, lw=3))

    # REPUMPER: F=2, m=-1 -> F'=2, m=-2 (Sigma-, pol=-1) -> RESONANT (Points to y_e[2])
    ax.annotate("", xy=(-1.85, y_e[2] - 2), xytext=(-0.95, y_f2+5),
                arrowprops=dict(arrowstyle="->", color=c_repump, lw=3))

    # 4. ANNOTATIONS
    # Energy differences & Shifts
    ax.annotate('', xy=(3.2, y_e[3]), xytext=(3.2, y_e[2]), arrowprops=dict(arrowstyle='<->', color='gray'))
    ax.text(3.3, (y_e[3]+y_e[2])/2, "266.65 MHz", color='gray', fontsize=9, va='center')
    
    # Ground State Hyperfine Splitting
    ax.annotate('', xy=(-3.5, y_f2), xytext=(-3.5, y_f1), arrowprops=dict(arrowstyle='<->', color='black'))
    ax.text(-3.6, (y_f1+y_f2)/2, "6.834 GHz", ha='right', va='center', fontweight='bold')

    # Stark Shift Indicators
    ax.annotate('', xy=(2.7, y_excited_base), xytext=(2.7, y_excited_base_unpert), 
                arrowprops=dict(arrowstyle='->', color='red', lw=1.5))
    ax.text(2.8, (y_excited_base + y_excited_base_unpert)/2, "+45 MHz\n(Stark)", color='red', fontsize=9, va='center')

    ax.annotate('', xy=(2.7, y_f2), xytext=(2.7, y_f2_unpert), 
                arrowprops=dict(arrowstyle='->', color='blue', lw=1.5))
    ax.text(2.8, (y_f2 + y_f2_unpert)/2, "-25 MHz\n(Stark)", color='blue', fontsize=9, va='center')

    # State Labels
    ax.text(-2.05, y_f2-35, r"$|m_F=-2\rangle$", color=c_probe, fontsize=11, ha='center', fontweight='bold')
    ax.text(-1.05, y_f1-35, r"$|m_F=-1\rangle$", color=c_coupler, fontsize=11, ha='center', fontweight='bold')
    ax.text(-0.95, y_f2-35, r"$|m_F=-1\rangle$", color=c_repump, fontsize=11, ha='center', fontweight='bold')
    ax.text(-2.1, y_e[2]+25, r"$|e\rangle$", color='black', fontsize=12, ha='center', fontweight='bold')

    # Laser Labels
    ax.text(-2.1, (y_f2+y_virtual)/2, "Probe ($\pi$)", color=c_probe, fontweight='bold', ha='right')
    ax.text(-1.15, (y_f1+y_virtual)/2 - 50, "Coupling ($\sigma^-$)", color=c_coupler, fontweight='bold', ha='right')
    ax.text(-0.85, (y_f2+y_e[2])/2, "Repump ($\sigma^-$)", color=c_repump, fontweight='bold', ha='left')

    # Formatting
    ax.set_title("87Rb EIT Level Diagram (with 1.2 mK AC Stark Shifts)", fontsize=16, pad=20)
    ax.set_xlabel("Magnetic Sublevel ($m_F$)", fontsize=12)
    ax.set_ylabel("Energy (Scaled)", fontsize=12)
    ax.set_xlim(-4, 4)
    ax.set_ylim(-120, 950)
    ax.grid(True, axis='y', alpha=0.2)
    
    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=c_probe, lw=3, label='Probe ($\Delta m = 0$, Detuned)'),
        Line2D([0], [0], color=c_coupler, lw=3, label='Coupler ($\Delta m = -1$, Detuned)'),
        Line2D([0], [0], color=c_repump, lw=3, label='Repumper ($\Delta m = -1$, Resonant)'),
        Line2D([0], [0], color=c_unpert, lw=1.5, ls='--', label='Unperturbed States (0 mW)'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=True)

    plt.tight_layout()
    output_path = os.path.join(images_dir, "eit_diagram_stark_shifted.png")
    plt.savefig(output_path, dpi=300)
    print(f"Diagram saved to: {output_path}")
    plt.show()

if __name__ == "__main__":
    plot_consistent_eit_diagram()