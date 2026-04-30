import matplotlib.pyplot as plt
import numpy as np
import os

# Get directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Create images folder inside it
images_dir = os.path.join(script_dir, "images")
os.makedirs(images_dir, exist_ok=True)

def plot_eit_cooling_levels():
    # Physical spacings from your simulation
    gamma = 1.0
    MHz = 1 / 6.067
    E_e3, E_e2, E_e1, E_e0 = 266.65 * MHz, 0.0, -156.94 * MHz, -229.16 * MHz
    
    fig = plt.figure(figsize=(16, 9))
    plt.style.use('seaborn-v0_8-whitegrid')
    
    ax = fig.add_subplot(111)

    # 1. DEFINE LEVELS
    y_f1 = 0
    y_f2 = 80
    y_excited_base = 450
    y_e = {
        0: y_excited_base + E_e0,
        1: y_excited_base + E_e1,
        2: y_excited_base + E_e2,
        3: y_excited_base + E_e3
    }

    # Draw F=1 (Leakage)
    for m in range(-1, 2):
        ax.hlines(y_f1, m-0.25, m+0.25, colors='gray', lw=2, alpha=0.5)
    ax.text(2.5, y_f1, "F=1 (Leakage)", color='gray', fontweight='bold')

    # Draw F=2 (Ground)
    for m in range(-2, 3):
        is_eit = (m == -2 or m == -1)
        color = 'darkblue' if is_eit else 'black'
        alpha = 1.0 if is_eit else 0.3
        ax.hlines(y_f2, m-0.25, m+0.25, colors=color, lw=4 if is_eit else 2, alpha=alpha)
    ax.text(2.5, y_f2, "F=2 (Ground)", color='black', fontweight='bold')

    # Draw Excited Manifold
    for f_prime in [0, 1, 2, 3]:
        for m in range(-f_prime, f_prime+1):
            is_eit_e = (f_prime == 2 and m == -2)
            color = 'darkred' if is_eit_e else 'black'
            alpha = 1.0 if is_eit_e else 0.2
            ax.hlines(y_e[f_prime], m-0.25, m+0.25, colors=color, lw=4 if is_eit_e else 2, alpha=alpha)
        ax.text(3.5, y_e[f_prime], f"F'={f_prime}", alpha=0.7)

    # Label key states
    ax.text(-2, y_f2 - 25, r"$|g_1\rangle$", color='blue', fontsize=15, fontweight='bold', ha='center')
    ax.text(-1, y_f2 - 25, r"$|g_2\rangle$", color='green', fontsize=15, fontweight='bold', ha='center')
    ax.text(-2, y_e[2] + 15, r"$|e\rangle$", color='red', fontsize=15, fontweight='bold', ha='center')

    # Frequency arrows
    def draw_freq(y1, y2, text, x=4.5):
        ax.annotate('', xy=(x, y1), xytext=(x, y2), 
                    arrowprops=dict(arrowstyle='<->', color='black', lw=1))
        ax.text(x+0.1, (y1+y2)/2, text, va='center', fontsize=10)

    draw_freq(y_e[3], y_e[2], "266.65 MHz")
    draw_freq(y_e[2], y_e[1], "156.94 MHz")

    # Ground splitting
    ax.annotate('', xy=(-3.5, y_f2), xytext=(-3.5, y_f1), 
                arrowprops=dict(arrowstyle='<->', color='blue', lw=1))
    ax.text(-3.6, (y_f1+y_f2)/2, "6.834 GHz", va='center', ha='right',
            fontsize=10, color='blue', fontweight='bold')

    # Lasers
    ax.annotate("", xy=(-2, y_e[2]-5), xytext=(-2, y_f2+5),
                arrowprops=dict(arrowstyle="->", color="blue", lw=3))
    ax.text(-2.2, (y_f2+y_e[2])/2, "Probe ($\\pi$)", color='blue',
            rotation=90, va='center', fontweight='bold')

    ax.annotate("", xy=(-1.95, y_e[2]-5), xytext=(-1, y_f2+5),
                arrowprops=dict(arrowstyle="->", color="green", lw=3))
    ax.text(-1.3, (y_f2+y_e[2])/2 + 50, "Coupler ($\\sigma^-$)",
            color='green', rotation=-70, fontweight='bold')

    ax.annotate("", xy=(0, y_e[2]-8), xytext=(0, y_f1+5),
                arrowprops=dict(arrowstyle="->", color="orange", lw=2, ls='--'))
    ax.text(0.1, y_f1 + 150, "Repumper\n(clears leakage)",
            color='orange', fontweight='bold', fontsize=10)

    # Detuning
    ax.hlines(y_e[2] + 40, -2.5, -0.5, colors='purple', ls='--', lw=2)
    ax.annotate('', xy=(-2.5, y_e[2] + 40), xytext=(-2.5, y_e[2]), 
                arrowprops=dict(arrowstyle='<->', color='purple', lw=2))
    ax.text(-3.5, y_e[2] + 20,
            r"$\Delta \approx 15\gamma$" + "\n(Blue Detuning)",
            color='purple', fontweight='bold')

    # Formatting
    ax.set_title("EIT Cooling $\Lambda$-System in $^{87}Rb$", fontsize=18)
    ax.set_xlabel("Magnetic Sublevel $m_F$", fontsize=14)
    ax.set_ylabel("Relative Energy (MHz scale)", fontsize=14)
    ax.set_xlim(-4, 5)
    ax.set_ylim(-50, 600)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='blue', alpha=0.7, label='|g1>: Cooling Transition ground'),
        Patch(facecolor='green', alpha=0.7, label='|g2>: Control Transition ground'),
        Patch(facecolor='red', alpha=0.7, label='|e>: Common Excited State')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=12)

    plt.tight_layout()

    # Save image in script-relative images folder
    output_path = os.path.join(images_dir, "eit_cooling_levels.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')

    print(f"Image saved to: {output_path}")

    plt.show()


plot_eit_cooling_levels()