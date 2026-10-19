import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.lines import Line2D

# Path Setup
script_dir = os.path.dirname(os.path.abspath(__file__))
images_dir = os.path.join(script_dir, "images")
os.makedirs(images_dir, exist_ok=True)

def plot_consistent_eit_diagram():
    # 1. PHYSICAL PARAMETERS
    y_f1 = 0
    y_f2 = 140
    y_excited_base = 550 
    e_offsets = {0: -229.16, 1: -156.94, 2: 0.0, 3: 266.65}
    y_e = {i: y_excited_base + e_offsets[i] for i in range(4)}
    
    # 2. COLORI (Campionati dalla tua immagine)
    c_control = "#1F497D"  # Dark Blue
    c_probe = "#A52A2A"    # Dark Red
    c_repump = "#DAA520"   # Goldenrod
    c_dark = "#2E8B57"     # Sea Green
    c_levels = "#111111"   # Black
    bg_color = "#F5F5F2"   # Off-white background
    
    fig, ax = plt.subplots(figsize=(11, 7), facecolor=bg_color)
    ax.set_facecolor(bg_color)

    # 3. DISEGNO DEI LIVELLI (MANIFOLDS)
    def draw_manifold(y, f_val, label):
        for m in range(-f_val, f_val + 1):
            ax.hlines(y, m-0.35, m+0.35, colors=c_levels, lw=2.5)
        ax.text(3.5, y, label, va='center', fontweight='bold', fontsize=11, color='#333')

    draw_manifold(y_f1, 1, "F=1")
    draw_manifold(y_f2, 2, "F=2") 
    for i in range(4):
        draw_manifold(y_e[i], i, f"F'={i}")

    # 4. MARKERS "DARK PAIR" (Cerchi Verdi)
    circle1 = plt.Circle((1, y_f1), 0.15, color=c_dark, fill=False, lw=3, zorder=3)
    circle2 = plt.Circle((1, y_f2), 0.15, color=c_dark, fill=False, lw=3, zorder=3)
    ax.add_patch(circle1)
    ax.add_patch(circle2)

    # 5. DISEGNO DEI LASER (Frecce)
    # PROBE: F=1, m=1 -> F'=2, m=2
    ax.annotate("", xy=(2, y_e[2]), xytext=(1, y_f1+10),
                arrowprops=dict(arrowstyle="-|>,head_width=0.4,head_length=0.6", color=c_probe, lw=2, zorder=2))
    
    # CONTROL: F=2, m=1 -> F'=2, m=2
    ax.annotate("", xy=(1.95, y_e[2]), xytext=(1.05, y_f2+10),
                arrowprops=dict(arrowstyle="-|>,head_width=0.4,head_length=0.6", color=c_control, lw=3, zorder=2))

    # REPUMPER: F=2, m=2 -> F'=2, m=1
    ax.annotate("", xy=(1, y_e[2]), xytext=(2, y_f2+10),
                arrowprops=dict(arrowstyle="-|>,head_width=0.4,head_length=0.6", color=c_repump, lw=2, ls='--', zorder=2))

    # 6. ANNOTAZIONI LATERALI
    ax.text(-3.4, y_f1 + (y_f2-y_f1)/2, "6.835\nGHz", ha='center', va='center', fontsize=10, color='#444')
    ax.annotate('', xy=(-3.0, y_f2), xytext=(-3.0, y_f1), arrowprops=dict(arrowstyle='<->', color='gray'))
    ax.text(-3.5, y_f1, "$5S_{1/2}$", ha='center', va='center', fontsize=11)
    ax.text(-3.5, y_e[2], "$5P_{3/2}$\n(F'=0..3,\n~0.5 GHz)", ha='center', va='center', fontsize=10)

    # 7. FORMATTAZIONE ASSI
    ax.set_xlabel("magnetic sublevel m", fontsize=12)
    ax.set_xlim(-4.2, 4.2)
    ax.set_ylim(-60, y_e[3] + 80)
    ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1.2)
    
    # 8. LEGENDA
    custom_lines = [
        Line2D([0], [0], color=c_control, lw=3),
        Line2D([0], [0], color=c_repump, lw=2, ls='--'),
        Line2D([0], [0], color=c_probe, lw=2),
        Line2D([0], [0], marker='o', color='w', markeredgecolor=c_dark, markerfacecolor='w', markersize=12, markeredgewidth=2.5)
    ]
    ax.legend(custom_lines, ['control (strong)', 'repump ($\sigma^-$)', 'probe (weak)', 'dark pair'],
              loc='upper left', ncol=2, frameon=True, facecolor='white', edgecolor='#CCC')

    plt.tight_layout()
    output_path = os.path.join(images_dir, "eit_diagram_updated.png")
    plt.savefig(output_path, dpi=300)
    print(f"Diagram saved to: {output_path}")
    plt.show()

if __name__ == "__main__":
    plot_consistent_eit_diagram()