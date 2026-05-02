import os
from qutip import qload
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# --- 1. PATH CONFIGURATION ---
script_dir = os.path.dirname(os.path.abspath(__file__))

# Data path (Input)
results_path = os.path.join(script_dir, 'results')
filename = 'eit_simulation_results'
full_load_path = os.path.join(results_path, filename)

# Plots path (Output)
plots_path = os.path.join(script_dir, 'plots')
if not os.path.exists(plots_path):
    os.makedirs(plots_path)
    print(f"Folder created: {plots_path}")

# --- 2. DATA LOADING ---
if os.path.exists(f"{full_load_path}.qu"):
    data = qload(full_load_path)
    print(f"Data loaded from: {full_load_path}.qu")
else:
    raise FileNotFoundError(f"File not found: {full_load_path}.qu")

t_list = data['t_list']
expect_n = data['expect_n']
states = data['states']
det_scan = data['det_scan']
spec = data['spec']
params = data['params']

Delta_p = params['Delta_p']
nu = params['nu']
N_vib = params['N_vib']

# --- 3. FIGURE SETUP ---
fig = plt.figure(figsize=(18, 5))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

def update(frame):
    # Skip frames to speed up animation (consistent with frames=len//3)
    step = frame * 3 
    if step >= len(t_list): 
        step = len(t_list) - 1
    
    ax1.clear()
    ax2.clear()
    
    # 1. Fock state histogram
    rho = states[step]
    if rho.type == 'ket': 
        rho = rho * rho.dag()
    n_pop = rho.ptrace(1).diag()
    n_avg = expect_n[step]
    
    ax1.bar(range(N_vib), n_pop.real, color='cyan', edgecolor='black')
    ax1.set_title(f"Step: {step} | <n> = {n_avg:.2f}")
    ax1.set_ylim(0, 1)
    ax1.set_xlabel("n (vibrational level)")
    ax1.set_ylabel("Population")
    
    # 2. Blue-detuned spectrum
    ax2.plot(det_scan, spec, 'k', lw=2)
    ax2.fill_between(det_scan, spec, color='blue', alpha=0.1)
    ax2.axvline(Delta_p, color='green', ls='--', label="Carrier")
    ax2.axvline(Delta_p + nu, color='red', lw=2, label="Red Sideband (Cooling)")
    ax2.axvline(Delta_p - nu, color='blue', lw=2, label="Blue Sideband (Heating)")
    ax2.set_title("EIT Profile")
    ax2.legend(loc='upper right', fontsize='small')
    #ax2.set_xlim(14.25, 15.80)
    ax2.set_xlabel("Detuning")
    
    # 3. Time evolution of <n>
    ax3.clear()
    ax3.plot(t_list, expect_n, 'k-', alpha=0.3)
    ax3.plot(t_list[step], n_avg, 'ro', markersize=8)
    ax3.set_title(r"Evolution $\langle n \rangle$")
    ax3.set_xlabel("Time")
    ax3.set_ylim(0, max(expect_n) + 1)

# --- 4. GENERATION AND SAVING ---

# Create animation
print("Generating animation...")
ani = FuncAnimation(fig, update, frames=len(t_list)//3, interval=50, repeat=False)

# 1. Save as GIF (requires 'pillow': pip install pillow)
gif_path = os.path.join(plots_path, 'cooling_evolution.gif')
try:
    ani.save(gif_path, writer='pillow', fps=20)
    print(f"Animation saved at: {gif_path}")
except Exception as e:
    print(f"Error saving GIF: {e}")

# 2. Save the last frame as a static PNG image
update(len(t_list)//3 - 1)  # Move figure to final state
img_path = os.path.join(plots_path, 'final_state.png')
plt.savefig(img_path, dpi=300)
print(f"Final image saved at: {img_path}")

# Show on screen
plt.tight_layout()
print("Opening interactive window...")
plt.show()