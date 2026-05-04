import os
from qutip import qload
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import config as cfg

def animate_time_evolution():
    # --- 1. PATH SETUP ---
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Path where the simulation saves the data
    results_dir = os.path.join(script_dir, "results_time")
    save_file = os.path.join(results_dir, f"time_evol_{cfg.RUN_NAME}")
    
    # Path to save the output images/animations
    images_dir = os.path.join(script_dir, "images", "time_evolution")
    os.makedirs(images_dir, exist_ok=True)

    # --- 2. LOAD DATA ---
    file_to_load = f"{save_file}.qu"
    if not os.path.exists(file_to_load):
        print(f"Error: The file {file_to_load} does not exist.")
        print("Make sure you have run 'simulation_time_rb.py' first.")
        return

    print(f"Loading data from: {file_to_load}...")
    data = qload(save_file)
    
    t_list = data['t_list']
    n_expect = data['n_expect']
    params = data.get('params', {}) # Retrieve parameters if saved
    
    print(f"Data loaded! Found {len(t_list)} time steps (from t=0 to t={t_list[-1]:.1f}).")

    # --- 3. FIGURE & AXES SETUP ---
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Set fixed axes limits based on the loaded data to prevent the plot from jumping
    ax.set_xlim(0, t_list[-1] + (t_list[-1] * 0.05)) # Add 5% padding to the right
    ax.set_ylim(min(n_expect)-0.5, max(n_expect) + 0.5)
    
    # Formatting and labels
    ax.set_title(f"EIT Cooling Dynamics - {cfg.RUN_NAME}", fontsize=14)
    ax.set_xlabel(r"Time ($\gamma^{-1}$)", fontsize=12)
    ax.set_ylabel(r"Average phonon number $\langle n \rangle$", fontsize=12)
    ax.grid(True, linestyle=':', alpha=0.7)
    
    # Display simulation parameters in a text box if available
    if params:
        param_text = f"$\\Delta_p = {params.get('dp', 0):.1f}$\n$\\eta = {params.get('eta', 0):.2f}$"
        ax.text(0.95, 0.5, param_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='center', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Initialize empty line and a leading point for the animation
    line, = ax.plot([], [], 'b-', linewidth=2, label=r'Evolution of $\langle n \rangle$')
    point, = ax.plot([], [], 'ro', markersize=6) # Red dot leading the line
    ax.legend(loc="upper right")
    
    plt.tight_layout()

    # --- 4. ANIMATION LOGIC ---
    # Target ~200 frames total to keep the GIF rendering fast and file size small
    target_frames = 200
    skip_step = max(1, len(t_list) // target_frames)
    total_frames = (len(t_list) // skip_step) + 1
    
    def init():
        """Initial state of the animation."""
        line.set_data([], [])
        point.set_data([], [])
        return line, point

    def update(frame):
        """Update function called for each frame."""
        # Calculate the actual index in the full data arrays
        idx = frame * skip_step
        if idx >= len(t_list):
            idx = len(t_list) - 1
            
        # Update the line up to the current index
        line.set_data(t_list[:idx+1], n_expect[:idx+1])
        # Update the leading point to the current coordinate
        point.set_data([t_list[idx]], [n_expect[idx]])
        
        return line, point

    print("Generating animation... This might take a few moments.")
    # Create the FuncAnimation object
    ani = FuncAnimation(fig, update, frames=total_frames,
                        init_func=init, blit=True, interval=30, repeat=False)

    # --- 5. SAVING & DISPLAY ---
    # 1. Save as GIF
    gif_path = os.path.join(images_dir, f"cooling_curve_{cfg.RUN_NAME}.gif")
    try:
        # Requires the 'pillow' library
        ani.save(gif_path, writer='pillow', fps=30)
        print(f"Animation successfully saved to: {gif_path}")
    except Exception as e:
        print(f"Failed to save GIF. Error: {e}")
        print("Note: Saving GIFs requires the 'pillow' library. Install it with: pip install pillow")
    
    # 2. Save the final static frame as a PNG for quick referencing
    update(total_frames - 1) # Force the plot to the last frame data
    png_path = os.path.join(images_dir, f"cooling_curve_final_{cfg.RUN_NAME}.png")
    plt.savefig(png_path, dpi=300)
    print(f"Final static frame saved to: {png_path}")

    # 3. Show interactive window
    print("Opening interactive viewer...")
    plt.show()

if __name__ == "__main__":
    animate_time_evolution()