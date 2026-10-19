import numpy as np
import config as cfg
from simulation import run_monte_carlo_evolution

# Define the range of intensity fractions based on the thermal spread
# Tr ~ 350 uK, U0 = 1200 uK -> atoms explore down to ~ 1 - (350/1200) = 0.71 of the max depth
min_intensity = 0.71
max_intensity = 1.0
num_spatial_classes = 5

# Create an array of intensity fractions to simulate
intensity_fractions = np.linspace(min_intensity, max_intensity, num_spatial_classes)

# Optional: Define statistical weights (e.g., a simple Boltzmann weight or Gaussian profile)
# For a uniform approximation:
weights = np.ones(num_spatial_classes) / num_spatial_classes 

print(f"--- Starting Inhomogeneous Broadening Simulation ---")
print(f"Running {num_spatial_classes} spatial classes from {min_intensity} to {max_intensity}\n")

for i, frac in enumerate(intensity_fractions):
    print(f"==================================================")
    print(f"Running Spatial Class {i+1}/{num_spatial_classes} | Intensity Fraction: {frac:.3f}")
    print(f"==================================================")
    
    # 1. Update the configuration parameter dynamically
    cfg.intensity_frac = frac
    
    # 2. Update the run name so results don't overwrite each other
    cfg.RUN_NAME = f"eit_axial_frac_{frac:.3f}"
    
    # 3. Run the Monte Carlo evolution for this specific spatial class
    run_monte_carlo_evolution()

print("\n--- All spatial classes simulated successfully! ---")
# After this, you would load the saved .qu files and compute the weighted average 
# of the populations to see the full broadened spectrum.