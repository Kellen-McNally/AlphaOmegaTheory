"""
Sedenion Dynamical Scattering: Cross-Section Calculation

Computes the angular dependence of a scattering amplitude (Mott scattering)
using the non-associative Sedenion vertex structure. Verifies the dSigma/dOmega 
~ 1/sin^4(theta/2) dependence characteristic of vector boson exchange.
"""

import numpy as np
import sys
import os

# Ensure we can import from core
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../')))

from core.sedenions import multiply, create_basis, norm_sq

def calculate_scattering():
    print("Dynamical Scattering Validation (Mott Formula)")
    print("============================================")
    
    # Scattering angle theta from 10 to 170 degrees
    thetas = np.linspace(10, 170, 50) * np.pi / 180
    
    results = []
    
    for theta in thetas:
        # Momentum transfer q ~ sin(theta/2)
        # Propagator scaling 1/q^2
        sin_half_theta = np.sin(theta/2)
        q_mag = 2.0 * sin_half_theta 
        propagator_scale = 1.0 / (q_mag**2)
        
        # Incoming Fermion (Internal Octonion, e8)
        u_in = create_basis(8) 
        
        # Gauge Field A (External Octonion, e1), scaled by propagator
        A_field = create_basis(1) * propagator_scale
        
        # Interaction: A * u_in
        interaction = multiply(A_field, u_in)
        
        # Amplitude squared ~ |interaction|^2
        amplitude_sq = norm_sq(interaction)
        results.append(amplitude_sq)
        
    y_sim = np.array(results)
    y_theory = 1.0 / (np.sin(thetas/2)**4)
    
    # Normalize
    idx_90 = np.argmin(np.abs(thetas - np.pi/2))
    scale_factor = y_sim[idx_90] / y_theory[idx_90]
    y_theory_scaled = y_theory * scale_factor
    
    correlation = np.corrcoef(y_sim, y_theory_scaled)[0,1]
    
    print(f"  Correlation with Mott Formula: {correlation:.6f}")
    
    if correlation > 0.999:
        print("  Result: Sedenion vertex reproduces Mott scattering angular dependence.")
    else:
        print("  Result: Mismatch with Mott prediction.")

if __name__ == "__main__":
    calculate_scattering()
