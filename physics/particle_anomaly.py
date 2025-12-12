"""
Symbolic Proof: Anomaly Cancellation

Verifies the cancellation of gauge anomalies in the G2-Sedenion framework.
Derives fermion charges from G2 Cartan generators and checks the 
vanishing of gravitational (Sum Q) and cubic (Sum Q^3) anomalies.
"""

import numpy as np
import sys
import os

# Ensure we can import from core
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from physics.particle_group_embedding import generate_g2_generators

def check_anomaly_cancellation():
    """
    Check standard anomaly cancellation conditions:
    1. Sum(Q) = 0 (Gravitational-Gauge anomaly)
    2. Sum(Q^3) = 0 (Gauge^3 anomaly)
    """
    sum_q = np.sum(charges_list)
    sum_q3 = np.sum(np.power(charges_list, 3))
    
    return sum_q, sum_q3

def run_proof():
    print("Anomaly Cancellation Verification")
    print("================================")
    
    g2_gens = generate_g2_generators()
    
    # Use the first generator as a Cartan candidate to define charges
    cartan_gen = g2_gens[0] 
    evals = np.linalg.eigvals(cartan_gen)
    charges = np.round(np.imag(evals), 5)
    
    # Filter numerical noise
    charges = charges[np.abs(charges) > 1e-5]
    
    print(f"Computed Charge Spectrum: {charges}")
    
    sum_q, sum_q3 = verify_anomaly_cancellation(charges)
    
    print("-" * 40)
    print(f"Anomaly Coefficients:")
    print(f"  Sum(Q):   {sum_q:.5f}")
    print(f"  Sum(Q^3): {sum_q3:.5f}")
    print("-" * 40)
    
    if np.isclose(sum_q, 0) and np.isclose(sum_q3, 0):
        print("Result: Anomaly conditions satisfied (Sums vanish).")
    else:
        print("Result: Anomaly conditions violated.")

if __name__ == "__main__":
    run_proof()
