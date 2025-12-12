"""
Yukawa Derivation: The Charge-Mass Connection

Models the dependence of fermion mass on Electric Charge (Q) via the Higgs 
mechanism. Defines the Charge Operator Q from the G2 Cartan subalgebra and 
computes the charge-mass splitting.
"""

import numpy as np
from core.constants import ALPHA_GUT
from physics.particle_yukawa_sedenion import SED_TABLE
from physics.particle_group_embedding import generate_g2_generators

def get_charge_eigenvalues():
    """
    Computes and returns the unique magnitudes of the charge eigenvalues from a G2 Cartan generator.
    """
    Q = get_charge_operator()
    vals, vecs = np.linalg.eig(Q[0:8, 0:8])
    eigen_charges = np.sort(np.abs(np.imag(vals)))
    unique_charges = sorted(list(set(np.round(eigen_charges, 4))))
    return [q for q in unique_charges if q > 1e-5]

def compute_charged_mass_splitting():
    print("Charge-Mass Interaction Energy")
    print("==============================")
    
    Q = get_charge_operator()
    
    charge_masses = []
    for k in range(16):
        v = np.zeros(16)
        v[k] = 1.0
        Qv = Q @ v
        m_Q = np.linalg.norm(Qv)**2
        charge_masses.append(m_Q)
        
    print("Charge Interaction Energies (Basis):")
    print(np.round(charge_masses, 4))
    
    cycle_A = [1, 2, 4] 
    cycle_B = [3, 6, 5] 
    
    mass_A = [charge_masses[i] for i in cycle_A]
    mass_B = [charge_masses[i] for i in cycle_B]
    
    print(f"Cycle A Charges: {np.round(mass_A, 4)}")
    print(f"Cycle B Charges: {np.round(mass_B, 4)}")
    
    unique_charge_magnitudes = get_charge_eigenvalues()
    
    print(f"\nUnique Charge Magnitudes: {unique_charge_magnitudes}")
    
    if len(unique_charge_magnitudes) >= 3:
        m1, m2, m3 = unique_charge_magnitudes[0], unique_charge_magnitudes[1], unique_charge_magnitudes[2]
        
        s1_val = np.exp(1/m1)
        s2_val = np.exp(1/m2)
        s3_val = np.exp(1/m3)
        print(f"\nExponential Scaling Factors (exp(1/Q)):")
        print(f"  s1 (e): {s1_val:.2e}")
        print(f"  s2 (u): {s2_val:.2e}")
        print(f"  s3 (t): {s3_val:.2e}")

if __name__ == "__main__":
    compute_charged_mass_splitting()