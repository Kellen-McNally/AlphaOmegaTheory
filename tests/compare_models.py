import sys
import os
import numpy as np

# Add repo root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from core.constants import (
    M_GUT, V_HIGGS_MEASURED, DIM_G2, TRIALITY
)
from physics.particle_neutrino_masses import mass_splittings_rigorous, DM21_SQ_EXP, DM31_SQ_EXP
from physics.particle_seesaw import rh_neutrino_masses, seesaw_summary
from physics.atomic_ionization import calculate_ionization_energy, calculate_ionization_energy_octonion, IE_EXPERIMENTAL
from physics.particle_higgs_vev import calculate_geometric_vev

def compare_neutrinos():
    print("\n--- NEUTRINO MASS SPLITTING COMPARISON ---")
    
    # 1. Current Model (Rigorous / 9/7 factor / Experimental VEV)
    current = mass_splittings_rigorous()
    curr_dm21 = current['prediction']['Dm21_sq_eV2']
    curr_dm31 = current['prediction']['Dm31_sq_eV2']
    curr_err21 = abs(curr_dm21 - DM21_SQ_EXP)/DM21_SQ_EXP * 100
    curr_err31 = abs(curr_dm31 - DM31_SQ_EXP)/DM31_SQ_EXP * 100
    
    print(f"Current Model (9/7, V_exp):")
    print(f"  Dm21 Error: {curr_err21:.2f}%")
    print(f"  Dm31 Error: {curr_err31:.2f}%")
    
    # 2. Alternative Model (Seesaw / 7/6 factor / Experimental VEV)
    # Re-implement logic from particle_seesaw.py but with flexibility
    # The seesaw module calculates masses but we need to check if it matches experimental VEV usage
    
    # From particle_seesaw.py:
    # M_R3 = M_GUT * (7/6)
    # M_R2 = M_R3 * (7/8)
    # M_R1 = M_R3 * (1/20)
    # m_i = v^2 * Y_i^2 / M_Ri
    
    # Geometric Yukawas (from particle_neutrino_masses.py logic)
    # Y3 = 13/11
    # Y2 = sqrt(13/11 * 1/20) ... wait, let's use the exact ones from the file
    
    # Let's just run seesaw_summary() which uses its internal logic
    alt = seesaw_summary()
    alt_dm21 = alt['comparison']['Dm21_sq_predicted']
    alt_dm31 = alt['comparison']['Dm31_sq_predicted']
    alt_err21 = alt['comparison']['Dm21_error_percent']
    alt_err31 = alt['comparison']['Dm31_error_percent']
    
    print(f"Alternative Model (7/6, V_exp):")
    print(f"  Dm21 Error: {alt_err21:.2f}%")
    print(f"  Dm31 Error: {alt_err31:.2f}%")
    
    # 3. Refined Model (7/6 factor / Calculated VEV)
    # Recalculate using V_pred
    v_pred = calculate_geometric_vev()
    ratio_v = (v_pred / V_HIGGS_MEASURED)**2 # Masses scale as v^2, so m_nu scales as v^2 (since M_R is fixed)
    # Actually m_nu ~ v^2. So Dm^2 ~ v^4.
    
    ref_dm21 = alt_dm21 * (ratio_v**2) # Scale Dm^2
    ref_dm31 = alt_dm31 * (ratio_v**2)
    
    ref_err21 = abs(ref_dm21 - DM21_SQ_EXP)/DM21_SQ_EXP * 100
    ref_err31 = abs(ref_dm31 - DM31_SQ_EXP)/DM31_SQ_EXP * 100
    
    print(f"Refined Model (7/6, V_pred={v_pred:.2f}):")
    print(f"  Dm21 Error: {ref_err21:.2f}%")
    print(f"  Dm31 Error: {ref_err31:.2f}%")


def compare_atomic():
    print("\n--- ATOMIC IONIZATION COMPARISON ---")
    
    z_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 18]
    
    # Current Model
    curr_errors = []
    for z in z_list:
        pred = calculate_ionization_energy(z)
        exp = IE_EXPERIMENTAL[z]
        curr_errors.append(abs(pred - exp)/exp * 100)
    curr_avg = sum(curr_errors) / len(curr_errors)
    
    print(f"Current Model (Standard):")
    print(f"  Mean Error: {curr_avg:.2f}%")
    
    # Alternative Model
    alt_errors = []
    for z in z_list:
        # Check if function exists/works
        try:
            pred = calculate_ionization_energy_octonion(z)
            exp = IE_EXPERIMENTAL[z]
            alt_errors.append(abs(pred - exp)/exp * 100)
        except Exception as e:
            print(f"  Failed for Z={z}: {e}")
            
    if alt_errors:
        alt_avg = sum(alt_errors) / len(alt_errors)
        print(f"Alternative Model (Advanced Octonion):")
        print(f"  Mean Error: {alt_avg:.2f}%")

if __name__ == "__main__":
    compare_neutrinos()
    compare_atomic()
