"""
Fermion Mass Spectrum: Quantitative Prediction

Derives charged lepton masses (electron, muon, tau) using geometrically 
derived Yukawa couplings and the calculated Higgs VEV.
"""

import numpy as np
import sys
import os

# Ensure we can import from core
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from core.constants import (
    M_ELECTRON_GEV_MEASURED, M_MUON_GEV_MEASURED, M_TAU_GEV_MEASURED,
    DIM_G2, TRIALITY, M_HIGGS_MEASURED
)
from physics.particle_fermion_masses import charged_lepton_yukawas
from physics.particle_higgs_vev import calculate_geometric_vev

# Experimental masses (in GeV)
M_ELECTRON_GEV_EXP = M_ELECTRON_GEV_MEASURED
M_MUON_GEV_EXP = M_MUON_GEV_MEASURED
M_TAU_GEV_EXP = M_TAU_GEV_MEASURED
M_HIGGS_GEV_EXP = M_HIGGS_MEASURED

def calculate_mass_spectrum():
    print("Fermion Mass Spectrum Derivation")
    print("===============================")
    
    yukawas = charged_lepton_yukawas()
    Y_e = yukawas['Y_e']
    Y_mu = yukawas['Y_mu']
    Y_tau = yukawas['Y_tau']
    
    print(f"Geometric Yukawa Couplings:")
    print(f"  Y_e:   {Y_e:.10f}")
    print(f"  Y_mu:  {Y_mu:.10f}")
    print(f"  Y_tau: {Y_tau:.10f}")
    
    v_higgs_pred = calculate_geometric_vev()
    print(f"\nPredicted Higgs VEV (v): {v_higgs_pred:.4f} GeV")

    lambda_bare = DIM_G2 + TRIALITY 
    manifold_vol = DIM_G2 * TRIALITY 
    phase_space = np.pi
    lambda_eff = lambda_bare / (manifold_vol * phase_space)
    
    m_h_pred = np.sqrt(2 * lambda_eff) * v_higgs_pred
    m_h_exp = M_HIGGS_GEV_EXP
    err_h = abs(m_h_pred - m_h_exp) / m_h_exp * 100
    
    print(f"Predicted Higgs Mass:    {m_h_pred:.4f} GeV (Exp: {m_h_exp})")
    print(f"Higgs Error:             {err_h:.2f}%")
    
    m_e_pred = Y_e * v_higgs_pred
    m_mu_pred = Y_mu * v_higgs_pred
    m_tau_pred = Y_tau * v_higgs_pred
    
    print(f"\nPredicted Masses (GeV):")
    print(f"  Electron: {m_e_pred:.10f} (Exp: {M_ELECTRON_GEV_EXP:.10f})")
    print(f"  Muon:     {m_mu_pred:.10f} (Exp: {M_MUON_GEV_EXP:.10f})")
    print(f"  Tau:      {m_tau_pred:.10f} (Exp: {M_TAU_GEV_EXP:.10f})")
    
    err_e = abs(m_e_pred - M_ELECTRON_GEV_EXP) / M_ELECTRON_GEV_EXP * 100
    err_mu = abs(m_mu_pred - M_MUON_GEV_EXP) / M_MUON_GEV_EXP * 100
    err_tau = abs(m_tau_pred - M_TAU_GEV_EXP) / M_TAU_GEV_EXP * 100
    
    print(f"\nErrors:")
    print(f"  Electron: {err_e:.2f}%")
    print(f"  Muon:     {err_mu:.2f}%")
    print(f"  Tau:      {err_tau:.2f}%")
    
    avg_error = (err_e + err_mu + err_tau) / 3
    print(f"  Average Error: {avg_error:.2f}%")
        
    return {
        'electron_mass_gev_pred': m_e_pred,
        'muon_mass_gev_pred': m_mu_pred,
        'tau_mass_gev_pred': m_tau_pred,
        'higgs_mass_gev_pred': m_h_pred,
        'higgs_error_pct': err_h,
        'electron_error_pct': err_e,
        'muon_error_pct': err_mu,
        'tau_error_pct': err_tau,
        'average_error_pct': avg_error,
        'v_higgs_pred': v_higgs_pred,
    }

if __name__ == "__main__":
    calculate_mass_spectrum()