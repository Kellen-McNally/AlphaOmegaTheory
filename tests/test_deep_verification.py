"""
Deep Verification Suite for αΩ Framework

This test suite performs "Double-Entry Bookkeeping":
1. It re-implements core formulas LOCALLY (independent of the physics modules).
2. It compares these independent calculations against the physics module outputs.
3. It validates that the generated paper LaTeX contains the correct values.

This ensures that:
- The physics modules are calculating exactly what we think they are.
- The paper accurately reflects the latest code.
"""

import pytest
import re
import os
import numpy as np
from pathlib import Path

# Import the modules to test
from core.constants import ALPHA_GUT, DIM_G2, RANK_G2, CASIMIR_C3_G2, TRIALITY, CASIMIR_C2_G2, V_HIGGS_MEASURED
from core.g2_clebsch_gordan import wolfenstein_rho
from physics.particle_quark_yukawas import heavy_quark_yukawas, light_quark_yukawas
from physics.particle_fermion_masses import fermion_mass_summary

class TestIndependentFormulas:
    """
    Re-calculates values using raw formulas to verify implementation correctness.
    """
    
    def test_down_quark_formula_integrity(self):
        """Verify Down Quark formula: Y_d = α³ * [(τ+rank)/dim] * [C2 / (1+α)]"""
        # 1. Get code value
        code_y_d = light_quark_yukawas()['down']['Y_d']
        
        # 2. Independent calculation
        # Re-implement the formula locally
        alpha = ALPHA_GUT
        term1 = alpha**3
        term2 = (TRIALITY + RANK_G2) / DIM_G2 # 5/14
        term3 = CASIMIR_C2_G2 / (1 + alpha)   # 4 / (1 + 1/42)
        
        independent_y_d = term1 * term2 * term3
        
        # 3. Assert equality (Double-entry check)
        assert independent_y_d == pytest.approx(code_y_d, rel=1e-9), \
            f"Down Yukawa Mismatch! Code: {code_y_d}, Independent: {independent_y_d}"

    def test_up_quark_formula_integrity(self):
        """Verify Up Quark formula: Y_u = α³ * [(C3-rank)/dim] * (1 + α/rank)"""
        # 1. Get code value
        code_y_u = light_quark_yukawas()['up']['Y_u']
        
        # 2. Independent calculation
        alpha = ALPHA_GUT
        term1 = alpha**3
        term2 = (CASIMIR_C3_G2 - RANK_G2) / DIM_G2 # 9/14
        term3 = (1 + alpha / RANK_G2)
        
        independent_y_u = term1 * term2 * term3
        
        # 3. Assert equality
        assert independent_y_u == pytest.approx(code_y_u, rel=1e-9), \
            f"Up Yukawa Mismatch! Code: {code_y_u}, Independent: {independent_y_u}"

    def test_bottom_quark_formula_integrity(self):
        """Verify Bottom Quark formula: Y_b = α * [(dim²-rank)/(2dim²)] * (1-α)"""
        # 1. Get code value
        code_y_b = heavy_quark_yukawas()['bottom']['Y_b_GUT']
        
        # 2. Independent calculation
        alpha = ALPHA_GUT
        dim_sq = DIM_G2**2
        term1 = alpha
        term2 = (dim_sq - RANK_G2) / (2 * dim_sq) # 97/196
        term3 = (1 - alpha)
        
        independent_y_b = term1 * term2 * term3
        
        # 3. Assert equality
        assert independent_y_b == pytest.approx(code_y_b, rel=1e-9), \
            f"Bottom Yukawa Mismatch! Code: {code_y_b}, Independent: {independent_y_b}"

    def test_wolfenstein_rho_formula(self):
        """Verify Wolfenstein Rho: ρ = 7 / (τ*dim + rank) = 7/44"""
        code_rho = wolfenstein_rho()
        independent_rho = 7.0 / (TRIALITY * DIM_G2 + RANK_G2) # 7/44 = 0.159...
        
        assert code_rho == independent_rho, f"Rho Mismatch! Code: {code_rho}, Independent: {independent_rho}"


class TestPaperConsistency:
    """
    Parses the generated LaTeX file to ensure it matches the code's output.
    This prevents 'stale paper' syndrome.
    """
    
    def test_paper_reflects_code_values(self):
        # 1. Get the authoritative values from code
        summary = fermion_mass_summary()
        m_d_code_mev = summary['masses_quarks']['m_d_GeV'] * 1000
        m_u_code_mev = summary['masses_quarks']['m_u_GeV'] * 1000
        m_b_code_gev = summary['masses_quarks']['m_b_GeV']
        
        # Format them exactly as the template engine does (2 decimal places usually)
        s_d = f"{m_d_code_mev:.2f}"
        s_u = f"{m_u_code_mev:.2f}"
        s_b = f"{m_b_code_gev:.2f}"
        
        # 2. Read the generated LaTeX file
        tex_path = Path("paper/aomega_paper_single.tex")
        if not tex_path.exists():
            pytest.skip("Paper not generated yet, skipping consistency check")
            
        tex_content = tex_path.read_text(encoding='utf-8')
        
        # 3. Scan for the values in the LaTeX table
        # We look for the specific strings that should appear in the table
        
        # Down quark check
        assert s_d in tex_content, f"Paper inconsistency! Down mass {s_d} MeV not found in LaTeX."
        
        # Up quark check
        assert s_u in tex_content, f"Paper inconsistency! Up mass {s_u} MeV not found in LaTeX."
        
        # Bottom quark check
        assert s_b in tex_content, f"Paper inconsistency! Bottom mass {s_b} GeV not found in LaTeX."
        
        # 4. Check for updated formulas in LaTeX text
        # We changed Y_d formula to include C_2. Let's verify the text description exists.
        assert r"C_2" in tex_content, "Paper missing C_2 reference in formulas."
        assert r"1+\alpha" in tex_content, "Paper missing (1+alpha) correction in formulas."

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
