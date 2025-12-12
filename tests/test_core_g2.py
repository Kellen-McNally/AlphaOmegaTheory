"""
Tests for G₂ module: Casimir operators, Clebsch-Gordan coefficients.
"""

import pytest
import numpy as np
from core import casimir
from core import g2_clebsch_gordan as clebsch_gordan


class TestCasimirOperators:
    """Test G₂ Casimir operators."""

    def test_second_casimir_value(self, tolerance):
        """Test C₂(G₂) = 4."""
        C2 = casimir.casimir_c2()
        assert abs(C2 - 4.0) < tolerance

    def test_third_casimir_value(self, tolerance):
        """Test C₃(G₂) = 11."""
        C3 = casimir.casimir_c3()
        assert abs(C3 - 11.0) < tolerance

    def test_triality_trace_value(self, tolerance):
        """Test Tr(τ³) = 17."""
        Tr_tau3 = casimir.triality_trace()
        assert abs(Tr_tau3 - 17.0) < tolerance

    def test_dark_energy_density(self, tolerance):
        """Test Ω_Λ = (Tr_τ³ - 6)/(dim + rank) = 11/16."""
        omega_lambda = casimir.dark_energy_density()
        expected = 11.0 / 16.0
        assert abs(omega_lambda - expected) < tolerance

    def test_weak_mixing_angle(self, tolerance):
        """Test sin²θ_W = 3/13."""
        sin2_theta_w = casimir.weak_mixing_angle()
        expected = 3.0 / 13.0
        assert abs(sin2_theta_w - expected) < tolerance


class TestClebschGordan:
    """Test G₂ Clebsch-Gordan coefficients."""

    def test_cabibbo_angle_range(self):
        """Test Cabibbo angle is in physical range."""
        theta_c = clebsch_gordan.cabibbo_angle()

        # Should be between 0 and π/2
        assert 0 < theta_c < np.pi / 2

        # Should be around 13 degrees
        theta_c_deg = np.degrees(theta_c)
        assert 12 < theta_c_deg < 14

    def test_cabibbo_formula(self):
        """Test θ_C is in physical range."""
        theta_c = clebsch_gordan.cabibbo_angle()
        expected = np.arcsin(np.sqrt(1.0 / 14.0))

        # Should be close to expected (within 20%)
        rel_error = abs(theta_c - expected) / expected
        assert rel_error < 0.2

    def test_ckm_v_us_computed_from_geometry(self, tolerance):
        """Test V_us = sin θ_C is computed from G₂."""
        ckm = clebsch_gordan.ckm_predictions()
        theta_c = ckm["theta_c"]
        V_us = ckm["V_us"]

        # V_us should equal sin θ_C
        expected = np.sin(theta_c)
        assert abs(V_us - expected) < tolerance

    def test_ckm_elements_in_range(self):
        """Test CKM elements are in physical range [0, 1]."""
        ckm = clebsch_gordan.ckm_predictions()

        assert 0 <= ckm["V_us"] <= 1
        assert 0 <= ckm["V_cb"] <= 1
        assert 0 <= ckm["V_ub"] <= 1

    def test_yukawa_ratios_physical(self):
        """Test Yukawa ratios are positive."""
        yukawa = clebsch_gordan.yukawa_ratios()

        assert yukawa["M_R_ratio_23"] > 0
        assert yukawa["Y_nu_3"] > 0
        assert yukawa["b_tau_unification"] > 0

    def test_yukawa_mr_ratio_7_over_8(self, tolerance):
        """Test M_R,2/M_R,3 = 7/8."""
        yukawa = clebsch_gordan.yukawa_ratios()
        ratio = yukawa["M_R_ratio_23"]
        expected = 7.0 / 8.0
        assert abs(ratio - expected) < tolerance

    def test_yukawa_y_nu3_13_over_11(self, tolerance):
        """Test Y_ν,3 = 13/11."""
        yukawa = clebsch_gordan.yukawa_ratios()
        Y_nu3 = yukawa["Y_nu_3"]
        expected = 13.0 / 11.0
        assert abs(Y_nu3 - expected) < tolerance
