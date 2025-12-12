"""
Tests for gauge physics: unification, beta functions, running.
"""

import pytest
import numpy as np
from physics import particle_gauge_unification as unification
from physics import particle_beta_functions as beta_functions
from core import constants


class TestGaugeUnification:
    """Test gauge coupling unification."""

    def test_alpha_gut_value(self, tolerance):
        """Test α_GUT = 1/42."""
        alpha_gut = unification.alpha_gut_prediction()
        expected = 1.0 / 42.0
        assert abs(alpha_gut - expected) < tolerance

    def test_alpha_gut_from_geometry(self, g2_constants, tolerance):
        """Test α_GUT = 1/(triality × dim)."""
        alpha_gut = unification.alpha_gut_prediction()
        expected = 1.0 / (g2_constants["triality"] * g2_constants["dim"])
        assert abs(alpha_gut - expected) < tolerance

    def test_run_to_gut_scale(self):
        """Test RG evolution to GUT scale."""
        # Initial couplings at M_Z
        alpha_mz = np.array([1/59.0, 1/29.6, 0.12])

        M_GUT, alpha_gut, mu_array, alpha_trajectory = unification.run_to_gut_scale(
            alpha_mz, m_z=91.2
        )

        # Check GUT scale is positive and large
        assert M_GUT > 1e10  # > 10^10 GeV

        # Check couplings array returned
        assert len(alpha_gut) == 3
        assert np.all(alpha_gut > 0)

        # Check trajectory returned
        assert len(mu_array) > 0
        assert alpha_trajectory.shape[1] == 3

    def test_unification_quality(self):
        """Test couplings approximately unify at GUT scale."""
        alpha_mz = np.array([1/59.0, 1/29.6, 0.12])

        M_GUT, alpha_gut, _, _ = unification.run_to_gut_scale(alpha_mz)

        # Check couplings are close to each other
        spread = np.max(alpha_gut) - np.min(alpha_gut)
        mean_coupling = np.mean(alpha_gut)

        # Relative spread should be small
        relative_spread = spread / mean_coupling
        assert relative_spread < 0.5  # Within 50%


class TestBetaFunctions:
    """Test RG beta functions."""

    def test_beta_functions_callable(self):
        """Test beta functions can be computed."""
        alpha = np.array([0.01, 0.02, 0.12])
        b_coeffs = beta_functions.beta_coefficients_sm()

        # b_coeffs is a tuple, convert to array
        b = np.array(b_coeffs)

        beta = beta_functions.beta_1loop(alpha, b)

        assert len(beta) == 3
        assert np.all(np.isfinite(beta))

    def test_rg_evolve_forward(self):
        """Test RG evolution runs forward in energy."""
        alpha_initial = np.array([1/59.0, 1/29.6, 0.12])
        t_initial = np.log(91.2)
        t_final = np.log(1e10)

        t_array, alpha_trajectory = beta_functions.rg_evolve(
            alpha_initial, t_initial, t_final, n_steps=100
        )

        # Check evolution completed
        assert len(t_array) == 100
        assert alpha_trajectory.shape == (100, 3)

        # All couplings should remain positive
        assert np.all(alpha_trajectory > 0)

    def test_beta_coefficients_exist(self):
        """Test SM beta function coefficients are defined."""
        b = beta_functions.beta_coefficients_sm()

        # Returns tuple of 3 coefficients
        assert len(b) == 3

        # All should be finite numbers
        assert all(np.isfinite(x) for x in b)
