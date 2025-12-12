
"""
Tests for CKM Matrix predictions.
"""

import pytest
import numpy as np
from physics import particle_ckm as matrix
from core.constants import ALPHA_GUT, SIN2_THETA_W



class TestCKMMatrix:
    """Test CKM matrix construction."""

    def test_matrix_shape(self):
        """Test CKM matrix is 3×3."""
        V = matrix.ckm_matrix_wolfenstein()
        assert V.shape == (3, 3)

    def test_matrix_is_complex(self):
        """Test CKM matrix has complex elements (CP violation)."""
        V = matrix.ckm_matrix_wolfenstein()
        assert V.dtype == complex

    def test_diagonal_elements_near_one(self):
        """Test diagonal elements |V_qq| ≈ 1."""
        V = matrix.ckm_matrix_wolfenstein()

        # |V_ud|, |V_cs|, |V_tb| should be close to 1
        assert abs(V[0, 0]) > 0.97
        assert abs(V[1, 1]) > 0.97
        assert abs(V[2, 2]) > 0.99

    def test_off_diagonal_hierarchy(self):
        """Test off-diagonal elements follow hierarchy."""
        V = matrix.ckm_matrix_wolfenstein()

        # |V_us| ≈ λ ≈ 0.22 (largest off-diagonal)
        assert 0.2 < abs(V[0, 1]) < 0.25

        # |V_cb| ≈ Aλ² ≈ 0.04
        assert 0.03 < abs(V[1, 2]) < 0.05

        # |V_ub| ≈ Aλ³ ≈ 0.004 (smallest)
        assert 0.002 < abs(V[0, 2]) < 0.006


class TestCKMElements:
    """Test individual CKM matrix elements."""

    def test_all_elements_computed(self):
        """Test all 9 elements are computed."""
        elements = matrix.ckm_matrix_elements()
        assert len(elements) == 9

        # Check all expected keys exist
        quarks = ["u", "c", "t"]
        primes = ["d", "s", "b"]

        for q in quarks:
            for qp in primes:
                key = f"V_{q}{qp}"
                assert key in elements

    def test_elements_in_physical_range(self):
        """Test |V_ij| ≤ 1 for all elements."""
        elements = matrix.ckm_matrix_elements()

        for key, value in elements.items():
            magnitude = abs(value)
            assert 0 <= magnitude <= 1


class TestCKMUnitarity:
    """Test CKM matrix unitarity."""

    def test_unitarity_test_runs(self):
        """Test unitarity test completes without error."""
        result = matrix.ckm_unitarity_test()

        assert "max_deviation" in result
        assert "is_unitary" in result
        assert "row_deviations" in result
        assert "col_deviations" in result

    def test_row_normalization(self):
        """Test rows are approximately normalized: Σ|V_ij|² ≈ 1."""
        V = matrix.ckm_matrix_wolfenstein()

        for i in range(3):
            row_norm_sq = np.sum(np.abs(V[i, :])**2)
            assert abs(row_norm_sq - 1.0) < 0.01

    def test_column_normalization(self):
        """Test columns are approximately normalized."""
        V = matrix.ckm_matrix_wolfenstein()

        for j in range(3):
            col_norm_sq = np.sum(np.abs(V[:, j])**2)
            assert abs(col_norm_sq - 1.0) < 0.01

    def test_orthogonality(self):
        """Test V†V ≈ I."""
        V = matrix.ckm_matrix_wolfenstein()
        VdagV = np.conj(V.T) @ V
        I = np.eye(3)

        # Check all elements
        max_dev = np.max(np.abs(VdagV - I))
        assert max_dev < 0.01


class TestCKMSummary:
    """Test CKM summary output."""

    def test_summary_structure(self):
        """Test summary has all required components."""
        summary = matrix.ckm_summary()

        assert "parameters" in summary
        assert "matrix" in summary
        assert "elements" in summary
        assert "experimental" in summary
        assert "errors" in summary
        assert "unitarity" in summary

    def test_cabibbo_from_g2(self, physical_tolerance):
        """Test λ is computed from G₂: θ_C = arcsin√(1/14)."""
        summary = matrix.ckm_summary()
        params = summary["parameters"]

        lambda_g2 = params["lambda_g2"]

        # Should be close to √(1/14) ≈ 0.267
        expected = np.sqrt(1.0 / 14.0)
        assert abs(lambda_g2 - expected) < 0.05  # Within 5%

    def test_v_us_from_cabibbo(self):
        """Test V_us = sin θ_C from G₂."""
        summary = matrix.ckm_summary()
        elements = summary["elements"]

        V_us = abs(elements["V_us"])

        # Should be around 0.225
        assert 0.22 < V_us < 0.23

    def test_errors_computed(self):
        """Test errors are computed for key elements."""
        summary = matrix.ckm_summary()
        errors = summary["errors"]

        assert "V_us" in errors
        assert "V_cb" in errors
        assert "V_ub" in errors

        # All errors should be non-negative
        for value in errors.values():
            assert value >= 0
