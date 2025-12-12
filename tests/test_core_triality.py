
"""
Tests for triality module: τ³ = 1 automorphism.
"""

import pytest
import numpy as np
from core import triality


class TestTriality:
    """Test triality automorphism properties."""


    def test_order(self):
        """Test triality has order 3: τ³ = 1."""
        order = triality.triality_order()
        assert order == 3

    def test_three_phases(self):
        """Test there are exactly 3 triality phases."""
        phases = triality.triality_phases()
        assert len(phases) == 3

    def test_phases_are_cube_roots_of_unity(self, tolerance):
        """Test ω³ = 1 for each phase."""
        phases = triality.triality_phases()

        for omega in phases:
            omega_cubed = omega ** 3
            # Should equal 1
            assert abs(omega_cubed - 1.0) < tolerance

    def test_first_phase_is_one(self, tolerance):
        """Test ω⁰ = 1."""
        phases = triality.triality_phases()
        assert abs(phases[0] - 1.0) < tolerance

    def test_phases_on_unit_circle(self, tolerance):
        """Test |ω| = 1 for all phases."""
        phases = triality.triality_phases()

        for omega in phases:
            magnitude = abs(omega)
            assert abs(magnitude - 1.0) < tolerance

    def test_phases_equally_spaced(self, tolerance):
        """Test phases are equally spaced by 2π/3."""
        phases = triality.triality_phases()

        # ω₁ = e^(i2π/3)
        expected_phase_1 = np.exp(1j * 2.0 * np.pi / 3.0)
        assert abs(phases[1] - expected_phase_1) < tolerance

        # ω₂ = e^(i4π/3)
        expected_phase_2 = np.exp(1j * 4.0 * np.pi / 3.0)
        assert abs(phases[2] - expected_phase_2) < tolerance

    def test_phases_sum_to_zero(self, tolerance):
        """Test ω⁰ + ω¹ + ω² = 0."""
        phases = triality.triality_phases()
        phase_sum = sum(phases)
        assert abs(phase_sum) < tolerance

    def test_generation_indices(self):
        """Test three generation indices (1, 2, 3)."""
        indices = triality.generation_indices()
        assert indices == (1, 2, 3)
        assert len(indices) == 3
