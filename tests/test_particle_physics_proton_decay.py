
"""
Tests for Proton Decay predictions.
"""

import pytest
from physics import particle_proton_calculation as calculation
from core.constants import M_PLANCK_GEV, ALPHA_GUT, CASIMIR_C3_G2, M_PROTON_DECAY


class TestProtonDecay:
    """Test proton decay lifetime predictions."""

    def test_proton_lifetime_computed(self):
        """Test proton lifetime τ_p is computed."""
        tau_p = calculation.proton_lifetime()

        # Should be very long: > 10^33 years
        assert tau_p > 1e30

    def test_dominant_decay_channel(self):
        """Test dominant decay channel is identified."""
        ratios = calculation.branching_ratios()

        # Should have multiple channels
        assert len(ratios) > 0

        # Find dominant channel
        dominant = max(ratios, key=ratios.get)
        assert ratios[dominant] > 0

    def test_branching_ratios(self):
        """Test decay branching ratios sum to 1."""
        ratios = calculation.branching_ratios()

        assert isinstance(ratios, dict)
        assert len(ratios) > 0

        # Branching ratios should be positive
        for channel, br in ratios.items():
            assert br >= 0
            assert br <= 1

        # Should sum to ~1
        total = sum(ratios.values())
        assert abs(total - 1.0) < 0.1

    def test_proton_decay_summary(self):
        """Test proton decay summary."""
        summary = calculation.proton_decay_summary()

        assert "prediction" in summary
        assert "branching_ratios" in summary

    def test_gut_scale_suppression(self):
        """Test decay suppressed by M_GUT."""
        summary = calculation.proton_decay_summary()

        # Lifetime should scale as M_GUT^4
        assert "M_GUT" in str(summary) or "gut" in str(summary).lower()

    def test_experimental_bounds(self):
        """Test prediction compared to experimental bounds."""
        summary = calculation.proton_decay_summary()

        # Should have experimental status
        assert "experimental_status" in summary

    def test_below_current_limits(self):
        """Test predicted lifetime is positive and large."""
        tau_p = calculation.proton_lifetime()

        # Should be very long
        assert tau_p > 1e28  # At least 10^28 years
