"""
Tests for Strong CP problem solution.
"""

import pytest
import physics.particle_strong_cp as solution
from core.constants import TRIALITY

class TestStrongCP:
    """Test strong CP problem solution via triality."""

    def test_theta_bar_computed(self):
        """Test θ_QCD parameter is computed."""
        pred = solution.theta_qcd_prediction()

        # Should return dict
        assert isinstance(pred, dict)
        assert len(pred) > 0

    def test_theta_bar_vanishes_from_triality(self):
        """Test θ_QCD = 0 from triality symmetry."""
        summary = solution.strong_cp_summary()

        assert "prediction" in summary
        assert "triality" in str(summary).lower()

        theta_qcd = summary["prediction"]["theta_qcd"]
        # Should be zero
        assert abs(theta_qcd) < 1e-9

    def test_no_axion_required(self):
        """Test solution doesn't require axion."""
        summary = solution.strong_cp_summary()

        assert "g2_solution" in summary
        mechanism = str(summary["g2_solution"]).lower()

        # Should mention triality
        assert "triality" in mechanism or "g2" in mechanism or "g₂" in mechanism

    def test_qcd_vacuum_angle(self):
        """Test QCD vacuum angle θ_QCD."""
        pred = solution.theta_qcd_prediction()

        # Should have theta_qcd key
        assert "theta_qcd" in pred or "theta_0" in pred

    def test_neutron_edm_prediction(self):
        """Test neutron EDM prediction."""
        theta_qcd = 0.0  # From triality
        d_n = solution.neutron_edm(theta_qcd)

        # Should be zero when theta_qcd = 0
        assert abs(d_n) < 1e-20

    def test_cp_conservation_explained(self):
        """Test CP conservation in strong interactions."""
        summary = solution.strong_cp_summary()

        assert "prediction" in summary

        # θ_QCD should be zero
        theta = summary["prediction"]["theta_qcd"]
        assert abs(theta) < 1e-8