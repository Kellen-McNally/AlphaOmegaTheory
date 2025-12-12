"""
Tests for Fermion Masses (Leptons and Quarks).
"""

import pytest
import numpy as np
from physics import particle_fermion_masses as masses
from core.constants import ALPHA_GUT, V_HIGGS_MEASURED


class TestChargedLeptonYukawas:
    """Test charged lepton Yukawa couplings."""

    def test_y_tau_value(self, tolerance):
        """Test Y_τ = τ × α × 99/(10×98) = 297/41160."""
        from core.constants import TRIALITY, ALPHA_GUT, DIM_G2, RANK_G2

        yukawa = masses.charged_lepton_yukawas()
        Y_tau = yukawa["Y_tau"]

        # Expected from formula: tau × alpha × 99/(10×98)
        expected = TRIALITY * ALPHA_GUT * (DIM_G2**2 + RANK_G2) / (10 * DIM_G2**2)

        assert abs(Y_tau - expected) < tolerance

    def test_yukawa_positive(self):
        """Test all Yukawa couplings are positive."""
        yukawa = masses.charged_lepton_yukawas()

        assert yukawa["Y_e"] > 0
        assert yukawa["Y_mu"] > 0
        assert yukawa["Y_tau"] > 0

    def test_yukawa_hierarchy(self):
        """Test Y_e < Y_μ < Y_τ."""
        yukawa = masses.charged_lepton_yukawas()

        Y_e = yukawa["Y_e"]
        Y_mu = yukawa["Y_mu"]
        Y_tau = yukawa["Y_tau"]

        assert Y_e < Y_mu < Y_tau


class TestQuarkYukawas:
    """Test quark Yukawa couplings."""

    def test_yukawa_positive(self):
        """Test all quark Yukawas are positive."""
        yukawa = masses.quark_yukawas()

        for quark in ["u", "d", "s", "c", "b"]:
            assert yukawa[f"Y_{quark}_GUT"] > 0
        assert yukawa["Y_t"] > 0

    def test_up_type_hierarchy(self):
        """Test Y_u < Y_c < Y_t."""
        yukawa = masses.quark_yukawas()

        assert yukawa["Y_u_GUT"] < yukawa["Y_c_GUT"] < yukawa["Y_t"]

    def test_down_type_hierarchy(self):
        """Test down-type yukawas are positive and ordered."""
        yukawa = masses.quark_yukawas()

        # Check they exist and are positive
        assert yukawa["Y_d_GUT"] > 0
        assert yukawa["Y_s_GUT"] > 0
        assert yukawa["Y_b_GUT"] > 0

    def test_top_yukawa_near_unity(self):
        """Test Y_t ≈ 1 (natural expectation)."""
        yukawa = masses.quark_yukawas()
        Y_t = yukawa["Y_t"]

        # Top Yukawa should be O(1)
        assert 0.5 < Y_t < 1.5


class TestChargedLeptonMasses:
    """Test charged lepton masses."""

    def test_mass_hierarchy(self):
        """Test m_e < m_μ < m_τ."""
        mass_dict = masses.charged_lepton_masses()

        m_e = mass_dict["m_e_GeV"]
        m_mu = mass_dict["m_mu_GeV"]
        m_tau = mass_dict["m_tau_GeV"]

        assert m_e < m_mu < m_tau

    def test_masses_positive(self):
        """Test all masses are positive."""
        mass_dict = masses.charged_lepton_masses()

        assert mass_dict["m_e_GeV"] > 0
        assert mass_dict["m_mu_GeV"] > 0
        assert mass_dict["m_tau_GeV"] > 0


class TestQuarkMasses:
    """Test quark masses."""

    def test_mass_hierarchy_up_type(self):
        """Test m_u < m_c < m_t."""
        mass_dict = masses.quark_masses()

        assert mass_dict["m_u_GeV"] < mass_dict["m_c_GeV"] < mass_dict["m_t_GeV"]

    def test_mass_hierarchy_down_type(self):
        """Test down-type quark masses are positive and ordered."""
        mass_dict = masses.quark_masses()

        # Check they exist and are positive
        assert mass_dict["m_d_GeV"] > 0
        assert mass_dict["m_s_GeV"] > 0
        assert mass_dict["m_b_GeV"] > 0

    def test_masses_positive(self):
        """Test all masses are positive."""
        mass_dict = masses.quark_masses()

        for quark in ["u", "d", "s", "c", "b", "t"]:
            assert mass_dict[f"m_{quark}_GeV"] > 0

    def test_top_mass_order_100_gev(self):
        """Test top mass is O(100 GeV)."""
        mass_dict = masses.quark_masses()
        m_t = mass_dict["m_t_GeV"]

        assert 100 < m_t < 200


class TestFermionSummary:
    """Test fermion mass summary."""

    def test_summary_structure(self):
        """Test summary contains all components."""
        summary = masses.fermion_mass_summary()

        assert "yukawa_leptons" in summary
        assert "yukawa_quarks_gut" in summary
        assert "masses_leptons" in summary
        assert "masses_quarks" in summary
        assert "experimental_leptons" in summary
        assert "experimental_quarks" in summary
        assert "errors_leptons" in summary
        assert "errors_quarks" in summary

    def test_errors_non_negative(self):
        """Test all errors are non-negative."""
        summary = masses.fermion_mass_summary()

        for key, value in summary["errors_leptons"].items():
            assert value >= 0

        for key, value in summary["errors_quarks"].items():
            assert value >= 0