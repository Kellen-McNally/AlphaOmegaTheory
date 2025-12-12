"""
Tests for quantum corrections module.
"""

import pytest
import numpy as np
from physics import particle_quantum_corrections as quantum_corrections


class TestAlphaS:
    """Test running strong coupling constant."""

    def test_alpha_s_at_mz(self):
        """Test α_s(M_Z) returns reference value."""
        alpha = quantum_corrections.alpha_s(quantum_corrections.M_Z)
        # Should return the reference value
        assert abs(alpha - quantum_corrections.ALPHA_S_MZ) < 1e-10

    def test_alpha_s_decreases_at_high_energy(self):
        """Test α_s decreases at higher energies (asymptotic freedom)."""
        alpha_low = quantum_corrections.alpha_s(quantum_corrections.M_Z)
        alpha_high = quantum_corrections.alpha_s(1000.0)  # 1 TeV
        assert alpha_high < alpha_low

    def test_alpha_s_positive(self):
        """Test α_s is always positive."""
        for mu in [10.0, 91.2, 1000.0, 10000.0]:
            assert quantum_corrections.alpha_s(mu) > 0


class TestAlphaEM:
    """Test running electromagnetic coupling constant."""

    def test_alpha_em_at_mz(self):
        """Test α_EM(M_Z) returns reference value."""
        alpha = quantum_corrections.alpha_em(quantum_corrections.M_Z)
        # Should return the reference value
        assert abs(alpha - quantum_corrections.ALPHA_EM_MZ) < 1e-10

    def test_alpha_em_increases_at_high_energy(self):
        """Test α_EM increases at higher energies."""
        alpha_low = quantum_corrections.alpha_em(quantum_corrections.M_Z)
        alpha_high = quantum_corrections.alpha_em(1000.0)  # 1 TeV
        assert alpha_high > alpha_low

    def test_alpha_em_positive(self):
        """Test α_EM is always positive."""
        for mu in [0.511e-3, 91.2, 1000.0]:
            assert quantum_corrections.alpha_em(mu) > 0


class TestOneLoopCorrections:
    """Test 1-loop quantum corrections."""

    def test_one_loop_lepton_correction(self):
        """Test 1-loop correction for leptons."""
        m_e = 0.511e-3  # GeV
        delta = quantum_corrections.one_loop_correction_lepton(m_e, 'electron')

        # Should be small (allow up to 2x the mass for large logs)
        assert abs(delta) < 2.0 * m_e

    def test_one_loop_quark_correction(self):
        """Test 1-loop correction for quarks."""
        m_c = 1.27  # Charm mass in GeV
        delta = quantum_corrections.one_loop_correction_quark(m_c, 'charm')

        # Should return a value
        assert delta is not None


class TestAnalyzeRemainingError:
    """Test error analysis function."""

    def test_analyze_runs_without_error(self):
        """Test analyze_remaining_error runs without error."""
        # This function prints output but returns None
        result = quantum_corrections.analyze_remaining_error()
        # Just verify it runs without crashing
        assert True

    def test_analyze_is_callable(self):
        """Test analyze_remaining_error is callable."""
        assert callable(quantum_corrections.analyze_remaining_error)


class TestEstimateQuantumCorrections:
    """Test quantum correction estimation."""

    def test_estimate_runs_without_error(self):
        """Test estimate_quantum_corrections runs without error."""
        # This function prints output but returns None
        result = quantum_corrections.estimate_quantum_corrections()
        # Just verify it runs without crashing
        assert True

    def test_estimate_is_callable(self):
        """Test estimate_quantum_corrections is callable."""
        assert callable(quantum_corrections.estimate_quantum_corrections)


class TestModuleStructure:
    """Test module structure and constants."""

    def test_module_has_constants(self):
        """Test module defines necessary constants."""
        assert hasattr(quantum_corrections, 'M_Z')
        assert hasattr(quantum_corrections, 'ALPHA_S_MZ')
        assert hasattr(quantum_corrections, 'ALPHA_EM_MZ')

    def test_constants_physical(self):
        """Test constants have physical values."""
        # Z mass around 91 GeV
        assert 90 < quantum_corrections.M_Z < 92

        # α_s(M_Z) around 0.118
        assert 0.11 < quantum_corrections.ALPHA_S_MZ < 0.13

        # α_EM(M_Z) around 1/128
        assert 0.007 < quantum_corrections.ALPHA_EM_MZ < 0.008

    def test_module_functions_callable(self):
        """Test all main functions are callable."""
        assert callable(quantum_corrections.alpha_s)
        assert callable(quantum_corrections.alpha_em)
        assert callable(quantum_corrections.one_loop_correction_lepton)
        assert callable(quantum_corrections.one_loop_correction_quark)
        assert callable(quantum_corrections.analyze_remaining_error)
        assert callable(quantum_corrections.estimate_quantum_corrections)
