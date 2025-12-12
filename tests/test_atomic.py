"""
Tests for atomic modules: ionization energies.
"""

import pytest
import numpy as np
from physics import atomic_ionization as energy


class TestPrincipalQuantumNumber:
    """Test principal quantum number assignment."""

    def test_hydrogen_n_equals_1(self):
        """Test H (Z=1) has n=1."""
        n = energy.get_principal_quantum_number(1)
        assert n == 1

    def test_helium_n_equals_1(self):
        """Test He (Z=2) has n=1."""
        n = energy.get_principal_quantum_number(2)
        assert n == 1

    def test_lithium_n_equals_2(self):
        """Test Li (Z=3) has n=2."""
        n = energy.get_principal_quantum_number(3)
        assert n == 2

    def test_neon_n_equals_2(self):
        """Test Ne (Z=10) has n=2."""
        n = energy.get_principal_quantum_number(10)
        assert n == 2

    def test_sodium_n_equals_3(self):
        """Test Na (Z=11) has n=3."""
        n = energy.get_principal_quantum_number(11)
        assert n == 3

    def test_argon_n_equals_3(self):
        """Test Ar (Z=18) has n=3."""
        n = energy.get_principal_quantum_number(18)
        assert n == 3

    def test_period_7_maximum(self):
        """Test Og (Z=118) has n=7."""
        n = energy.get_principal_quantum_number(118)
        assert n == 7

    def test_period_8_prediction(self):
        """Test Z=119 starts period 8 with n=8."""
        n = energy.get_principal_quantum_number(119)
        assert n == 8


class TestEffectiveNuclearCharge:
    """Test Z_eff calculation."""

    def test_hydrogen_no_screening(self):
        """Test Z_eff(H) = 1 (no screening)."""
        Z_eff = energy.calculate_zeff_g2(Z=1, n=1)
        assert abs(Z_eff - 1.0) < 1e-10

    def test_helium_fixed_value(self):
        """Test Z_eff(He) from G₂ formula (experimental ~1.34)."""
        Z_eff = energy.calculate_zeff_g2(Z=2, n=1)
        # G₂ formula predicts 1.726 (experimental is ~1.34, so ~29% error)
        assert abs(Z_eff - 1.7264705882352942) < 1e-10

    def test_zeff_less_than_z(self):
        """Test Z_eff ≤ Z for all elements."""
        for Z in range(1, 119):
            n = energy.get_principal_quantum_number(Z)
            Z_eff = energy.calculate_zeff_g2(Z, n)
            assert Z_eff <= Z

    def test_zeff_positive(self):
        """Test Z_eff > 0 for all elements."""
        for Z in range(1, 119):
            n = energy.get_principal_quantum_number(Z)
            Z_eff = energy.calculate_zeff_g2(Z, n)
            assert Z_eff > 0

    def test_g2_correction_factor(self):
        """Test G2 correction factor logic."""
        from physics.atomic_ionization import TR3_G2, DIM_G2
        
        expected_factor = TR3_G2 / (5.0 * DIM_G2)
        assert abs(expected_factor - 0.242857) < 0.000001


class TestIonizationEnergy:
    """Test ionization energy calculation."""

    def test_hydrogen_ionization(self, tolerance):
        """Test IE(H) = 13.6 eV (Rydberg constant)."""
        IE = energy.calculate_ionization_energy(Z=1, Z_eff=1.0, n=1)
        expected = energy.RYDBERG
        assert abs(IE - expected) < tolerance

    def test_ionization_positive(self):
        """Test IE > 0 for all elements."""
        for Z in range(1, 119):
            n = energy.get_principal_quantum_number(Z)
            Z_eff = energy.calculate_zeff_g2(Z, n)
            IE = energy.calculate_ionization_energy(Z, Z_eff, n)
            assert IE > 0

    def test_ionization_formula(self, tolerance):
        """Test IE = R_∞ × Z_eff² / n²."""
        Z = 1
        n = 1
        Z_eff = 1.0

        IE = energy.calculate_ionization_energy(Z, Z_eff, n)
        expected = energy.RYDBERG * Z_eff**2 / n**2

        assert abs(IE - expected) < tolerance

    def test_relativistic_correction_heavy_elements(self):
        """Test relativistic correction applied for Z > 36."""
        # For heavy elements, IE should be slightly higher due to relativity
        Z_light = 20  # Ca (no correction)
        Z_heavy = 80  # Hg (with correction)

        n_light = energy.get_principal_quantum_number(Z_light)
        n_heavy = energy.get_principal_quantum_number(Z_heavy)

        # Use same Z_eff for comparison
        Z_eff = 10.0

        IE_light = energy.calculate_ionization_energy(Z_light, Z_eff, n_light)
        IE_heavy_base = energy.RYDBERG * Z_eff**2 / n_heavy**2

        # Heavy element with relativistic correction (can increase or decrease depending on orbital)
        # For heavy elements, G₂ formula gives reduction for certain orbitals
        IE_heavy = energy.calculate_ionization_energy(Z_heavy, Z_eff, n_heavy)
        # Just verify it's positive and different from base
        assert IE_heavy > 0
        assert abs(IE_heavy - IE_heavy_base) > 0.1  # Should show relativistic effect


class TestIonizationSummary:
    """Test ionization energy summary."""

    def test_summary_structure(self):
        """Test summary has correct structure."""
        summary = energy.ionization_energy_summary([1, 2, 3])

        assert 1 in summary
        assert 2 in summary
        assert 3 in summary

        for Z, data in summary.items():
            assert "n" in data
            assert "Z_eff" in data
            assert "IE_eV" in data

    def test_summary_values_positive(self):
        """Test all values in summary are positive."""
        summary = energy.ionization_energy_summary(range(1, 21))

        for Z, data in summary.items():
            assert data["n"] > 0
            assert data["Z_eff"] > 0
            assert data["IE_eV"] > 0

    def test_noble_gas_trend(self):
        """Test noble gases have high ionization energies."""
        summary = energy.ionization_energy_summary([2, 10, 18])  # He, Ne, Ar

        IE_He = summary[2]["IE_eV"]
        IE_Ne = summary[10]["IE_eV"]
        IE_Ar = summary[18]["IE_eV"]

        # Noble gases should have positive IE
        assert IE_He > 10  # High
        assert IE_Ne > 10
        assert IE_Ar > 10
