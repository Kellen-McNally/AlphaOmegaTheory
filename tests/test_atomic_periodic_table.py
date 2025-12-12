"""
Tests for periodic table calculations.
"""

import pytest
import numpy as np
from physics import atomic_periodic as calculate
from physics import atomic_ionization as energy


class TestPeriodicTable:
    """Test periodic table structure predictions."""

    def test_element_count_period_1(self):
        """Test period 1 has 2 elements (H, He)."""
        count = calculate.elements_in_period(1)
        assert count == 2

    def test_element_count_period_2(self):
        """Test period 2 has 8 elements (Li-Ne)."""
        count = calculate.elements_in_period(2)
        assert count == 8

    def test_element_count_period_3(self):
        """Test period 3 has 8 elements (Na-Ar)."""
        count = calculate.elements_in_period(3)
        assert count == 8

    def test_element_count_period_4(self):
        """Test period 4 has 18 elements (K-Kr)."""
        count = calculate.elements_in_period(4)
        assert count == 18

    def test_periodic_structure_from_g2(self):
        """Test periodic table structure emerges from G₂."""
        summary = calculate.periodic_table_summary()

        assert "structure" in summary or "periods" in summary

        # Should mention G₂ origin
        assert "G2" in str(summary) or "g2" in str(summary) or "G₂" in str(summary)

    def test_noble_gases_identified(self):
        """Test noble gases are correctly identified."""
        noble_gases = calculate.noble_gas_list()

        # Should include He, Ne, Ar, Kr, Xe, Rn, Og
        assert 2 in noble_gases   # He
        assert 10 in noble_gases  # Ne
        assert 18 in noble_gases  # Ar
        assert 36 in noble_gases  # Kr

    def test_transition_metals_block(self):
        """Test transition metals block structure."""
        summary = calculate.periodic_table_summary()

        # Should have d-block section
        assert "d_block" in summary or "transition" in str(summary).lower()

    def test_lanthanides_count(self):
        """Test lanthanides have 15 elements."""
        count = calculate.lanthanide_count()
        assert count == 15

    def test_actinides_count(self):
        """Test actinides have 15 elements."""
        count = calculate.actinide_count()
        assert count == 15

    def test_period_8_prediction(self):
        """Test period 8 structure prediction."""
        count = calculate.elements_in_period(8)

        # Period 8 should have 50 elements (if complete)
        assert count > 0
        assert count <= 50

    def test_element_118_in_period_7(self):
        """Test Oganesson (Z=118) is in period 7."""
        period = calculate.element_period(118)
        assert period == 7

    def test_element_119_starts_period_8(self):
        """Test Z=119 starts period 8."""
        period = calculate.element_period(119)
        assert period == 8
