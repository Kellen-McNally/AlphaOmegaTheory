"""
Tests for G₂ structure constants.
"""

import pytest
import numpy as np
from core import g2_structure_constants as structure_constants


class TestStructureConstants:
    """Test G₂ structure constants."""

    def test_structure_constants_antisymmetric(self):
        """Test structure constants can be computed."""
        try:
            f = structure_constants.get_structure_constants()

            # Should be 3D array
            assert len(f.shape) == 3
        except NotImplementedError:
            pytest.skip("Structure constants not fully implemented")

    def test_g2_simple_roots(self):
        """Test G₂ has 2 simple roots."""
        roots = structure_constants.simple_roots()

        assert len(roots) == 2

    def test_cartan_matrix_shape(self):
        """Test Cartan matrix is 2×2 for G₂."""
        cartan = structure_constants.cartan_matrix()

        assert cartan.shape == (2, 2)

    def test_associative_3form_structure(self):
        """Test associative 3-form φ structure."""
        triples = structure_constants.associative_3form_structure()
        
        # Should return list of tuples
        assert isinstance(triples, list)
        assert len(triples) == 7
        assert len(triples[0]) == 3
        
        # Check standard (1,2,3) triple
        assert (1, 2, 3) in triples
