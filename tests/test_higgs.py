"""
Tests for Higgs Potential and VEV calculations.
"""

import subprocess
import pytest
import numpy as np
import os
import sys

def test_higgs_potential_explanation_runs():
    """Verify that particle_higgs_potential.py runs without errors and articulates the Higgs potential."""
    try:
        result = subprocess.run(
            [sys.executable, 'physics/particle_higgs_potential.py'],
            capture_output=True, text=True, check=True
        )
        assert "Conclusion: Incorporating these quantum and symmetry considerations, the effective potential of the Higgs field φ takes the form:" in result.stdout
        assert result.returncode == 0
    except subprocess.CalledProcessError as e:
        pytest.fail(f"particle_higgs_potential.py failed with error: {e.stderr}")

def test_coleman_weinberg_stability():
    """Verify that particle_coleman_weinberg.py confirms vacuum stability with positive total quartic coupling."""
    try:
        result = subprocess.run(
            [sys.executable, 'physics/particle_coleman_weinberg.py'],
            capture_output=True, text=True, check=True
        )
        assert "SUCCESS: Quantum corrections stabilize the vacuum!" in result.stdout
        assert "Total Quartic Coupling: 20.3171" in result.stdout # Check for specific value
        assert result.returncode == 0
    except subprocess.CalledProcessError as e:
        pytest.fail(f"particle_coleman_weinberg.py failed with error: {e.stderr}")

from physics.particle_higgs_vev import calculate_geometric_vev

def test_calculate_vev_accuracy():
    """Verify that calculate_vev.py predicts the Higgs VEV with acceptable accuracy."""
    try:
        predicted_vev = calculate_geometric_vev() # Directly call the function
        
        expected_vev = 245.1739 # GeV, based on improved geometric constants
        tolerance = 0.1 # 0.1 GeV
        
        assert np.isclose(predicted_vev, expected_vev, atol=tolerance), f"Predicted VEV {predicted_vev} not close to expected {expected_vev}"
    except Exception as e:
        pytest.fail(f"calculate_geometric_vev() failed with error: {e}")
