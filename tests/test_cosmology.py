"""
Tests for Cosmology related derivations.
"""

import subprocess
import pytest
import os
import numpy as np
# from unittest.mock import MagicMock, patch # No longer needed for this test

# Test for physics/cosmology/dark_matter_dynamics.py
def test_dark_matter_dynamics_runs():
    """Verify that dark_matter_dynamics.py runs without errors and handles matplotlib gracefully."""
    try:
        # Clean up plot file if it was created in a previous run.
        if os.path.exists('galaxy_rotation_curve.png'):
            os.remove('galaxy_rotation_curve.png') 
            
        result = subprocess.run(
            ['python3', 'physics/cosmo_dark_matter.py'],
            capture_output=True, text=True, check=True,
            env={**os.environ, 'PYTHONPATH': os.environ.get('PYTHONPATH', '') + ':' + '/home/kfm0/.local/lib/python3.12/site-packages'}
        )
        assert "SUCCESS: Galaxy rotation curve simulated." in result.stdout
        
        # Check if plot was skipped or generated based on stderr output
        if "Skipping plot generation due to missing matplotlib." in result.stderr:
            assert not os.path.exists('galaxy_rotation_curve.png'), "Plot file should not exist if matplotlib was skipped."
        else:
            assert "Output saved to 'galaxy_rotation_curve.png'." in result.stdout, "Plot success message not found."
            assert os.path.exists('galaxy_rotation_curve.png'), "Plot file should exist if matplotlib was not skipped."
        
        assert result.returncode == 0
    except subprocess.CalledProcessError as e:
        pytest.fail(f"cosmo_dark_matter.py failed with error: {e.stderr}")

# Test for physics/cosmo_gravity.py
def test_einstein_hilbert_runs():
    """Verify that cosmo_gravity.py runs without errors and articulates the GR derivation."""
    try:
        result = subprocess.run(
            ['python3', 'physics/cosmo_gravity.py'],
            capture_output=True, text=True, check=True,
            env={**os.environ, 'PYTHONPATH': os.environ.get('PYTHONPATH', '') + ':' + '/home/kfm0/.local/lib/python3.12/site-packages'}
        )
        assert "SUCCESS: The Einstein-Hilbert action emerges from the Sedenion curvature term." in result.stdout
        assert result.returncode == 0
    except subprocess.CalledProcessError as e:
        pytest.fail(f"cosmo_gravity.py failed with error: {e.stderr}")

# Test for physics/cosmo_equilibrium.py
def test_equilibrium_state_runs_and_values():
    """Verify that cosmo_equilibrium.py runs and calculates correct H_eq and T_eq."""
    try:
        result = subprocess.run(
            ['python3', 'physics/cosmo_equilibrium.py'],
            capture_output=True, text=True, check=True,
            env={**os.environ, 'PYTHONPATH': os.environ.get('PYTHONPATH', '') + ':' + '/home/kfm0/.local/lib/python3.12/site-packages'}
        )
        assert "SUCCESS: Equilibrium state properties (H_eq, T_eq) calculated based on geometric constants." in result.stdout
        
        # Extract H_eq and T_eq and assert their values
        output_lines = result.stdout.splitlines()
        h_eq_line = [line for line in output_lines if "Predicted Equilibrium Hubble Constant (H_eq):" in line]
        t_eq_line = [line for line in output_lines if "Predicted Equilibrium Temperature (T_eq):" in line]
        
        assert h_eq_line, "H_eq line not found in output."
        assert t_eq_line, "T_eq line not found in output."
        
        h_eq_str = h_eq_line[0].split(':')[1].strip().split(' ')[0]
        t_eq_str = t_eq_line[0].split(':')[1].strip().split(' ')[0]
        
        predicted_h_eq = float(h_eq_str)
        predicted_t_eq = float(t_eq_str)
        
        expected_h_eq = 55.885
        expected_t_eq = 2.20e-30
        tolerance_h = 0.001
        tolerance_t = 1e-32 # For scientific notation comparison
        
        assert np.isclose(predicted_h_eq, expected_h_eq, atol=tolerance_h), f"Predicted H_eq {predicted_h_eq} not close to expected {expected_h_eq}"
        assert np.isclose(predicted_t_eq, expected_t_eq, atol=tolerance_t), f"Predicted T_eq {predicted_t_eq} not close to expected {expected_t_eq}"
        assert result.returncode == 0
    except subprocess.CalledProcessError as e:
        pytest.fail(f"equilibrium_state.py failed with error: {e.stderr}")