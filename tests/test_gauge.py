"""
Tests for Gauge theory related derivations.
"""

import subprocess
import pytest
import os
import sys
import numpy as np

# Test for physics/particles/gauge/field_mapping.py
def test_field_mapping_runs():
    """Verify that field_mapping.py runs without errors and produces a field map."""
    try:
        result = subprocess.run(
            [sys.executable, 'physics/particle_field_mapping.py'],
            capture_output=True, text=True, check=True
        )
        assert "Proposed Field Map:" in result.stdout
        assert "Gluons (g1-g8)" in result.stdout
        assert "Weak Bosons (W+, W-, Z)" in result.stdout
        assert "Photon (gamma)" in result.stdout
        assert "Higgs (h)" in result.stdout
        assert result.returncode == 0
    except subprocess.CalledProcessError as e:
        pytest.fail(f"particle_field_mapping.py failed with error: {e.stderr}")

# Test for physics/particle_scattering.py
def test_scattering_amplitude_runs():
    """Verify that particle_scattering.py runs without errors and articulates the vertex structure."""
    try:
        result = subprocess.run(
            [sys.executable, 'physics/particle_scattering.py'],
            capture_output=True, text=True, check=True
        )
        assert "Deriving Interaction Vertex from Cubic Term..." in result.stdout
        assert "Status: Verified that External units act as operators on Internal units." in result.stdout
        assert result.returncode == 0
    except subprocess.CalledProcessError as e:
        pytest.fail(f"particle_scattering.py failed with error: {e.stderr}")

# Test for physics/particle_lagrangian.py
def test_lagrangian_reduction_runs():
    """Verify that particle_lagrangian.py runs without errors and confirms the reduction."""
    try:
        result = subprocess.run(
            [sys.executable, 'physics/particle_lagrangian.py'],
            capture_output=True, text=True, check=True
        )
        assert "SUCCESS: Cubic geometric term -> Yang-Mills cubic interaction." in result.stdout
        assert "SUCCESS: Geometric kinetic term -> Dirac Lagrangian." in result.stdout
        assert result.returncode == 0
    except subprocess.CalledProcessError as e:
        pytest.fail(f"particle_lagrangian.py failed with error: {e.stderr}")

# Test for physics/particle_symbolic_sedenion.py
def test_symbolic_sedenion_runs():
    """Verify that particle_symbolic_sedenion.py runs without errors and confirms non-associativity."""
    try:
        result = subprocess.run(
            [sys.executable, 'physics/particle_symbolic_sedenion.py'],
            capture_output=True, text=True, check=True
        )
        assert "SUCCESS: Associator is non-zero (Non-associative algebra confirmed)." in result.stdout
        assert "Result: 0 for imaginary A. (Alternative property holds)" in result.stdout
        assert result.returncode == 0
    except subprocess.CalledProcessError as e:
        pytest.fail(f"symbolic_sedenion.py failed with error: {e.stderr}")
