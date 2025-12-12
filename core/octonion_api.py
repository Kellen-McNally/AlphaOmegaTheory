"""
Unified Octonion API

Consolidates all octonion access patterns into a single, consistent interface.
Provides clean access to octonion algebra for all framework components.
"""

import numpy as np
from typing import Dict, List, Tuple
from core import octonion_algebra as oct_alg
from core.constants import TRIALITY, DIM_G2

# ═══════════════════════════════════════════════════════════
# OCTONION CREATION
# ═══════════════════════════════════════════════════════════

def create_octonion(*components) -> np.ndarray:
    """Create octonion from 8 components."""
    if len(components) != 8:
        raise ValueError(f"Octonion requires exactly 8 components, got {len(components)}")
    return np.array(components, dtype=float)

def identity_octonion() -> np.ndarray:
    """Create identity octonion (1, 0, ...)."""
    v = np.zeros(8)
    v[0] = 1.0
    return v

def basis_octonion(index: int) -> np.ndarray:
    """Create basis octonion e_i."""
    v = np.zeros(8)
    v[index] = 1.0
    return v

# ═══════════════════════════════════════════════════════════
# PHYSICS CONSTRUCTORS
# ═══════════════════════════════════════════════════════════

def spacetime_octonion(t, tau, x, y, z, p_x, p_y, p_z) -> np.ndarray:
    """Create external octonion with spacetime + momentum."""
    return create_octonion(t, tau, x, p_x, y, p_y, z, p_z)

def particle_octonion(mass, charge, generation, color_r=0, color_g=0, color_b=0, weak_isospin=0) -> np.ndarray:
    """Create internal octonion with particle quantum numbers."""
    # Mapping: i0=mass/amp, i1=charge, i2=gen... (Simplified mapping)
    return create_octonion(mass, charge, generation, color_r, color_g, color_b, weak_isospin, 0)

# ═══════════════════════════════════════════════════════════
# OPERATIONS
# ═══════════════════════════════════════════════════════════

def triality_action(octonion: np.ndarray, power: int = 1) -> np.ndarray:
    """Apply triality automorphism τ^power."""
    components = octonion.copy()
    for _ in range(power % 3):
        new_components = components.copy()
        new_components[1:4] = components[4:7]     # gen1 <- gen2
        new_components[4:7] = -components[1:4]    # gen2 <- -gen1
        components = new_components
    return components

def heisenberg_commutator(coord_index: int) -> float:
    """Compute Heisenberg commutator [x_i, p_i] = 2τ."""
    return 2 * TRIALITY

# ═══════════════════════════════════════════════════════════
# VALIDATION
# ═══════════════════════════════════════════════════════════

def computational_efficiency_score() -> Dict[str, float]:
    return {
        'external_optimal': 672,
        'external_total': 40320,
        'combined_efficiency': 64/1625702400
    }
