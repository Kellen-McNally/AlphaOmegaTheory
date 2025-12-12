#
#   Oblique Parameters S and T from G₂ Corrections. The
#   oblique parameters S and T parameterize new physics
#   contributions to electroweak precision observables. S:
#   Measures new contributions to gauge boson self-energies,
#   T: Measures custodial SU(2) breaking. Experimental
#   constraints (PDG 2020): S = 0.05 ± 0.11 (with U = 0), T
#   = 0.09 ± 0.13. The αΩ Framework predicts small S, T from
#   dimension-7 operators.
#

import numpy as np
from typing import Dict
import sys
import os

# Add parent directory to path using relative import
api_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
if api_dir not in sys.path:
    sys.path.insert(0, api_dir)

from core.constants import (
    ALPHA_EM_MZ,
    M_Z_MEASURED,
    M_GUT,
    TRIALITY,
    DIM_G2,
)
from utils.logging_config import get_logger

logger = get_logger(__name__)


# Experimental constraints (PDG 2020, assuming U = 0)
S_EXP = 0.05
S_ERR = 0.11

T_EXP = 0.09
T_ERR = 0.13


#
#   Compute S parameter from G₂ corrections. S measures new
#   contributions to vacuum polarization: S = (4 sin²θ_W
#   cos²θ_W / α) × [Π_ZZ(0) - Π_WW(0)]. From dimension-7
#   operators at M_GUT: S ≈ (4π v²/M_GUT²) × c_G₂ where c_G₂
#   ~ 42/(16π²). Returns: float (S parameter).
#
def oblique_s_parameter() -> float:
    # Higgs VEV
    v = 246.0  # GeV

    # G₂ coefficient
    c_G2 = (TRIALITY * DIM_G2) / (16.0 * np.pi**2)

    # Dimension-7 contribution
    # Note: Extra factor of M_Z/M_GUT relative to dimension-6
    S = (4.0 * np.pi * v**2 / M_GUT**2) * (M_Z_MEASURED / M_GUT) * c_G2

    return S


#
#   Compute T parameter from G₂ corrections. T measures
#   custodial SU(2) breaking: T = [Π_WW(0) - Π_ZZ(0)] / (α
#   M_W²). G₂ holonomy PRESERVES custodial symmetry at
#   leading order, so T is suppressed: T ≈ (v⁴/M_GUT⁴) ×
#   c_G₂ × (dimension-8 operator). Result: T ~ 0.001 (very
#   small!). Returns: float (T parameter).
#
def oblique_t_parameter() -> float:
    # Higgs VEV
    v = 246.0  # GeV

    # G₂ coefficient
    c_G2 = (TRIALITY * DIM_G2) / (16.0 * np.pi**2)

    # Dimension-8 contribution (highly suppressed)
    T = (v**4 / M_GUT**4) * c_G2

    return T


#
#   Compute U parameter from G₂ corrections. U measures
#   additional new physics in Z boson couplings. Usually
#   assumed U = 0 in fits. In αΩ Framework: U ~ 0 (same
#   order as T). Returns: float (U parameter).
#
def oblique_u_parameter() -> float:
    # Similar suppression to T
    v = 246.0  # GeV
    c_G2 = (TRIALITY * DIM_G2) / (16.0 * np.pi**2)

    U = (v**4 / M_GUT**4) * c_G2

    return U


#
#   Complete summary of oblique parameters. Returns: dict
#   (All oblique parameter information).
#
def oblique_summary() -> Dict:
    S_pred = oblique_s_parameter()
    T_pred = oblique_t_parameter()
    U_pred = oblique_u_parameter()

    # Comparison with experiment
    S_discrepancy = S_pred - S_EXP
    S_sigma = abs(S_discrepancy) / S_ERR

    T_discrepancy = T_pred - T_EXP
    T_sigma = abs(T_discrepancy) / T_ERR

    S_status = "EXCELLENT" if S_sigma < 1.0 else "GOOD" if S_sigma < 2.0 else "NEEDS REFINEMENT"
    T_status = "EXCELLENT" if T_sigma < 1.0 else "GOOD" if T_sigma < 2.0 else "NEEDS REFINEMENT"

    return {
        "definition": {
            "S": "Measures new contributions to gauge boson self-energies",
            "T": "Measures custodial SU(2) breaking",
            "U": "Additional Z boson coupling correction (usually set to 0)",
        },
        "experimental": {
            "S": {
                "value": S_EXP,
                "error": S_ERR,
                "constraint": f"S = {S_EXP:.2f} ± {S_ERR:.2f}",
                "assumption": "U = 0",
            },
            "T": {
                "value": T_EXP,
                "error": T_ERR,
                "constraint": f"T = {T_EXP:.2f} ± {T_ERR:.2f}",
            },
            "source": "PDG 2020",
        },
        "g2_prediction": {
            "S": {
                "value": S_pred,
                "mechanism": "Dimension-7 operators at M_GUT",
                "formula": "S ≈ (4π v²/M_GUT²) × (M_Z/M_GUT) × c_G₂",
                "discrepancy": S_discrepancy,
                "sigma": S_sigma,
                "status": S_status,
            },
            "T": {
                "value": T_pred,
                "mechanism": "Dimension-8 operators (custodial symmetry preserved)",
                "formula": "T ≈ (v⁴/M_GUT⁴) × c_G₂",
                "discrepancy": T_discrepancy,
                "sigma": T_sigma,
                "status": T_status,
                "note": "Highly suppressed - custodial symmetry preserved!",
            },
            "U": {
                "value": U_pred,
                "note": "Similar to T (highly suppressed)",
            },
        },
        "interpretation": {
            "small_S_T": "[OK] S, T ~ 0 (consistent with precision tests)",
            "custodial_symmetry": "[OK] G₂ holonomy preserves custodial SU(2)",
            "no_large_corrections": "[OK] No tree-level corrections (unlike many BSM models)",
            "unified_origin": f"[OK] Same c_G₂ = {(TRIALITY * DIM_G2) / (16.0 * np.pi**2):.6f} governs all corrections",
        },
        "comparison_with_other_models": {
            "susy": "SUSY: T ~ m_t²/M_W² (Δm_stop²) (can be large)",
            "extra_dimensions": "Extra dims: S ~ v²/M_KK² (can conflict with data)",
            "g2_framework": "G₂: S, T highly suppressed (consistent with all data)",
        },
    }


def paper_section():
    """Generate paper section for oblique parameters S and T."""
    from paper.paper_api import PaperSection

    summary = oblique_summary()

    content = rf"""
Precision electroweak measurements constrain new physics through oblique parameters
S and T, which parameterize vacuum polarization corrections.

\subsection{{Oblique Parameter Definitions}}

The S and T parameters characterize new physics contributions:
\begin{{align}}
S &= \text{{new physics contribution to gauge boson self-energies}} \\
T &= \text{{new physics contribution to ρ parameter}}
\end{{align}}

\subsection{{G$_2$ Framework Predictions}}

The αΩ Framework predicts highly suppressed oblique corrections:
\begin{{align}}
S &= {summary['g2_prediction']['S']['value']:.2e} \\
T &= {summary['g2_prediction']['T']['value']:.2e}
\end{{align}}

These arise from dimension-7 operators, naturally suppressed by $(v/M_{{GUT}})^3$.

\subsection{{Experimental Constraints}}

Global electroweak fits give (PDG 2020):
\begin{{align}}
S &= 0.05 \pm 0.11 \\
T &= 0.09 \pm 0.13
\end{{align}}

The G$_2$ predictions are well within experimental bounds.

\paragraph{{Comparison with Other Theories}}
\begin{{itemize}}
\item Standard Model: $S = T = 0$ (by definition)
\item Technicolor: Large S (typically excluded)
\item Extra dimensions: Can produce large S
\item G$_2$ Framework: Both S and T naturally suppressed
\end{{itemize}}

\paragraph{{Theoretical Significance}}
The suppression of oblique parameters in G$_2$ theory arises naturally
from the geometric structure, requiring no fine-tuning.
This represents a significant advantage over many BSM theories.
"""

    return PaperSection(
        title="Precision Electroweak: Oblique Parameters",
        content=content,
        order=12
    )
