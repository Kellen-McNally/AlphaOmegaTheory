"""
Derivation of the Grand Unified Coupling Constant (α_GUT)

Calculates the geometric value of the unification coupling constant 
based on G2 group invariants.
"""

import numpy as np
from core.constants import TRIALITY, DIM_G2, ALPHA_GUT
from utils.logging_config import get_logger

logger = get_logger(__name__)

# Experimental value from RG unification
ALPHA_GUT_EXP = 0.024
ALPHA_GUT_ERR = 0.002

def alpha_gut_rigorous() -> dict:
    """
    Calculate α_GUT = 1/(τ × dim(G₂)).
    """
    # Pure geometric formula
    denominator = TRIALITY * DIM_G2
    alpha_gut_predicted = 1.0 / denominator

    # Experimental comparison
    error_abs = abs(alpha_gut_predicted - ALPHA_GUT_EXP)
    error_rel = error_abs / ALPHA_GUT_EXP
    error_pct = error_rel * 100
    sigma = error_abs / ALPHA_GUT_ERR

    return {
        'formula': {
            'tau': TRIALITY,
            'dim': DIM_G2,
            'denominator': denominator,
            'value': alpha_gut_predicted,
        },
        'experiment': {
            'value': ALPHA_GUT_EXP,
            'uncertainty': ALPHA_GUT_ERR,
        },
        'comparison': {
            'error_percent': error_pct,
            'sigma': sigma,
        }
    }

if __name__ == '__main__':
    from utils.logging_config import configure_for_cli
    configure_for_cli(verbose=True)

    result = alpha_gut_rigorous()
    form = result['formula']
    comp = result['comparison']

    logger.info("Grand Unified Coupling (α_GUT)")
    logger.info("==============================")
    logger.info(f"Formula: 1 / (τ × dim) = 1 / ({form['tau']} × {form['dim']})")
    logger.info(f"Predicted Value: {form['value']:.6f} (1/{form['denominator']})")
    logger.info(f"Experimental Value: {result['experiment']['value']} ± {result['experiment']['uncertainty']}")
    logger.info(f"Error: {comp['error_percent']:.2f}% ({comp['sigma']:.1f}σ)")
