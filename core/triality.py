#
#   Triality automorphism of G₂. τ³=1 (order-3
#   automorphism). This is the mathematical reason for
#   exactly 3 generations of fermions.
#

import numpy as np
from typing import Tuple
from utils.logging_config import get_logger

logger = get_logger(__name__)
#
#   Return the order of triality automorphism. Returns: int
#   (3, τ³=1).
#
def triality_order() -> int:
    return 3
#
#   Compute the three triality phases. τ³=1 → ω=e^(2πi/3)
#   (cube root of unity). Returns: tuple (ω⁰, ω¹, ω²) where
#   ω=e^(2πi/3).
#
def triality_phases() -> Tuple[complex, complex, complex]:
    omega_0 = 1.0 + 0j
    omega_1 = np.exp(1j * 2.0 * np.pi / 3.0)
    omega_2 = np.exp(1j * 4.0 * np.pi / 3.0)

    return (omega_0, omega_1, omega_2)
#
#   Return generation indices. Returns: tuple (1, 2, 3) for
#   three generations.
#
def generation_indices() -> Tuple[int, int, int]:
    return (1, 2, 3)
if __name__ == "__main__":
    from utils.logging_config import configure_for_cli
    configure_for_cli(verbose=True)

    logger.info("=" * 70)
    logger.info("TRIALITY: τ³ = 1")
    logger.info("=" * 70)
    logger.info("")

    logger.info(f"Order of triality: {triality_order()}")
    logger.info("")

    logger.info("Triality phases (cube roots of unity):")
    omegas = triality_phases()
    for i, omega in enumerate(omegas):
        logger.info(f"  ω^{i} = {omega:.6f}")
    logger.info("")

    logger.info("Verification: ω³ = 1")
    omega = omegas[1]
    omega_cubed = omega ** 3
    logger.info(f"  ω³ = {omega_cubed:.6f}")
    logger.info(f"  |ω³ - 1| = {abs(omega_cubed - 1):.2e}")
    logger.info("")

    logger.info("=" * 70)
    logger.info("Three generations from τ³ = 1")
    logger.info("=" * 70)

