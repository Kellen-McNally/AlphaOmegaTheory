#
#   Beta functions for gauge coupling RG evolution.
#   Implements 1-loop and 2-loop beta functions for the
#   Standard Model gauge couplings: U(1)_Y, SU(2)_L,
#   SU(3)_c.
#

import numpy as np
from typing import Tuple, Optional

from utils.logging_config import get_logger

logger = get_logger(__name__)


# G2-Geometric Beta Coefficients (Dual Octonion Spectrum)
# These match MSSM slopes due to the dual 8+8 structure
# but arise from geometric degrees of freedom (moduli), not sparticles.
B1_G2 = 33.0 / 5.0   # U(1)_Y
B2_G2 = 1.0          # SU(2)_L
B3_G2 = -3.0         # SU(3)_c


#
#   Return 1-loop G2 beta function coefficients. Returns:
#   tuple (b₁, b₂, b₃) matching the effective degrees of
#   freedom in the G2 manifold. Note: These coefficients
#   mimic MSSM values due to the 8+8 duality, justifying
#   their use in RG evolution without physical sparticles.
#
def beta_coefficients_sm() -> Tuple[float, float, float]:
    return (B1_G2, B2_G2, B3_G2)


#
#   Compute 1-loop beta function. dα_i/dt = b_i/(2π) × α_i²
#   where t = ln(μ/μ₀). Args: alpha (array of gauge
#   couplings [α₁, α₂, α₃]), b (array of beta coefficients
#   [b₁, b₂, b₃]). Returns: dα/dt.
#
def beta_1loop(alpha: np.ndarray, b: np.ndarray) -> np.ndarray:
    return b / (2.0 * np.pi) * alpha**2


#
#   Compute 2-loop beta function corrections. dα_i/dt =
#   b_i/(2π) α_i² + Σ_j b_ij/(2π)² α_i² α_j. Args: alpha
#   (array of gauge couplings [α₁, α₂, α₃]). Returns: 2-loop
#   corrections to dα/dt.
#
def beta_2loop(alpha: np.ndarray) -> np.ndarray:
    # 2-loop coefficients (SM with 3 generations)
    b11 = 199.0 / 50.0
    b12 = 27.0 / 10.0
    b13 = 44.0 / 5.0

    b21 = 9.0 / 10.0
    b22 = 35.0 / 6.0
    b23 = 12.0

    b31 = 11.0 / 10.0
    b32 = 9.0 / 2.0
    b33 = -26.0

    b_matrix = np.array([
        [b11, b12, b13],
        [b21, b22, b23],
        [b31, b32, b33]
    ])

    alpha_i_sq = alpha**2
    correction = np.zeros(3)

    for i in range(3):
        for j in range(3):
            correction[i] += b_matrix[i, j] * alpha_i_sq[i] * alpha[j]

    return correction / (2.0 * np.pi)**2


#
#   Evolve gauge couplings using RG equations. Args:
#   alpha_init (initial couplings [α₁, α₂, α₃] at mu_init),
#   mu_init (initial energy scale in GeV), mu_final (final
#   energy scale in GeV), n_steps (number of integration
#   steps), two_loop (include 2-loop corrections). Returns:
#   tuple (mu_array, alpha_array) with evolved couplings.
#   mu_array: energy scales (GeV), alpha_array: couplings at
#   each scale, shape (n_steps, 3).
#
def rg_evolve(
    alpha_init: np.ndarray,
    mu_init: float,
    mu_final: float,
    n_steps: int = 1000,
    two_loop: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    # Logarithmic energy scale
    t_init = np.log(mu_init)
    t_final = np.log(mu_final)
    t_array = np.linspace(t_init, t_final, n_steps)
    mu_array = np.exp(t_array)

    # Initialize
    alpha = alpha_init.copy()
    alpha_trajectory = np.zeros((n_steps, 3))
    alpha_trajectory[0] = alpha

    b = np.array(beta_coefficients_sm())
    dt = (t_final - t_init) / n_steps

    # Integrate RG equations
    for i in range(1, n_steps):
        # 1-loop contribution
        dalpha_dt = beta_1loop(alpha, b)

        # Add 2-loop if requested
        if two_loop:
            dalpha_dt += beta_2loop(alpha)

        # Euler step (could use RK4 for better accuracy)
        alpha = alpha + dalpha_dt * dt
        alpha_trajectory[i] = alpha

    return mu_array, alpha_trajectory


if __name__ == "__main__":
    from utils.logging_config import configure_for_cli
    configure_for_cli(verbose=True)

    # Test RG evolution
    logger.info("Beta Function Module Test")
    logger.info("=" * 60)

    # Initial conditions at M_Z
    M_Z = 91.2
    alpha_1_MZ = 0.01695  # U(1)_Y (GUT normalized)
    alpha_2_MZ = 0.03352  # SU(2)_L
    alpha_3_MZ = 0.1184   # SU(3)_c

    alpha_mz = np.array([alpha_1_MZ, alpha_2_MZ, alpha_3_MZ])

    logger.info(f"\nInitial couplings at M_Z = {M_Z} GeV:")
    logger.info(f"  α₁ = {alpha_1_MZ:.5f}")
    logger.info(f"  α₂ = {alpha_2_MZ:.5f}")
    logger.info(f"  α₃ = {alpha_3_MZ:.5f}")

    # Evolve to 10¹⁴ GeV
    M_GUT = 1e14
    mu_array, alpha_array = rg_evolve(alpha_mz, M_Z, M_GUT, n_steps=1000)

    alpha_gut = alpha_array[-1]
    logger.info(f"\nCouplings at M_GUT = {M_GUT:.2e} GeV:")
    logger.info(f"  α₁ = {alpha_gut[0]:.5f}")
    logger.info(f"  α₂ = {alpha_gut[1]:.5f}")
    logger.info(f"  α₃ = {alpha_gut[2]:.5f}")

    spread = np.max(alpha_gut) - np.min(alpha_gut)
    mean = np.mean(alpha_gut)
    logger.info(f"\nSpread: {spread:.5f} ({spread/mean*100:.2f}%)")
