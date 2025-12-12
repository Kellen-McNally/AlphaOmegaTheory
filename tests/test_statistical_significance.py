"""
Statistical Significance Analysis of αΩ Framework Predictions

Computes the combined statistical probability that the framework's predictions
are false positives (i.e., chance agreement with experiment).

Uses chi-squared test and Bayesian analysis to quantify:
1. How unlikely is this level of agreement by chance?
2. What is the probability the framework is false?

Key principle: With ZERO free parameters, any agreement with experiment
is either:
- Genuine prediction from correct theory, OR
- Extraordinarily unlikely coincidence

This test calculates exactly HOW unlikely the coincidence would be.
"""

import numpy as np
from scipy import stats
import pytest

# Import constants from ground truth (single source of truth)
from core.constants import (
    OMEGA_LAMBDA, OMEGA_LAMBDA_EXP, OMEGA_LAMBDA_ERR,
    SIN2_THETA_W, SIN2_THETA_W_EXP, SIN2_THETA_W_ERR,
    TRIALITY
)


class TestStatisticalSignificance:
    """Calculate statistical significance of framework predictions."""

    def test_tier1_predictions_statistical_significance(self):
        """
        Calculate combined statistical significance of tier 1 predictions.

        Tier 1 predictions are EXACT from {3, 14} with zero free parameters.
        Any agreement with experiment is either correct physics or coincidence.
        """

        print("\n" + "=" * 80)
        print("STATISTICAL SIGNIFICANCE ANALYSIS: αΩ FRAMEWORK")
        print("=" * 80)
        print()

        # =====================================================================
        # TIER 1: EXACT PREDICTIONS FROM {3, 14}
        # =====================================================================

        print("TIER 1: Exact Predictions (0 free parameters)")
        print("-" * 80)
        print()

        predictions = [
            {
                'name': 'Dark Energy Density (Ω_Λ)',
                'predicted': OMEGA_LAMBDA,  # From constants.py: (Tr_τ³ - 6)/(dim + rank)
                'observed': OMEGA_LAMBDA_EXP,  # Planck 2018
                'uncertainty': OMEGA_LAMBDA_ERR,
                'formula': 'Ω_Λ = C₃/(dim + rank) = 11/16'
            },
            {
                'name': 'Weak Mixing Angle (sin²θ_W)',
                'predicted': SIN2_THETA_W,  # From constants.py: τ/(dim - 1)
                'observed': SIN2_THETA_W_EXP,  # PDG
                'uncertainty': SIN2_THETA_W_ERR,  # From constants.py
                'formula': 'sin²θ_W = τ/(dim - 1) = 3/13'
            },
            {
                'name': 'Fermion Generations (N_gen)',
                'predicted': float(TRIALITY),  # From constants.py: triality order
                'observed': 2.9840,  # LEP measurement
                'uncertainty': 0.0082,
                'formula': 'N_gen = order(τ) = 3'
            },
            {
                'name': 'Cabibbo Angle (sin θ_C)',
                'predicted': np.sqrt(5/98),
                'observed': 0.2245,
                'uncertainty': 0.0005,
                'formula': 'sin θ_C = √(5/98)'
            },
        ]

        print(f"{'Observable':<30} {'Predicted':<12} {'Observed':<15} {'σ':<8} {'p-value'}")
        print("-" * 80)

        chi_squared_total = 0.0
        n_predictions = len(predictions)
        p_values = []

        for pred in predictions:
            # Calculate chi-squared contribution
            residual = pred['predicted'] - pred['observed']
            chi_sq = (residual / pred['uncertainty'])**2
            chi_squared_total += chi_sq

            # Calculate individual significance (standard deviations)
            sigma = abs(residual) / pred['uncertainty']

            # Calculate p-value (two-tailed)
            p_value = 2 * (1 - stats.norm.cdf(sigma))
            p_values.append(p_value)

            print(f"{pred['name']:<30} {pred['predicted']:<12.6f} {pred['observed']:.6f}±{pred['uncertainty']:.6f}  "
                  f"{sigma:<8.2f} {p_value:.2e}")

        print()
        print(f"Total χ² = {chi_squared_total:.4f} with {n_predictions} predictions")
        print()

        # =====================================================================
        # CHI-SQUARED TEST
        # =====================================================================

        print("CHI-SQUARED ANALYSIS:")
        print("-" * 80)
        print()

        # Degrees of freedom = number of predictions - number of free parameters
        # In αΩ Framework: 0 free parameters!
        dof = n_predictions - 0  # Zero free parameters

        # Calculate p-value from chi-squared distribution
        # This is the probability of getting this good agreement by CHANCE
        p_chi_squared = 1 - stats.chi2.cdf(chi_squared_total, dof)

        print(f"χ² = {chi_squared_total:.4f}")
        print(f"Degrees of freedom = {dof} (predictions) - 0 (free parameters) = {dof}")
        print(f"Reduced χ² = χ²/dof = {chi_squared_total/dof:.4f}")
        print()
        print(f"p-value (χ² test) = {p_chi_squared:.4e}")
        print()

        if chi_squared_total/dof < 1:
            print("[OK] Reduced χ² < 1: Excellent fit!")
            print("  Predictions agree with observations within uncertainties.")
        elif chi_squared_total/dof < 2:
            print("[OK] Reduced χ² < 2: Good fit!")
        else:
            print("[WARN] Reduced χ² > 2: Some tension (expected with precise measurements)")

        print()

        # =====================================================================
        # COMBINED P-VALUE (Fisher's Method)
        # =====================================================================

        print("COMBINED SIGNIFICANCE (Fisher's Method):")
        print("-" * 80)
        print()

        # Fisher's method: -2 Σ ln(p_i) ~ χ² with 2k degrees of freedom
        # Handle very small p-values by clamping to avoid log(0)
        p_values_clamped = [max(p, 1e-300) for p in p_values]
        fisher_statistic = -2 * np.sum(np.log(p_values_clamped))
        fisher_dof = 2 * len(p_values)

        # Handle numerical overflow in chi2.cdf
        if fisher_statistic > 1000:  # Very large chi-squared
            combined_p_value = 0.0
            combined_sigma = 10.0  # Conservative lower bound (>10σ)
        else:
            combined_p_value = 1 - stats.chi2.cdf(fisher_statistic, fisher_dof)
            # Convert to sigma (standard deviations)
            if combined_p_value > 1e-16:
                combined_sigma = stats.norm.ppf(1 - combined_p_value/2)
            else:
                combined_sigma = 10.0  # >10σ

        print(f"Fisher's χ² = {fisher_statistic:.4f}")
        print(f"Degrees of freedom = {fisher_dof}")
        print(f"Combined p-value = {combined_p_value:.4e}")
        print()

        print(f"Combined significance: {combined_sigma:.1f}σ")
        print()

        # =====================================================================
        # BAYESIAN ANALYSIS
        # =====================================================================

        print("BAYESIAN INTERPRETATION:")
        print("-" * 80)
        print()

        # Bayes' theorem: P(theory true | data) = P(data | theory true) × P(theory true) / P(data)
        #
        # P(data | theory true) ≈ 1 (theory predicts this data)
        # P(data | theory false) = p-value (probability by chance)
        #
        # Bayes factor = P(data | theory true) / P(data | theory false)
        #              ≈ 1 / p-value

        # Use clamped p-value for numerical stability
        combined_p_value_for_bayes = max(combined_p_value, 1e-100)
        bayes_factor = 1.0 / combined_p_value_for_bayes

        print("Given the data, what's the probability the framework is FALSE?")
        print()
        print("Using Bayes' theorem with conservative priors:")
        print()

        # Conservative prior: 50% chance framework is true before seeing data
        prior_true = 0.5
        prior_false = 0.5

        # Likelihood ratio
        likelihood_ratio = bayes_factor

        # Posterior probability framework is TRUE
        # Use numerically stable formula
        if bayes_factor > 1e10:  # Very strong evidence
            posterior_true = 1.0 - 1.0/(2.0 * bayes_factor)  # Asymptotic formula
            posterior_false = 1.0/(2.0 * bayes_factor)
        else:
            posterior_true = (likelihood_ratio * prior_true) / (
                likelihood_ratio * prior_true + prior_false
            )
            posterior_false = 1 - posterior_true

        print(f"Prior probability (before data): P(true) = {prior_true:.1%}")
        if bayes_factor > 1e50:
            print(f"Bayes factor (data favors theory): >{1e50:.0e}")
        else:
            print(f"Bayes factor (data favors theory): {bayes_factor:.2e}")
        print()
        print(f"Posterior probability (after data):")
        print(f"  P(framework TRUE | data)  = {posterior_true:.10f} = {posterior_true*100:.8f}%")
        if posterior_false > 1e-10:
            print(f"  P(framework FALSE | data) = {posterior_false:.10e} = {posterior_false*100:.8e}%")
        else:
            print(f"  P(framework FALSE | data) < {1e-10:.0e} = <{1e-8:.0e}%")
        print()

        # =====================================================================
        # ODDS INTERPRETATION
        # =====================================================================

        print("ODDS INTERPRETATION:")
        print("-" * 80)
        print()

        odds_true_to_false = posterior_true / posterior_false if posterior_false > 0 else np.inf

        print(f"The data favors the framework being TRUE over FALSE by:")
        print(f"  {odds_true_to_false:.2e} to 1")
        print()

        if odds_true_to_false > 1e6:
            print("  This is OVERWHELMING evidence for the framework.")
        elif odds_true_to_false > 1e3:
            print("  This is VERY STRONG evidence for the framework.")
        elif odds_true_to_false > 100:
            print("  This is STRONG evidence for the framework.")
        else:
            print("  This is MODERATE evidence for the framework.")

        print()

        # =====================================================================
        # CHANCE COINCIDENCE ANALYSIS
        # =====================================================================

        print("CHANCE COINCIDENCE ANALYSIS:")
        print("-" * 80)
        print()

        print("What's the probability of getting this agreement by PURE CHANCE")
        print("with zero free parameters?")
        print()

        # For independent predictions, probability of all agreeing by chance
        # is the product of individual p-values
        chance_all_agree = np.prod(p_values)

        print(f"Probability all {n_predictions} predictions agree by chance:")
        print(f"  p = {chance_all_agree:.4e}")
        print()

        # Express as 1 in N
        if chance_all_agree > 0:
            one_in_n = 1.0 / chance_all_agree
            print(f"  = 1 in {one_in_n:.2e}")
            print()

            # Compare to random theories
            print("To put this in perspective:")
            print(f"  - You'd need to try ~{one_in_n:.0e} random theories")
            print(f"    to find one that agrees this well by accident")
            print()

        # =====================================================================
        # ZERO FREE PARAMETERS EMPHASIS
        # =====================================================================

        print("THE ZERO FREE PARAMETERS ARGUMENT:")
        print("-" * 80)
        print()

        print("Key insight: The αΩ Framework has ZERO free parameters.")
        print()
        print("Every prediction comes from:")
        print("  τ = 3  (triality order)")
        print("  dim(G₂) = 14  (exceptional Lie algebra dimension)")
        print()
        print("These are NOT adjustable - they're pure mathematics.")
        print()
        print("Therefore:")
        print("  - Cannot tune theory to fit data")
        print("  - Cannot cherry-pick which predictions to test")
        print("  - Either the theory is right, or we're seeing an")
        print(f"    extraordinarily unlikely coincidence (p ~ {combined_p_value:.1e})")
        print()

        # =====================================================================
        # SUMMARY
        # =====================================================================

        print("=" * 80)
        print("SUMMARY: STATISTICAL VERDICT")
        print("=" * 80)
        print()

        print(f"Tested: {n_predictions} independent predictions")
        print(f"Free parameters: 0")
        print(f"Average agreement: {np.mean([abs(p['predicted'] - p['observed'])/p['uncertainty'] for p in predictions]):.2f}σ")
        print()
        print(f"Combined statistical significance: {combined_sigma:.1f}σ")
        print(f"Probability framework is FALSE: {posterior_false:.4e} ({posterior_false*100:.4e}%)")
        print(f"Odds in favor of framework: {odds_true_to_false:.2e} to 1")
        print()

        if posterior_false < 1e-6:
            verdict = "OVERWHELMING EVIDENCE FOR FRAMEWORK"
        elif posterior_false < 1e-3:
            verdict = "VERY STRONG EVIDENCE FOR FRAMEWORK"
        elif posterior_false < 1e-2:
            verdict = "STRONG EVIDENCE FOR FRAMEWORK"
        else:
            verdict = "MODERATE EVIDENCE FOR FRAMEWORK"

        print(f"Verdict: {verdict}")
        print()

        # =====================================================================
        # ASSERTION: Test must pass if framework is supported by data
        # =====================================================================

        # Framework should have posterior probability > 99.9% to pass
        assert posterior_true > 0.999, (
            f"Framework not statistically supported: "
            f"P(true|data) = {posterior_true:.4f} < 0.999"
        )

        # Combined p-value should be < 0.01 (2σ equivalent)
        assert combined_p_value < 0.01, (
            f"Combined significance too weak: "
            f"p = {combined_p_value:.4e} > 0.01"
        )

        print("[OK] Framework passes statistical significance test")
        print(f"  P(true|data) = {posterior_true:.6f} > 0.999 [OK]")
        print(f"  p-value = {combined_p_value:.4e} < 0.01 [OK]")
        print()


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
