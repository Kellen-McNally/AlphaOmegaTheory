"""
Paper Section API

[WARN] This defines the API for modules to contribute sections to the paper.
[WARN] Paper is generated via: ./run_tests.sh (not by running this file)

Standard interface for modules to contribute to the auto-generated paper.

KEY PRINCIPLE: All numbers in the paper text are computed procedurally.
Never hardcode numerical values - always compute them from the canonical functions.

Each module can define a paper_section() function that returns a dict with:
  - title: Section title
  - content: LaTeX content (can use f-strings with computed values)
  - subsections: Optional list of subsection dicts
  - equations: Key equations (will be formatted consistently)
  - results: Numerical results with comparisons to experiment
  - order: Integer ordering (lower = earlier in paper)
"""

from typing import List, Optional, Callable, Any
from dataclasses import dataclass, field


@dataclass
class PaperEquation:
    """Represents a key equation in the paper."""

    label: str  # LaTeX label (e.g., "alpha_gut")
    equation: str  # LaTeX equation (e.g., r"\alpha_{GUT} = \frac{1}{42}")
    description: str  # Plain text description
    exact_form: Optional[str] = None  # Exact form (e.g., "1/42")
    compute: Optional[Callable[[], Any]] = None  # Function to compute the value


@dataclass
class PaperResult:
    """
    Represents a prediction/result with experimental comparison.

    IMPORTANT: predicted value is COMPUTED, not hardcoded.
    """

    name: str  # Observable name
    symbol: str  # LaTeX symbol
    compute_predicted: Callable[[], float]  # Function that computes predicted value
    exact_form: str  # Exact form (e.g., "11/16")
    observed: Optional[float] = None  # Experimental value
    uncertainty: Optional[float] = None  # Experimental uncertainty
    units: str = ""  # Units (if applicable)
    status: str = "prediction"  # "prediction", "postdiction", "explanation"

    @property
    def predicted(self) -> float:
        """Compute predicted value on demand."""
        return self.compute_predicted()

    @property
    def error_percent(self) -> Optional[float]:
        """Compute error percentage if observed value available."""
        if self.observed is None:
            return None
        return abs(self.predicted - self.observed) / self.observed * 100


@dataclass
class PaperSection:
    """
    Standard paper section structure.

    All modules should return this format from their paper_section() function.

    CRITICAL: content can be a callable that returns LaTeX string with
    all numbers computed procedurally.
    """

    title: str
    content: str | Callable[[], str]  # LaTeX content or function that generates it
    order: int  # Section ordering (10, 20, 30, ...)

    # Optional fields
    subsections: List["PaperSection"] = field(default_factory=list)
    equations: List[PaperEquation] = field(default_factory=list)
    results: List[PaperResult] = field(default_factory=list)
    figures: List[str] = field(default_factory=list)  # Figure filenames
    citations: List[str] = field(default_factory=list)  # BibTeX keys

    def get_content(self) -> str:
        """Get content, computing it if it's a callable."""
        if callable(self.content):
            return self.content()
        return self.content

    def to_latex(self, level: int = 1) -> str:
        """
        Convert to LaTeX format.

        Args:
            level: Section level (1=section, 2=subsection, 3=subsubsection)

        Returns:
            str: LaTeX formatted section
        """
        if level == 1:
            section_cmd = "section"
        elif level == 2:
            section_cmd = "subsection"
        elif level == 3:
            section_cmd = "subsubsection"
        else:
            section_cmd = "paragraph"

        latex = f"\\{section_cmd}{{{self.title}}}\n"
        latex += f"\\label{{sec:{self.title.lower().replace(' ', '_')}}}\n\n"

        # Get content (compute if callable)
        latex += self.get_content() + "\n\n"

        # Add equations if present
        if self.equations:
            for eq in self.equations:
                latex += "\\begin{equation}\n"
                latex += f"  {eq.equation}\n"
                latex += f"  \\label{{eq:{eq.label}}}\n"
                latex += "\\end{equation}\n"
                latex += f"{eq.description}"

                # Add computed value if available
                if eq.compute is not None:
                    computed_value = eq.compute()
                    latex += f" (numerically: ${computed_value:.10f}$)"

                latex += "\n\n"

        # Add results table if present
        if self.results:
            latex += self._format_results_table()

        # Add subsections
        for subsec in self.subsections:
            latex += subsec.to_latex(level=level + 1)

        return latex

    def _format_results_table(self) -> str:
        """Format results as a LaTeX table with computed values."""
        latex = "\\begin{table}[ht]\n"
        latex += "\\centering\n"
        latex += "\\begin{tabular}{lcccc}\n"
        latex += "\\hline\n"
        latex += "Observable & Predicted & Exact & Observed & Error \\\\\n"
        latex += "\\hline\n"

        for result in self.results:
            obs_str = f"${result.symbol}$"

            # Compute predicted value procedurally
            predicted_val = result.predicted
            pred_str = f"{predicted_val:.6g}"
            exact_str = f"${result.exact_form}$"

            if result.observed is not None:
                obs_val_str = f"{result.observed:.6g}"
                if result.uncertainty is not None:
                    obs_val_str += f" $\\pm$ {result.uncertainty:.6g}"

                # Compute error procedurally
                error_pct = result.error_percent
                error_str = f"{error_pct:.2f}\\%"
            else:
                obs_val_str = "---"
                error_str = "---"

            latex += f"{obs_str} & {pred_str} & {exact_str} & {obs_val_str} & {error_str} \\\\\n"

        latex += "\\hline\n"
        latex += "\\end{tabular}\n"
        latex += f"\\caption{{Predictions from {self.title}}}\n"
        latex += "\\end{table}\n\n"

        return latex

    def to_markdown(self, level: int = 1) -> str:
        """
        Convert to Markdown format (for README).

        Args:
            level: Heading level

        Returns:
            str: Markdown formatted section
        """
        header = "#" * level
        md = f"{header} {self.title}\n\n"

        # Get content (compute if callable)
        content = self.get_content()

        # Convert LaTeX to plain markdown (simple conversion)
        content = content.replace(r"\textbf{", "**").replace("}", "**")
        content = content.replace(r"\textit{", "*").replace("}", "*")
        content = content.replace(r"$", "`")

        md += content + "\n\n"

        # Add equations
        if self.equations:
            md += "**Key Equations:**\n\n"
            for eq in self.equations:
                eq_text = f"- {eq.description}: `{eq.exact_form or eq.equation}`"
                if eq.compute is not None:
                    computed = eq.compute()
                    eq_text += f" = {computed:.10g}"
                md += eq_text + "\n"
            md += "\n"

        # Add results
        if self.results:
            md += "**Results:**\n\n"
            md += "| Observable | Predicted | Observed | Error |\n"
            md += "|------------|-----------|----------|-------|\n"

            for result in self.results:
                # Compute values procedurally
                pred_val = result.predicted
                obs_str = f"{result.observed:.6g}" if result.observed else "---"
                error_str = (
                    f"{result.error_percent:.2f}%" if result.error_percent else "---"
                )

                md += f"| {result.name} | {pred_val:.6g} | {obs_str} | {error_str} |\n"
            md += "\n"

        # Add subsections
        for subsec in self.subsections:
            md += subsec.to_markdown(level=level + 1)

        return md


def create_section(
    title: str, content: str | Callable[[], str], order: int, **kwargs
) -> PaperSection:
    """
    Convenience function to create a PaperSection.

    Args:
        title: Section title
        content: LaTeX content (string or callable that returns string)
        order: Section ordering
        **kwargs: Additional PaperSection fields

    Returns:
        PaperSection instance
    """
    return PaperSection(title=title, content=content, order=order, **kwargs)


# ============================================================================
# EXAMPLE: How to implement paper_section() in a module
# ============================================================================


def example_paper_section() -> PaperSection:
    """
    Example of how a module should implement paper_section().

    KEY: All numbers are computed procedurally, never hardcoded!
    """
    # Import the canonical functions
    from core.constants import ALPHA_GUT, TRIALITY, DIM_G2

    # Define content generation function that computes all values
    def generate_content() -> str:
        # Compute values procedurally
        alpha_value = ALPHA_GUT
        tau = TRIALITY
        dim = DIM_G2

        # Generate LaTeX with computed values
        product = tau * dim
        return rf"""
The grand unified coupling is predicted exactly from G₂ geometry.
The triality order $\tau = {tau}$ and dimension $\dim(G_2) = {dim}$
are fixed by the Lie algebra structure, giving
$\alpha_{{GUT}} = 1/{product} = {alpha_value:.10f}$.

No free parameters are required. This is pure geometry.
"""

    return PaperSection(
        title="Grand Unified Coupling",
        order=10,
        content=generate_content,  # Pass the function, not the result!
        equations=[
            PaperEquation(
                label="alpha_gut",
                equation=r"\alpha_{GUT} = \frac{1}{\tau \times \dim(G_2)} = \frac{1}{42}",
                description="Grand unified coupling from G₂ triality and dimension",
                exact_form="1/42",
                compute=lambda: ALPHA_GUT,  # Compute procedurally
            )
        ],
        results=[
            PaperResult(
                name="Grand Unified Coupling",
                symbol=r"\alpha_{GUT}",
                compute_predicted=lambda: ALPHA_GUT,  # COMPUTED, not hardcoded!
                exact_form="1/42",
                observed=None,
                status="prediction",
            )
        ],
    )


# ============================================================================
# EXAMPLE: More complex content with multiple computed values
# ============================================================================


def example_complex_section() -> PaperSection:
    """Example with multiple computed values in text."""
    from core.constants import OMEGA_LAMBDA, CASIMIR_C3_G2, DIM_G2, RANK_G2
    from core.constants import OMEGA_LAMBDA_EXP, OMEGA_LAMBDA_ERR

    def generate_content() -> str:
        # Compute all values
        omega_pred = OMEGA_LAMBDA
        omega_obs = OMEGA_LAMBDA_EXP
        omega_err = OMEGA_LAMBDA_ERR
        c3 = CASIMIR_C3_G2
        dim = DIM_G2
        rank = RANK_G2

        # Compute error
        error_pct = abs(omega_pred - omega_obs) / omega_obs * 100

        return rf"""
The dark energy density is determined by the G₂ Casimir invariants:
$$\Omega_\Lambda = \frac{{C_3}}{{\dim + \text{{rank}}}} = \frac{{{c3}}}{{{dim} + {rank}}} = {omega_pred:.6f}$$

This prediction agrees with Planck 2018 measurements:
$\Omega_\Lambda = {omega_obs:.4f} \pm {omega_err:.4f}$
with error {error_pct:.2f}\%.
"""

    return PaperSection(
        title="Dark Energy Density",
        order=20,
        content=generate_content,
        results=[
            PaperResult(
                name="Dark Energy Density",
                symbol=r"\Omega_\Lambda",
                compute_predicted=lambda: OMEGA_LAMBDA,
                exact_form="11/16",
                observed=OMEGA_LAMBDA_EXP,
                uncertainty=OMEGA_LAMBDA_ERR,
                status="postdiction",
            )
        ],
    )
