#!/usr/bin/env python3
"""
New Paper Build System

Replaces the old scattered paper compilation with a clean,
centralized approach using a unified context and structure definition.
"""

import sys
import os
import re
import argparse
from pathlib import Path
from typing import List

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# New Architecture Imports
from core.context_builder import build_paper_context
from paper.structure import SECTIONS
from paper.template_engine import render_template
from paper.paper_api import PaperSection
from paper.citations import citation_db, initialize_attributions

def _escape_unicode_for_latex(text: str) -> str:
    """Escape Unicode characters to LaTeX commands."""
    replacements = {
        "α": r"$\alpha$",
        "β": r"$\beta$",
        "γ": r"$\gamma$",
        "δ": r"$\delta$",
        "Λ": r"$\Lambda$",
        "Δ": r"$\Delta$",
        "μ": r"$\mu$",
        "ν": r"$\nu$",
        "π": r"$\pi$",
        "ρ": r"$\rho$",
        "σ": r"$\sigma$",
        "τ": r"$\tau$",
        "Ω": r"$\Omega$",
        "Ω": r"$\Omega$",
        "θ": r"$\theta$",
        "−": r"-",
        "×": r"$\times$",
        "±": r"$\pm$",
        "°": r"$^\circ$",
        "É": r"'E",
        "ℏ": r"$\hbar$",
        "₀": r"$_0$",
        "₁": r"$_1$",
        "₂": r"$_2$",
        "₃": r"$_3$",
        "₄": r"$_4$",
    }
    for char, replacement in replacements.items():
        text = text.replace(char, replacement)
    return text


def _wrap_subsections(content: str) -> str:
    """
    Wraps subsection titles in \texorpdfstring to fix hyperref warnings.
    Matches \subsection{Title} and replaces with \subsection{\texorpdfstring{Title}{CleanTitle}}.
    """
    def repl(match):
        cmd = match.group(1) # subsection or subsubsection
        title = match.group(2)
        
        # Generate clean title for PDF
        clean = title
        clean = re.sub(r'\$([^\$]+)\$', r'\1', clean) # Remove $..$ but keep content
        clean = clean.replace('\\', '')
        clean = clean.replace('{', '').replace('}', '')
        clean = clean.replace('^', '')
        clean = clean.replace('_', '')
        
        # Replace symbols
        replacements = {
            "alpha": "Alpha", "beta": "Beta", "gamma": "Gamma", "zeta": "Zeta",
            "pi": "Pi", "tau": "Tau", "sigma": "Sigma", "lambda": "Lambda",
            "mu": "Mu", "nu": "Nu", "theta": "Theta", "delta": "Delta",
            "infty": "inf", "to": "->", "approx": "~",
            "neq": "!=", "le": "<=", "ge": ">=",
        }
        for k, v in replacements.items():
            clean = re.sub(r'\\b' + k + r'\\b', v, clean)
            
        # Use format to avoid f-string backslash hell
        # We want: \subsection{\texorpdfstring{Title}{Clean}}
        # Correctly formatted string for output (single backslashes for commands)
        return r"\{0}{{\texorpdfstring{{{1}}}{{{2}}}}}".format(cmd, title, clean)

    # Match \subsection{...} or \subsubsection{...}
    pattern = r'\\((?:sub)*section)\{((?:[^\{\}]|\{[^\}]*\})*)\}'
    return re.sub(pattern, repl, content)


def _latex_references() -> str:
    """Generate references section in Physical Review style."""
    latex = r"""
\begin{thebibliography}{99}

"""
    # Sort citations by year
    sorted_cites = sorted(citation_db.citations.values(), key=lambda c: c.year)

    for cite in sorted_cites:
        if len(cite.authors) > 1:
            authors_formatted = (
                ", ".join(cite.authors[:-1]) + ", and " + cite.authors[-1]
            )
        else:
            authors_formatted = cite.authors[0]

        title_escaped = _escape_unicode_for_latex(cite.title).replace("&", r"\&")
        journal_escaped = _escape_unicode_for_latex(cite.journal).replace("&", r"\&")

        latex += r"\bibitem{"+ cite.key + "}\n"
        latex += f"{authors_formatted},\n"
        latex += r"\textit{"+ title_escaped + "},\n"
        latex += f"{journal_escaped} ({cite.year}).\n\n"

    latex += r"""
\end{thebibliography}

"""
    return latex


def compile_latex(output_filename: str = "aomega_paper.tex") -> str:
    """
    Compile all sections into a single LaTeX document using unified context.
    """
    # Initialize citations
    initialize_attributions()

    # 1. Build Global Context
    print("Building unified paper context...")
    context = build_paper_context()
    print(f"Context loaded with {len(context)} keys.")

    # LaTeX document header (read from template)
    header_path = Path(project_root) / "paper" / "templates" / "header.tex"
    latex_content = header_path.read_text(encoding="utf-8")

    # 2. Dynamic Abstract
    abstract_content = render_template("abstract.tex", context)
    latex_content += abstract_content

    # Footer for abstract (part of footer.tex)
    footer_path = Path(project_root) / "paper" / "templates" / "footer.tex"
    footer_content = footer_path.read_text(encoding="utf-8")

    # Define part boundaries
    part_boundaries = {
        "foundation": "Part I: Mathematical Foundation",
        "triality": "Part II: Core Physical Predictions",
        "temporal_antimatter": "Part III: Extended Predictions",
        "proton_decay": r"Part IV: Cosmology \& Gravity",
        "complete_action": "Part V: Meta-Analysis and Validation",
    }

    # 3. Add each section
    for section_def in SECTIONS:
        sec_id = section_def["id"]
        title = section_def["title"]
        template = section_def["template"]
        
        # Add part header
        if sec_id in part_boundaries:
            latex_content += "\n"
            latex_content += r"\vspace{1cm}" + "\n"
            latex_content += r"\begin{center}" + "\n"
            latex_content += (
                r"\Large\textbf{"+ part_boundaries[sec_id] + r"}" + "\n"
            )
            latex_content += r"\end{center}" + "\n"
            latex_content += r"\vspace{0.5cm}" + "\n\n"

        # Inject \appendix command
        if sec_id == "mathematical_appendices":
            latex_content += r"\appendix" + "\n"

        # Render Content
        raw_content = render_template(template, context)

        # Escape unicode characters in content
        processed_content = _escape_unicode_for_latex(raw_content)
        
        # Wrap subsections in texorpdfstring
        processed_content = _wrap_subsections(processed_content)

        # Title processing
        title_latex = _escape_unicode_for_latex(title)
        
        # Clean title for PDF bookmarks
        title_pdf = title
        pdf_replacements = {
            "α": "Alpha", "β": "Beta", "γ": "Gamma", "δ": "Delta",
            "Ω": "Omega", "Λ": "Lambda", "τ": "Tau", "θ": "Theta",
            "μ": "Mu", "ν": "Nu", "π": "Pi", "G₂": "G2",
            "²": "2", "³": "3", "₀": "0", "₁": "1"
        }
        for k, v in pdf_replacements.items():
            title_pdf = title_pdf.replace(k, v)
        title_pdf = title_pdf.replace("$", "").replace("\\", "")

        # Generate label
        slug = re.sub(r"[^a-zA-Z0-9]+", "_", title.lower()).strip("_")
        label = f"sec:{slug}"

        latex_content += "\n"
        # Correctly formatted section command with single backslashes
        latex_content += r"\section{\texorpdfstring{"+ title_latex + r"}{"+ title_pdf + r"}}" + "\n"
        latex_content += r"\label{"+ label + r"}" + "\n\n"
        latex_content += processed_content + "\n"

    # 4. Add references
    latex_content += _latex_references()

    # 5. Add final footer
    latex_content += footer_content

    # Write to file
    output_path = Path(__file__).parent / output_filename
    with open(output_path, "w") as f:
        f.write(latex_content)

    return latex_content

def build_paper(format_type: str = "single") -> str:
    """Build the paper in specified format."""
    print(f"Building αΩ Framework paper ({format_type}-column)...")

    # Compile LaTeX
    filename = f"aomega_paper_{format_type}.tex"
    latex_content = compile_latex(filename)

    print(f"Generated {filename}")
    print(f"LaTeX document: {len(latex_content)} characters")

    return str(Path(__file__).parent / filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build αΩ Framework paper"
    )
    parser.add_argument("--latex", action="store_true", help="Generate LaTeX files")
    parser.add_argument("--markdown", action="store_true", help="Generate Markdown (Placeholder)")

    args = parser.parse_args()

    if not args.latex and not args.markdown:
        args.latex = True

    if args.latex:
        single_path = build_paper("single")
        print("\nPaper compilation complete!")
        print(f"Single-column: {single_path}")
