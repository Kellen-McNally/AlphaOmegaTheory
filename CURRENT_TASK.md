# Current Task: Restore Comprehensive Summary Table in Paper Generation

## Context
The αΩ Framework project uses a Python-based pipeline (`paper/build.py`) to generate a LaTeX paper (`paper/aomega_paper_single.pdf`) from the underlying physics code. A critical component of this paper is the "Comprehensive Summary" table, which must list over 50 derived physical observables to demonstrate the framework's breadth.

## Problem
The summary table is currently broken or truncated. Previous attempts to fix it by modifying `core/context_builder.py` (the script responsible for aggregating data for the LaTeX templates) have introduced regression errors:
1.  **Missing Data Keys**: Code blocks responsible for formatting Lepton and Quark masses (e.g., `electron_mass_pred_mev`) were accidentally deleted during file rewrites, causing `KeyError`s during build.
2.  **Syntax Errors**: Manual string construction for the table rows resulted in Python syntax errors (unterminated strings, bad escaping of `%` characters).
3.  **Truncated Output**: Prior to these errors, the table only showed ~19 entries instead of the required >50.

## Objective
The immediate goal is to **restore the full >50 item summary table** and ensure the paper compiles successfully.

### Action Plan
1.  **Repair `core/context_builder.py`**:
    *   Re-insert the missing code sections for "Lepton Formatting" and "Quark Formatting" (converting values to MeV/GeV and calculating error percentages).
    *   Correct the "Higgs VEV" row which is currently missing its numerical error value in the formatted string.
    *   Verify all LaTeX special characters (like `%`) are properly escaped as `\\%` in the Python strings.
2.  **Verify Data Generation**:
    *   Execute `core/context_builder.py` directly as a script to confirm it runs without `KeyError` or `SyntaxError`.
    *   Verify it returns the expected number of keys (approx 166+).
3.  **Build Paper**:
    *   Run `./run_tests.sh` to compile the LaTeX document.
4.  **Verification**:
    *   Confirm the PDF is generated.
    *   Confirm the summary table covers all physics sectors (Constants, Cosmology, Gauge, Leptons, Quarks, Neutrinos, CKM, Atomic, Exotics).

## Status
*   [ ] `core/context_builder.py` fixed (missing blocks restored).
*   [ ] Context generation verification passed.
*   [ ] Paper build passed.
