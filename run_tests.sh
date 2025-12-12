#!/bin/bash
################################################################################
# αΩ Framework: Master Build and Test Script
#
# ** THIS IS THE ENTRY POINT FOR ALL OPERATIONS **
#
# What this script does:
#   1. Runs ALL tests (validates code correctness)
#   2. Generates predictions summary (if tests pass)
#   3. Compiles paper from code (if tests pass)
#
# Why use this instead of running scripts directly?
#   - Ensures paper reflects PASSING, VALIDATED code
#   - Prevents drift between paper and implementation
#   - Maintains contract that paper = working code
#
# Usage:
#   ./run_tests.sh --test   # Run tests + generate paper (default)
#   ./run_tests.sh --paper  # Generate paper only (skip tests, use cached results)
#   ./run_tests.sh          # Same as --test (backward compatible)
#
# Advanced options:
#   ./run_tests.sh --test --thermodynamic_samples 1000000000  # Production mode
#
# DO NOT run individual scripts like:
#   python paper/compiler.py    # [FAIL] WRONG (bypasses tests)
#   pytest api/                 # [FAIL] WRONG (no paper generation)
#
# See: .claude-workflow and CONTRIBUTING.md for details
################################################################################

cd "$(dirname "$0")"

# Parse command line arguments
MODE="test"  # Default: run tests + paper
THERMODYNAMIC_SAMPLES="100000"  # Default: fast demo

while [[ $# -gt 0 ]]; do
    case $1 in
        --test)
            MODE="test"
            shift
            ;;
        --paper)
            MODE="paper"
            shift
            ;;
        --thermodynamic_samples)
            THERMODYNAMIC_SAMPLES="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--test|--paper] [--thermodynamic_samples N]"
            echo ""
            echo "Modes:"
            echo "  --test   Run tests + generate paper (default)"
            echo "  --paper  Generate paper only (skip tests)"
            exit 1
            ;;
    esac
done

# Export for use by paper generation
export AOMEGA_THERMODYNAMIC_SAMPLES="$THERMODYNAMIC_SAMPLES"

# Activate venv
source .venv/bin/activate

# Initialize test exit code
TEST_EXIT=0

# Run tests if in test mode
if [ "$MODE" = "test" ]; then
    echo "=========================================="
    echo "αΩ Framework Complete Test Suite"
    echo "=========================================="
    echo ""
    echo "Mode: TEST (running tests + generating paper)"
    echo "Thermodynamic samples: ${THERMODYNAMIC_SAMPLES}"
    echo ""

    echo "Running all tests (excluding octonions which needs torch)..."
    echo ""

    # Run tests with clear output
    # Note: Test summary is now handled by pytest hook in tests/conftest.py
    # which displays key formulas, statistical significance, and results
    .venv/bin/python -m pytest tests/ \
        --ignore=tests/test_core_octonions.py \
        -v \
        --tb=line

    # Capture test exit code
    TEST_EXIT=$?
    echo ""
fi

# Generate paper if tests passed (or if in paper-only mode)
if [ $TEST_EXIT -eq 0 ]; then
    if [ "$MODE" = "paper" ]; then
        echo "=========================================="
        echo "αΩ Framework Paper Generation"
        echo "=========================================="
        echo ""
        echo "Mode: PAPER (using cached test results)"
        echo "Thermodynamic samples: ${THERMODYNAMIC_SAMPLES}"
        echo ""
    fi
    echo "Generating Paper from Code"
    echo "=========================================="
    echo ""

    # Generate README.md from code
    # Set environment variable so compiler knows it was called from run_tests.sh
    FROM_RUN_TESTS=1 .venv/bin/python paper/build.py --markdown
    echo ""

    # Generate LaTeX
    echo "Generating LaTeX..."
    export FROM_RUN_TESTS=1
    .venv/bin/python paper/build.py --latex
    
    # Compile PDF
    echo "Compiling PDF..."
    cd paper
    pdflatex -interaction=nonstopmode aomega_paper_single.tex > /dev/null
    bibtex aomega_paper_single > /dev/null
    pdflatex -interaction=nonstopmode aomega_paper_single.tex > /dev/null
    pdflatex -interaction=nonstopmode aomega_paper_single.tex
    cd ..
else
    echo ""
    echo "=========================================="
    echo "[FAIL] Tests Failed - Paper NOT Generated"
    echo "=========================================="
    echo ""
    echo "Paper generation skipped because tests did not pass."
    echo "Fix the failing tests and run again:"
    echo "  ./run_tests.sh --test"
    echo ""
fi

exit $TEST_EXIT
