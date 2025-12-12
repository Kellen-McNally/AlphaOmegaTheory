#!/bin/bash
################################################################################
# αΩ Framework Setup Script
#
# This script sets up the Python virtual environment and installs all
# dependencies needed to run tests and generate the paper.
#
# Usage:
#   bash setup.sh
#
################################################################################

set -e  # Exit on any error

echo "════════════════════════════════════════════════════════════════"
echo "  αΩ Framework Setup"
echo "════════════════════════════════════════════════════════════════"
echo ""

# Navigate to script directory
cd "$(dirname "$0")"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check Python version
echo -e "${BLUE}[1/5]${NC} Checking Python version..."
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}[FAIL] python3 could not be found${NC}"
    exit 1
fi
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
echo "      Found: Python $PYTHON_VERSION"

# Check for python3-venv
echo ""
echo -e "${BLUE}[2/5]${NC} Checking for python3-venv package..."
if ! python3 -m venv --help > /dev/null 2>&1; then
    echo -e "${RED}[FAIL] python3-venv is not installed${NC}"
    echo ""
    echo "Please install python3-venv (e.g., sudo apt install python3-venv)"
    exit 1
fi
echo -e "      ${GREEN}[OK] python3-venv is available${NC}"

# Create/Update virtual environment
echo ""
echo -e "${BLUE}[3/5]${NC} Setting up virtual environment (.venv)..."
if [ ! -d ".venv" ]; then
    python3 -m venv .venv
    echo -e "      ${GREEN}[OK] Virtual environment created${NC}"
else
    echo -e "      ${GREEN}[OK] Virtual environment already exists${NC}"
fi

# Activate virtual environment
echo ""
echo -e "${BLUE}[4/5]${NC} Activating virtual environment..."
source .venv/bin/activate
echo -e "      ${GREEN}[OK] Virtual environment activated${NC}"

# Upgrade pip and install dependencies
echo ""
echo -e "${BLUE}[5/5]${NC} Installing dependencies..."
pip install --quiet --upgrade pip

if [ -f "pyproject.toml" ]; then
    pip install --quiet -e .[dev]
    echo -e "      ${GREEN}[OK] Project and dev dependencies installed (pyproject.toml)${NC}"
elif [ -f "requirements.txt" ]; then
    pip install --quiet -r requirements.txt
    if [ -f "requirements-test.txt" ]; then
        pip install --quiet -r requirements-test.txt
    fi
    echo -e "      ${GREEN}[OK] Dependencies installed from requirements files${NC}"
else
    echo -e "      ${YELLOW}[WARN] No dependency file found. Installing core packages...${NC}"
    pip install --quiet numpy scipy pytest matplotlib sympy
    echo -e "      ${GREEN}[OK] Core dependencies installed${NC}"
fi

echo ""
echo "════════════════════════════════════════════════════════════════"
echo -e "  ${GREEN}[OK] Setup Complete!${NC}"
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "To activate the environment manually, run:"
echo -e "  ${GREEN}source .venv/bin/activate${NC}"
echo ""