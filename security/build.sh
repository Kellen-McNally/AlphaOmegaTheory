#!/bin/bash
# Build the Security Whitepaper

echo "Building AlphaOmega Security Whitepaper..."
pdflatex -interaction=nonstopmode alphaomega_security_whitepaper.tex
pdflatex -interaction=nonstopmode alphaomega_security_whitepaper.tex

echo "Done. Output: security/alphaomega_security_whitepaper.pdf"
