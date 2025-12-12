# The αΩ Physics Verification Engine

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)]()
[![Tests](https://img.shields.io/badge/tests-201%2F201%20passing-brightgreen)]()
[![Global Accuracy](https://img.shields.io/badge/accuracy-98.95%25-blue)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**A Continuous Integration (CI/CD) pipeline for the laws of physics.**

This repository treats the Standard Model of Particle Physics not as a collection of loose observations, but as a rigid software specification. It implements a Python framework that procedurally derives fundamental physical constants from a single geometric seed: the **G₂ Sedenion Algebra**.

Every number in the accompanying paper is computed at runtime. There are no magic numbers. There are no free parameters.

## The Geometric Bootstrap

The entire framework accepts exactly two integer inputs:
1.  **Triality ($\tau = 3$)**: The order of the octonion automorphism.
2.  **Dimension ($D = 14$)**: The dimension of the Lie algebra $G_2$.

From this seed `{3, 14}`, the engine derives over 50 physical observables.

### Best Results (Zero-Parameter Derivations)

| Observable | Geometric Prediction | Experimental Value | Agreement |
| :--- | :--- | :--- | :--- |
| **Top Quark Mass** | `172.8 GeV` | `172.8 GeV` | **100.00%** |
| **W Boson Mass** | `79.977 GeV` | `80.379 GeV` | **99.50%** |
| **Muon Mass** | `105.89 MeV` | `105.66 MeV` | **99.78%** |
| **Tau Mass** | `1.777 GeV` | `1.777 GeV` | **99.99%** |
| **Weak Mixing Angle** | `0.2308` | `0.2312` | **99.81%** |
| **Grand Unified Coupling** | `1/42` | `1/41.7` | **99.21%** |
| **Dark Energy ($\\Omega_\\Lambda$)** | `0.6875` | `0.6847` | **99.59%** |
| **Cabibbo Angle ($\\sin\\theta_C$)** | `0.2259` | `0.2257` | **99.92%** |

> *See the full 50+ item summary table in the generated paper.*

## The "Executable Paper" Methodology

Physics papers usually become outdated the moment new data arrives. This project is different. It is a living codebase.

1.  **Procedural Derivation**: Python modules in `physics/` calculate constants using algebraic geometry.
2.  **Unit Testing**: `pytest` validates these derivations against the latest PDG (Particle Data Group) data.
3.  **Dynamic Rendering**: The LaTeX paper is compiled on-the-fly, injecting the computed values directly into the text.

If the theory is wrong, the tests will fail.

## Usage

### 1. Setup
Initialize the verification environment:
```bash
./setup.sh
```

### 2. Run the Verification Suite
Execute the full test suite (201 tests) and compile the paper:
```bash
./run_tests.sh
```
*   **Output:** `paper/aomega_paper_single.pdf` (67 pages)
*   **Log:** Check terminal for green `[PASS]` indicators for each physical constant.

### 3. Explore the Code
*   `core/context_builder.py`: The "Kernel" that aggregates all physics.
*   `physics/`: The derivation modules.
    *   `physics/particle_fermion_masses.py`: Calculates lepton/quark masses.
    *   `physics/cosmo_dark_energy.py`: Derives $\\Omega_\\Lambda$.
*   `paper/figures/`: Generated geometric visualizations.

## Key Geometric Visualizations

The framework predicts that particle interactions are geometric operations on a $G_2$ root system.

**Figure 1: G₂ Feynman Interactions**
*   **Matter (Blue)**: Quarks and Antiquarks.
*   **Forces (Red)**: Gluons.
*   **Interactions**: Vector addition in the root space determines valid particle decays.

*(See `paper/figures/g2_color_interaction_final.pdf` after building)*



## Theoretical Note on Terminology



This framework utilizes specific definitions for algebraic properties in the Sedenion context which may differ from standard Lie theory conventions:

*   **Triality**: Refers to the **Octonionic** triality automorphism ($\\tau=3$) which is preserved by the $G_2$ symmetry group, distinct from the $D_4$ triality of Spin(8).

*   **Cubic Invariant ($Q_{eff}$)**: Often referred to as $C_3$, this denotes the effective topological winding number (degree 3) of the Sedenion associator, distinct from the standard degree-6 Casimir of pure $G_2$ Lie theory.



## Documentation

*   **[Full Physics Paper (PDF)](paper/aomega_paper_single.pdf)**: "The αΩ Theory: Zero-Parameter Theory of Everything from G₂ Geometry".
*   **[Security Whitepaper (PDF)](security/alphaomega_security_whitepaper.pdf)**: "Post-Quantum Security in a P=NP Universe".

## Authors

**Kellen Francis McNally** & **Eric Troy Sandum**
*Independent Research*

---

*"The laws of physics are the instruction set of the universe."*