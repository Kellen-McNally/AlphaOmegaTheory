#
#   Periodic table calculations without visualization
#   dependencies. This module provides the core periodic
#   table calculation functions without requiring matplotlib
#   or other plotting libraries.
#

import numpy as np
from core.constants import DIM_G2, TRIALITY, TRACE_TRIALITY as TR3_G2, RYDBERG_EV, ALPHA_EM_MEASURED

# Rydberg constant (imported from constants.py - single source of truth)
RYDBERG = RYDBERG_EV


#
#   Get principal quantum number n for element Z. Args: Z
#   (Atomic number). Returns: int (Principal quantum number
#   n, shell).
#
def get_principal_quantum_number(Z: int) -> int:
    if Z <= 2:
        return 1
    elif Z <= 10:
        return 2
    elif Z <= 18:
        return 3
    elif Z <= 36:
        return 4
    elif Z <= 54:
        return 5
    elif Z <= 86:
        return 6
    elif Z <= 118:
        return 7
    elif Z <= 168:
        return 8
    else:
        return 9


#
#   Get period number for element Z. Args: Z (Atomic
#   number). Returns: int (Period number).
#
def get_period(Z: int) -> int:
    return get_principal_quantum_number(Z)


#
#   Number of elements in a given period. Args: period
#   (Period number 1-9). Returns: int (Number of elements in
#   that period).
#
def elements_in_period(period: int) -> int:
    period_elements = {
        1: 2,   # H, He
        2: 8,   # Li-Ne
        3: 8,   # Na-Ar
        4: 18,  # K-Kr
        5: 18,  # Rb-Xe
        6: 32,  # Cs-Rn
        7: 32,  # Fr-Og
        8: 50,  # 119-168 (g-block begins)
        9: 50   # 169-218
    }
    return period_elements.get(period, 0)


#
#   Get the period that element Z belongs to. Args: Z
#   (Atomic number). Returns: int (Period number).
#
def element_period(Z: int) -> int:
    return get_period(Z)


#
#   Return list of noble gas atomic numbers. Returns: list
#   (Atomic numbers of noble gases).
#
def noble_gas_list() -> list:
    return [2, 10, 18, 36, 54, 86, 118]


#
#   Number of lanthanide elements. Returns: int (Count of
#   lanthanides, 15).
#
def lanthanide_count() -> int:
    return 15


#
#   Number of actinide elements. Returns: int (Count of
#   actinides, 15).
#
def actinide_count() -> int:
    return 15


#
#   Summary statistics about the periodic table structure.
#   Returns: dict (Periodic table structure info from G₂
#   geometry).
#
def periodic_table_summary() -> dict:
    return {
        "total_periods": 9,
        "elements_known": 118,
        "elements_predicted": 208,
        "structure": {
            "s_block": "2 columns",
            "p_block": "6 columns",
            "d_block": "10 columns (transition metals)",
            "f_block": "14 columns (lanthanides, actinides)",
            "g_block": "18 columns (period 8+, predicted)"
        },
        "g2_origin": "G₂ Lie algebra with dim=14, triality=3",
        "shell_formula": "2n² electrons per shell",
        "noble_gases": noble_gas_list(),
        "periods": {
            1: elements_in_period(1),
            2: elements_in_period(2),
            3: elements_in_period(3),
            4: elements_in_period(4),
            5: elements_in_period(5),
            6: elements_in_period(6),
            7: elements_in_period(7),
            8: elements_in_period(8),
            9: elements_in_period(9)
        }
    }


#
#   Calculate effective nuclear charge using G₂ geometry.
#   Args: Z (Atomic number), n (Principal quantum number,
#   optional). Returns: float (Effective nuclear charge
#   Z_eff).
#
def calculate_zeff_g2(Z: int, n: int = None) -> float:
    if n is None:
        n = get_principal_quantum_number(Z)

    # Special cases
    if Z == 1:
        return 1.0
    if Z == 2:
        return 1.34  # Measured value

    # Slater's screening constant with G₂ correction
    # G₂ correction factor: (Tr₃ - dim) / Tr₃ = (17 - 14) / 17 = 3/17
    g2_correction = (TR3_G2 - DIM_G2) / TR3_G2

    # Simple Slater approximation
    S = 0.30 * (Z - 1)  # Screening by other electrons

    # Apply G₂ correction
    Z_eff = Z - S * (1 - g2_correction)

    # Physical bounds
    if Z_eff < 1:
        Z_eff = 1.0
    if Z_eff > Z:
        Z_eff = float(Z)

    return Z_eff


#
#   Calculate ionization energy from G₂-derived Z_eff. Args:
#   Z (Atomic number), Z_eff (Effective nuclear charge,
#   optional), n (Principal quantum number, optional).
#   Returns: float (Ionization energy in eV).
#
def calculate_ionization_energy(Z: int, Z_eff: float = None, n: int = None) -> float:
    if n is None:
        n = get_principal_quantum_number(Z)
    if Z_eff is None:
        Z_eff = calculate_zeff_g2(Z, n)

    # Base ionization energy
    IE = RYDBERG * (Z_eff ** 2) / (n ** 2)

    # Relativistic correction for heavy elements
    if Z > 36:
        alpha_fine = ALPHA_EM_MEASURED
        rel_correction = 1.0 + (alpha_fine * Z) ** 2
        IE *= rel_correction

    return IE


#
#   Calculate ionization energies for multiple elements.
#   Args: Z_list (List of atomic numbers, optional).
#   Returns: dict (Summary of ionization energies with
#   predictions and comparisons).
#
def ionization_energy_summary(Z_list: list = None) -> dict:
    if Z_list is None:
        Z_list = list(range(1, 21))

    results = {}
    for Z in Z_list:
        n = get_principal_quantum_number(Z)
        Z_eff = calculate_zeff_g2(Z, n)
        IE = calculate_ionization_energy(Z, Z_eff, n)

        results[Z] = {
            "n": n,
            "Z_eff": Z_eff,
            "IE_eV": IE
        }

    return results
