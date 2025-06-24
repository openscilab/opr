# -*- coding: utf-8 -*-
"""OPR functions."""
from typing import Tuple
import math
from .params import A_WEIGHT, T_WEIGHT, C_WEIGHT, G_WEIGHT, ANHYDROUS_MOLECULAR_WEIGHT_CONSTANT
from .params import BASE_EXTINCTION_COEFFICIENTS, NN53_EXTINCTION_COEFFICIENTS
from .params import NN_PARAMS, DNA_COMPLEMENT_MAP

def molecular_weight_calc(sequence: str) -> float:
    """
    Calculate molecular weight and return it.

    :param sequence: primer nucleotides sequence
    """
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    return (a_count * A_WEIGHT) + (t_count * T_WEIGHT) + (c_count * C_WEIGHT) + \
        (g_count * G_WEIGHT) - ANHYDROUS_MOLECULAR_WEIGHT_CONSTANT


def basic_melting_temperature_calc(sequence: str) -> float:
    """
    Calculate basic melting temperature and return it.

    :param sequence: primer nucleotides sequence
    """
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    if len(sequence) <= 13:
        melting_temperature = (a_count + t_count) * 2 + (g_count + c_count) * 4
    else:
        melting_temperature = 64.9 + 41 * ((g_count + c_count - 16.4) / (a_count + t_count + g_count + c_count))
    return melting_temperature


def salt_adjusted_melting_temperature_calc(sequence: str, salt: float) -> float:
    """
    Calculate the salt-adjusted melting temperature (Tm) of a primer sequence and return it.

    :param sequence: Primer nucleotides sequence
    :param salt: Sodium ion concentration in moles (unit mM)
    """
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    seq_length = len(sequence)
    if seq_length <= 13:
        salt_adjustment = 16.6 * (math.log10(salt) - 3) - 16.6 * math.log10(0.050)
        tm = (a_count + t_count) * 2 + (g_count + c_count) * 4 + salt_adjustment
    else:
        tm = (
            100.5 + (41 * (g_count + c_count) / seq_length)
            - (820 / seq_length)
            + 16.6 * (math.log10(salt) - 3)
        )
    return tm


def calculate_thermodynamics_constants(sequence: str) -> Tuple[float, float]:
    """
    Calculate ΔH (in kcal/mol) and ΔS (in kcal/mol·K) for a primer sequence and return it.

    :param sequence: Primer nucleotides sequence
    """
    delta_h = 0.0
    delta_s = 0.0

    # Sum ΔH and ΔS from nearest neighbors
    for i in range(len(sequence) - 1):
        pair = sequence[i:i+2]
        if pair in NN_PARAMS:
            dh, ds = NN_PARAMS[pair]
        else:
            rev_comp_pair = ''.join(DNA_COMPLEMENT_MAP[base] for base in reversed(pair))
            dh, ds = NN_PARAMS[rev_comp_pair]
        delta_h += dh
        delta_s += ds
    return delta_h, delta_s


def nearest_neighbor_melting_temperature_calc(sequence: str, na_salt: float, thermodynamic_constants: Tuple[float, float]) -> float:
    """
    Calculate the Nearest neighbor melting temperature (Tm) of a primer sequence and return it.

    :param sequence: Primer nucleotides sequence
    :param na_salt: Sodium ion concentration in millimoles (unit mM)
    :param thermodynamic_constants: (ΔH, ΔS)
    """
    # Convert salt from mM to M
    na_conc = na_salt / 1000.0
    # Ensure uppercase sequence
    sequence = sequence.upper()
    delta_h , delta_s = thermodynamic_constants

    # Constants for Nearest Neighbors formula
    A = -0.0108         # kcal / (K·mol), helix initiation constant
    R = 0.00199         # kcal / (K·mol), gas constant
    C = 0.5e-6          # mol/L = 0.5 µM oligonucleotide concentration

    # Compute melting temperature using the formula
    denominator = A + delta_s + (R * math.log(C / 4))
    tm = (delta_h / denominator) - 273.15 + (16.6 * math.log10(na_conc))
    return tm


def gc_clamp_calc(sequence: str) -> int:
    """
    Calculate GC clamp, number of guanine (G) or cytosine (C) bases in the last 5 bases of the primer, and return it.

    :param sequence: primer sequence
    """
    if len(sequence) < 5:
        return 0
    return sequence[-5:].count('G') + sequence[-5:].count('C')


def e260_ssnn_calc(sequence: str) -> float:
    """
    Calculate the extinction coefficient for a primer and return it.

    It uses nearest-neighbor model for single strand DNA (ss-nn) based on https://www.sigmaaldrich.com/US/en/technical-documents/technical-article/genomics/pcr/quantitation-of-oligos.

    :param sequence: primer sequence
    """
    e260 = 0
    for i in range(len(sequence) - 1):
        e260 += NN53_EXTINCTION_COEFFICIENTS[sequence[i]][sequence[i + 1]]
    for base in sequence[1:-1]:
        e260 -= BASE_EXTINCTION_COEFFICIENTS[base]
    return e260
