# -*- coding: utf-8 -*-
"""OPR functions."""
import math
from .params import A_WEIGHT, T_WEIGHT, C_WEIGHT, G_WEIGHT, ANHYDROUS_MOLECULAR_WEIGHT_CONSTANT


def molecular_weight_calc(sequence):
    """
    Calculate molecular weight.

    :param sequence: primer nucleotides sequence
    :type sequence: str
    :return: molecular weight as float
    """
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    return (a_count * A_WEIGHT) + (t_count * T_WEIGHT) + (c_count * C_WEIGHT) + \
        (g_count * G_WEIGHT) - ANHYDROUS_MOLECULAR_WEIGHT_CONSTANT


def basic_melting_temperature_calc(sequence):
    """
    Calculate basic melting temperature.

    :param sequence: primer nucleotides sequence
    :type sequence: str
    :return: melting temperature as float
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


def salt_adjusted_melting_temperature_calc(sequence, salt):
    """
    Calculate the salt-adjusted melting temperature (Tm) of a primer sequence.

    :param sequence: Primer nucleotides sequence
    :type sequence: str
    :param salt: Sodium ion concentration in moles (unit mM)
    :type salt: float
    :return: Salt-adjusted melting temperature as float
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


def gc_clamp_calc(sequence):
    """
    Calculate GC clamp.

    :param sequence: primer sequence
    :type sequence: str
    :return: number of guanine (G) or cytosine (C) bases in the last 5 bases of the primer
    """
    if len(sequence) < 5:
        return 0
    return sequence[-5:].count('G') + sequence[-5:].count('C')
