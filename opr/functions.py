# -*- coding: utf-8 -*-
"""OPR functions."""
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


def single_run_length(sequence, base):
    """
    Calculate the maximum consecutive occurrence of a Nucleic acid (base) in a sequence.

    :param sequence: primer sequence
    :type sequence: str
    :param base: target Nucleic acid
    :type base: str
    :return: the length of maximum consecutive occurrence
    """
    max_length = 0
    current_length = 0
    for c in sequence:
        if c == base:
            current_length += 1
            max_length = max(max_length, current_length)
        else:
            current_length = 0
    if max_length == 1:
        return 0
    return max_length
