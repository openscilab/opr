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
    return (a_count * A_WEIGHT) + (t_count * T_WEIGHT) + (c_count * C_WEIGHT) + (g_count * G_WEIGHT) - ANHYDROUS_MOLECULAR_WEIGHT_CONSTANT