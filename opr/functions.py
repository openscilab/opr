# -*- coding: utf-8 -*-
"""OPR functions."""
from .params import A_WEIGHT, T_WEIGHT, C_WEIGHT, G_WEIGHT
from .params import ANHYDROUS_MOLECULAR_WEIGHT_CONSTANT
from .params import CHEMICAL_FORMULA_FORMAT, CHEMICAL_FORMULA_FORMAT_SHORT
from .params import CHEMICAL_FORMULA_BASES, CHEMICAL_FORMULA_WATER, CHEMICAL_FORMULA_PHOSPHODIESTER


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


def gc_content_calc(sequence):
    """
    Calculate gc content.

    :param sequence: primer nucleotides sequence
    :type sequence: str
    :return: gc content as float
    """
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence)


def chemical_formula_calc(sequence):
    """
    Calculate the chemical formula.

    :param sequence: primer nucleotides sequence
    :type sequence: str
    :return: chemical formula as dict
    """
    count_mapping = {
        'A': sequence.count('A'),
        'T': sequence.count('T'),
        'C': sequence.count('C'),
        'G': sequence.count('G'),
    }
    n = len(sequence)

    carbon_count = sum([count_mapping[x] * y['C'] for x, y in CHEMICAL_FORMULA_BASES.items()])
    hydrogen_count = sum([count_mapping[x] * y['H'] for x, y in CHEMICAL_FORMULA_BASES.items()])
    nitrogen_count = sum([count_mapping[x] * y['N'] for x, y in CHEMICAL_FORMULA_BASES.items()])
    oxygen_count = sum([count_mapping[x] * y['O'] for x, y in CHEMICAL_FORMULA_BASES.items()])
    # A water is removed from the formula for each phosphodiester bond
    hydrogen_count -= (n - 1) * CHEMICAL_FORMULA_WATER['H']
    hydrogen_count += (n - 1) * CHEMICAL_FORMULA_PHOSPHODIESTER['H']
    oxygen_count -= (n - 1) * CHEMICAL_FORMULA_WATER['O']
    oxygen_count += (n - 1) * CHEMICAL_FORMULA_PHOSPHODIESTER['O']
    phosphor_count = (n - 1) * CHEMICAL_FORMULA_PHOSPHODIESTER['P']

    if len(sequence) == 1:
        return CHEMICAL_FORMULA_FORMAT_SHORT.format(carbon_count, hydrogen_count, nitrogen_count, oxygen_count)
    return CHEMICAL_FORMULA_FORMAT.format(carbon_count, hydrogen_count, nitrogen_count, oxygen_count, phosphor_count)
