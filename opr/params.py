# -*- coding: utf-8 -*-
"""OPR parameters and constants."""
OPR_VERSION = "0.1"
VALID_BASES = set('ATCG')
DNA_COMPLEMENT_MAP = {"A": "T", "C": "G", "G": "C", "T": "A"}
CHEMICAL_FORMULA_BASES = {
    'A': {
        'C': 10,
        'H': 13,
        'N': 5,
        'O': 3,
    },
    'T': {
        'C': 10,
        'H': 14,
        'N': 2,
        'O': 5,
    },
    'C': {
        'C': 9,
        'H': 13,
        'N': 3,
        'O': 4,
    },
    'G': {
        'C': 10,
        'H': 13,
        'N': 5,
        'O': 4,
    },
}
CHEMICAL_FORMULA_WATER = {
    'H': 2,
    'O': 1,
}
CHEMICAL_FORMULA_FORMAT = "C{0}H{1}N{2}O{3}P{4}"
CHEMICAL_FORMULA_FORMAT_SHORT = "C{0}H{1}N{2}O{3}"

PRIMER_LOWER_LENGTH = 18
PRIMER_HIGHEST_LENGTH = 30
PRIMER_LOWEST_GC_RANGE = 0.3
PRIMER_HIGHEST_GC_RANGE = 0.8

A_WEIGHT = 313.21
T_WEIGHT = 304.2
C_WEIGHT = 289.18
G_WEIGHT = 329.21
ANHYDROUS_MOLECULAR_WEIGHT_CONSTANT = 61.96

PRIMER_SEQUENCE_TYPE_ERROR = "Primer sequence should be a string variable."
PRIMER_SEQUENCE_LENGTH_WARNING = "The recommended range for primer length is between {0} and {1}.".format(
    PRIMER_LOWER_LENGTH, PRIMER_HIGHEST_LENGTH)
PRIMER_SEQUENCE_VALID_BASES_ERROR = "Primer sequence should only contain the nucleotide bases A, T, C, and G."
PRIMER_SEQUENCE_VALID_GC_CONTENT_RANGE_WARNING = "The recommended range for GC content is between {0}% and {1}%.".format(
    PRIMER_LOWEST_GC_RANGE * 100, PRIMER_HIGHEST_GC_RANGE * 100)
PRIMER_READ_ONLY_ATTRIBUTE_ERROR = "This attribute is read-only."
PRIMER_NOT_REMOVABLE_ATTRIBUTE_ERROR = "This attribute is not removable."

PRIMER_ADDITION_ERROR = "You can only add two Primer objects."
PRIMER_MULTIPLICATION_ERROR = "The primer sequence can only be multiplied by an integer."

PRIMER_MELTING_TEMPERATURE_NOT_IMPLEMENTED_ERROR = "This method for calculating melting temperature has not been implemented."
