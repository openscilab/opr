# -*- coding: utf-8 -*-
"""OPR parameters and constants."""
OPR_VERSION = "0.3"
VALID_BASES = set('ATCG')
DNA_COMPLEMENT_MAP = {"A": "T", "C": "G", "G": "C", "T": "A"}

PRIMER_LOWER_LENGTH = 18
PRIMER_HIGHEST_LENGTH = 30
PRIMER_LOWEST_GC_RANGE = 0.3
PRIMER_HIGHEST_GC_RANGE = 0.8

A_WEIGHT = 313.21
T_WEIGHT = 304.2
C_WEIGHT = 289.18
G_WEIGHT = 329.21
ANHYDROUS_MOLECULAR_WEIGHT_CONSTANT = 61.96

BASE_EXTINCTION_COEFFICIENTS = {  # L ⋅ mmol-1 ⋅ cm-1
    "A": 15.4,
    "C": 7.4,
    "G": 11.5,
    "T": 8.7
}
NN53_EXTINCTION_COEFFICIENTS = {  # L ⋅ mmol-1 ⋅ cm-1
    "A": {
        "A": 27.4,
        "C": 21.2,
        "G": 25.0,
        "T": 22.8
    },
    "C": {
        "A": 21.2,
        "C": 14.6,
        "G": 18.0,
        "T": 15.2
    },
    "G": {
        "A": 25.2,
        "C": 17.6,
        "G": 21.6,
        "T": 20.0
    },
    "T": {
        "A": 23.4,
        "C": 16.2,
        "G": 19.0,
        "T": 16.8
    }
}

DEFAULT_PRIMER_NAME = "unknown"

PRIMER_SEQUENCE_TYPE_ERROR = "Primer sequence should be a string variable."
PRIMER_SEQUENCE_LENGTH_WARNING = "The recommended range for primer length is between 18 and 30."
PRIMER_SEQUENCE_VALID_BASES_ERROR = "Primer sequence should only contain the nucleotide bases A, T, C, and G."
PRIMER_SEQUENCE_VALID_GC_CONTENT_RANGE_WARNING = "The recommended range for GC content is between 30% and 80%."

PRIMER_ADDITION_ERROR = "You can only add two Primer objects."
PRIMER_MULTIPLICATION_ERROR = "The primer sequence can only be multiplied by an integer."

PRIMER_MELTING_TEMPERATURE_NOT_IMPLEMENTED_ERROR = "This method for calculating melting temperature has not been implemented."

PRIMER_ATTRIBUTE_NOT_COMPUTABLE_ERROR = "This attribute either doesn't exist or cannot be computed/cached."
