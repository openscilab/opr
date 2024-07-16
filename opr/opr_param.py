# -*- coding: utf-8 -*-
"""Parameters and constants."""
OPR_VERSION = "0.1"
VALID_BASES = set('ATCG')
DNA_COMPLEMENT_MAP = {"A":"T", "C":"G", "G":"C", "T":"A"}

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
PRIMER_SEQUENCE_LENGTH_ERROR = "Primer length should be between 18 and 30 nucleotides."
PRIMER_SEQUENCE_VALID_BASES_ERROR = "Primer sequence should only contain the nucleotide bases A, T, C, and G."
PRIMER_SEQUENCE_VALID_GC_CONTENT_RANGE_ERROR = "Primer GC content should be between 30% and 80%."
PRIMER_READ_ONLY_ATTRIBUTE_ERROR = "This attribute is read-only."
PRIMER_NOT_REMOVABLE_ATTRIBUTE_ERROR = "This attribute is not removable."

