# -*- coding: utf-8 -*-
"""OPR parameters and constants."""
OPR_VERSION = "0.5"
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

CODONS_TO_AMINO_ACIDS_SHORT = {
    "UUU": "F", "UUC": "F",
    "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S",
    "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y",
    "UAA": "*", "UAG": "*",
    "UGU": "C", "UGC": "C",
    "UGA": "*", "UGG": "W",

    "CUU": "L", "CUC": "L",
    "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P",
    "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R",
    "CGA": "R", "CGG": "R",

    "AUU": "I", "AUC": "I",
    "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T",
    "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S",
    "AGA": "R", "AGG": "R",

    "GUU": "V", "GUC": "V",
    "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A",
    "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G",
    "GGA": "G", "GGG": "G",
}

CODONS_TO_AMINO_ACIDS_LONG = {
    "UUU": "Phe", "UUC": "Phe",
    "UUA": "Leu", "UUG": "Leu",
    "UCU": "Ser", "UCC": "Ser",
    "UCA": "Ser", "UCG": "Ser",
    "UAU": "Tyr", "UAC": "Tyr",
    "UAA": "Stop","UAG": "Stop",
    "UGU": "Cys", "UGC": "Cys",
    "UGA": "Stop","UGG": "Trp",

    "CUU": "Leu", "CUC": "Leu",
    "CUA": "Leu", "CUG": "Leu",
    "CCU": "Pro", "CCC": "Pro",
    "CCA": "Pro", "CCG": "Pro",
    "CAU": "His", "CAC": "His",
    "CAA": "Gln", "CAG": "Gln",
    "CGU": "Arg", "CGC": "Arg",
    "CGA": "Arg", "CGG": "Arg",

    "AUU": "lle", "AUC": "lle", 
    "AUA": "lle", "AUG": "Met",
    "ACU": "Thr", "ACC": "Thr",
    "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn",
    "AAA": "Lys", "AAG": "Lys",
    "AGU": "Ser", "AGC": "Ser",
    "AGA": "Arg", "AGG": "Arg",

    "GUU": "Val", "GUC": "Val",
    "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala",
    "GCA": "Ala", "GCG": "Ala",
    "GAU": "Asp", "GAC": "Asp",
    "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly",
    "GGA": "Gly", "GGG": "Gly",
}

FRAME_ERROR = "Parameter `frame` must be 1, 2, or 3."
DEFAULT_PRIMER_NAME = "unknown"

PRIMER_SEQUENCE_TYPE_ERROR = "Primer sequence should be a string variable."
PRIMER_SEQUENCE_LENGTH_WARNING = "The recommended range for primer length is between 18 and 30."
PRIMER_SEQUENCE_VALID_BASES_ERROR = "Primer sequence should only contain the nucleotide bases A, T, C, and G."
PRIMER_SEQUENCE_VALID_GC_CONTENT_RANGE_WARNING = "The recommended range for GC content is between 30% and 80%."

PRIMER_ADDITION_ERROR = "You can only add two Primer objects."
PRIMER_MULTIPLICATION_ERROR = "The primer sequence can only be multiplied by an integer."

PRIMER_MELTING_TEMPERATURE_NOT_IMPLEMENTED_ERROR = "This method for calculating melting temperature has not been implemented."

PRIMER_ATTRIBUTE_NOT_COMPUTABLE_ERROR = "This attribute either doesn't exist or cannot be computed/cached."

# For DNA
# Nearest-neighbor parameters (ΔH in kcal/mol, ΔS in kcal/K·mol)
# ref: https://www.sigmaaldrich.com/CA/en/technical-documents/protocol/genomics/pcr/oligos-melting-temp
NN_PARAMS = {
    'AA': (-9.1, -0.0240),
    'AT': (-8.6, -0.0239), 
    'TA': (-6.0, -0.0169),
    'CA': (-5.8, -0.0129), 
    'GT': (-6.5, -0.0173),
    'CT': (-7.8, -0.0208),
    'GA': (-5.6, -0.0135),
    'CG': (-11.9, -0.0278), 
    'GC': (-11.1, -0.0267),
    'GG': (-11.0, -0.0266), 
}
