# -*- coding: utf-8 -*-
"""OPR modules."""
from .opr_error import OPRBaseError
from .opr_param import VALID_BASES
from .opr_param import PRIMER_SEQUENCE_TYPE_ERROR, PRIMER_SEQUENCE_LENGTH_ERROR, PRIMER_SEQUENCE_VALID_BASES_ERROR, PRIMER_SEQUENCE_VALID_GC_CONTENT_RANGE_ERROR
from .opr_param import PRIMER_LOWER_LENGTH, PRIMER_HIGHEST_LENGTH, PRIMER_LOWEST_GC_RANGE, PRIMER_HIGHEST_GC_RANGE
from .opr_param import PRIMER_READ_ONLY_ATTRIBUTE_ERROR, PRIMER_NOT_REMOVABLE_ATTRIBUTE_ERROR
from .opr_param import A_WEIGHT, T_WEIGHT, C_WEIGHT, G_WEIGHT, ANHYDROUS_MOLECULAR_WEIGHT_CONSTANT

class Primer:
    """
    The Primer class facilitates working with the primer sequence.

    >>> oprimer = Primer("ATCGATCGATCGATCGAT")
    >>> oprimer.molecular_weight
    """

    def __init__(self, primer_sequence):
        """
        Initialize the Primer instance.

        :param primer_sequence: primer nucleotides sequence
        :type primer_sequence: str
        :return: an instance of the Primer class
        """
        self._sequence = Primer.validate_primer(primer_sequence)
        self._molecular_weight = None

    @staticmethod
    def validate_primer(primer_sequence):
        """
        Validate the given primer sequence.

        :param primer_sequence: primer nucleotides sequence
        :type primer_sequence: any
        :return: an uppercased primer sequence
        """
        if not isinstance(primer_sequence, str):
            raise OPRBaseError(PRIMER_SEQUENCE_TYPE_ERROR)
        primer_sequence = primer_sequence.upper()

        if len(primer_sequence) < PRIMER_LOWER_LENGTH or len(primer_sequence) > PRIMER_HIGHEST_LENGTH:
            raise OPRBaseError(PRIMER_SEQUENCE_LENGTH_ERROR)

        if not all(base in VALID_BASES for base in primer_sequence):
            raise OPRBaseError(PRIMER_SEQUENCE_VALID_BASES_ERROR)

        gc_count = primer_sequence.count('G') + primer_sequence.count('C')
        gc_content = gc_count / len(primer_sequence)

        if gc_content < PRIMER_LOWEST_GC_RANGE or gc_content > PRIMER_HIGHEST_GC_RANGE:
            raise OPRBaseError(PRIMER_SEQUENCE_VALID_GC_CONTENT_RANGE_ERROR)
        return primer_sequence

    @property
    def sequence(self):
        """
        Return the primer sequence.

        :return: primer sequence
        """
        return self._sequence

    @sequence.setter
    def sequence(self, _):
        raise OPRBaseError(PRIMER_READ_ONLY_ATTRIBUTE_ERROR)

    @sequence.deleter
    def sequence(self, _):
        raise OPRBaseError(PRIMER_NOT_REMOVABLE_ATTRIBUTE_ERROR)

    @property
    def molecular_weight(self):
        """
        Calculate(if needed) the molecular weight.

        :return: molecular weight
        """
        if self._molecular_weight:
            return self._molecular_weight
        a_count = self._sequence.count('A')
        t_count = self._sequence.count('T')
        c_count = self._sequence.count('C')
        g_count = self._sequence.count('G')
        # Anhydrous Molecular Weight = (An x 313.21) + (Tn x 304.2) + (Cn x 289.18) + (Gn x 329.21) - 61.96
        self._molecular_weight = (a_count * A_WEIGHT) + (t_count * T_WEIGHT) + (c_count * C_WEIGHT) + (g_count * G_WEIGHT) - ANHYDROUS_MOLECULAR_WEIGHT_CONSTANT
        return self._molecular_weight    

    @molecular_weight.setter
    def molecular_weight(self, _):
        raise OPRBaseError(PRIMER_READ_ONLY_ATTRIBUTE_ERROR)

    @molecular_weight.deleter
    def molecular_weight(self, _):
        self._molecular_weight = None
