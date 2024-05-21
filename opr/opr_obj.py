# -*- coding: utf-8 -*-
"""OPR modules."""
from .opr_error import OPRBaseError
from .opr_param import VALID_BASES


class OPR:
    """
    TODO

    >>> TODO
    """

    def __init__(self, primer_sequence):
        """
        Initialize the Optimized Primer(OPR) instance.

        :param primer_sequence: primer nucleotides sequence
        :type primer_sequence: str
        :return: an instance of the OPR class
        """
        self._sequence = OPR.validate_primer(primer_sequence)
        self._molecular_weight = None

    @staticmethod
    def validate_primer(primer_sequence):
        if not isinstance(primer_sequence, str):
            raise OPRBaseError("Primer sequence should be a string variable.")
        primer_sequence = primer_sequence.upper()

        if len(primer_sequence) < 18 or len(primer_sequence) > 22:
            raise OPRBaseError("Primer length should be between 18 and 22 nucleotides.")

        if not all(base in VALID_BASES for base in primer_sequence):
            raise OPRBaseError("Primer sequence should only contain the nucleotide bases A, T, C, and G.")

        gc_count = primer_sequence.count('G') + primer_sequence.count('C')
        gc_content = gc_count / len(primer_sequence)

        if gc_content < 0.4 or gc_content > 0.6:
            raise OPRBaseError("Primer GC content should be between 40% and 60%.")
        return primer_sequence
    
    @property
    def sequence(self):
        return self._sequence
    
    @sequence.setter
    def sequence(self, any):
        raise OPRBaseError("sequence attribute is read-only.")

    @sequence.deleter
    def sequence(self, any):
        raise OPRBaseError("This attribute is not removable.")

    @property
    def molecular_weight(self):
        if self._molecular_weight:
            return self._molecular_weight
        a_count = self._sequence.count('A')
        t_count = self._sequence.count('T')
        c_count = self._sequence.count('C')
        g_count = self._sequence.count('G')
        # Anhydrous Molecular Weight = (An x 313.21) + (Tn x 304.2) + (Cn x 289.18) + (Gn x 329.21) - 61.96
        self._molecular_weight = (a_count * 313.21) + (t_count * 304.2) + (c_count * 289.18) + (g_count * 329.21) - 61.96
        return self._molecular_weight    

    @molecular_weight.setter
    def molecular_weight(self, any):
        raise OPRBaseError("molecular_weight attribute is read-only.")

    @molecular_weight.deleter
    def molecular_weight(self, any):
        self._molecular_weight = None
