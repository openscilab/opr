# -*- coding: utf-8 -*-
"""OPR primer."""
from enum import Enum
from warnings import warn
from .errors import OPRBaseError
from .params import VALID_BASES
from .params import PRIMER_SEQUENCE_TYPE_ERROR, PRIMER_SEQUENCE_LENGTH_WARNING, PRIMER_SEQUENCE_VALID_BASES_ERROR, PRIMER_SEQUENCE_VALID_GC_CONTENT_RANGE_WARNING
from .params import PRIMER_LOWER_LENGTH, PRIMER_HIGHEST_LENGTH, PRIMER_LOWEST_GC_RANGE, PRIMER_HIGHEST_GC_RANGE
from .params import PRIMER_READ_ONLY_ATTRIBUTE_ERROR, PRIMER_NOT_REMOVABLE_ATTRIBUTE_ERROR
from .params import A_WEIGHT, T_WEIGHT, C_WEIGHT, G_WEIGHT, ANHYDROUS_MOLECULAR_WEIGHT_CONSTANT
from .params import DNA_COMPLEMENT_MAP
from .params import PRIMER_ADDITION_ERROR, PRIMER_MULTIPICATION_ERROR
from .params import PRIMER_SUPPORTED_MELTING_TEMPERATURE_CALCULATIONS


class MeltingTemperatureMode(Enum):
    """Mode used to calculate the Melting Temperature of the Primer accordingly."""

    BASIC = 1
    SALT_ADJUSTED = 2
    NEAREST_NEIGHBOR = 3


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
        self._gc_content = None
        self._melting_temperature = {
            MeltingTemperatureMode.BASIC: None,
            MeltingTemperatureMode.SALT_ADJUSTED: None,
            MeltingTemperatureMode.NEAREST_NEIGHBOR: None,
        }

    def __len__(self):
        """
        Return the length of the Primer sequence.

        :return: length of the Primer sequence
        """
        return len(self._sequence)

    def __add__(self, other_primer):
        """
        Concatenate the sequences of the current Primer with another one.

        :param other_primer: another Primer to concat its sequence to the current Primer
        :type other_primer: Primer
        :return: new Primer with concatenated sequence
        """
        if isinstance(other_primer, Primer):
            return Primer(self._sequence + other_primer._sequence)
        raise OPRBaseError(PRIMER_ADDITION_ERROR)

    def __mul__(self, number):
        """
        Multiply the Primer sequence `number` times.

        :param number: times to concat the Primer sequence to itself
        :type number: int
        :return: new Primer with multiplied sequence
        """
        if isinstance(number, int):
            return Primer(self._sequence * number)
        raise OPRBaseError(PRIMER_MULTIPICATION_ERROR)

    def reverse(self, inplace=False):
        """
        Reverse sequence.

        :param inplace: inplace flag
        :type inplace: bool
        :return: new Primer object or None
        """
        new_seq = self._sequence[::-1]
        if inplace:
            self._sequence = new_seq
        else:
            return Primer(primer_sequence=new_seq)

    def complement(self, inplace=False):
        """
        Complement sequence.

        :param inplace: inplace flag
        :type inplace: bool
        :return: new Primer object or None
        """
        new_seq = ""
        for item in self._sequence:
            new_seq += DNA_COMPLEMENT_MAP[item]
        if inplace:
            self._sequence = new_seq
        else:
            return Primer(primer_sequence=new_seq)

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
            warn(PRIMER_SEQUENCE_LENGTH_WARNING, RuntimeWarning)

        if not all(base in VALID_BASES for base in primer_sequence):
            raise OPRBaseError(PRIMER_SEQUENCE_VALID_BASES_ERROR)
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
        self._molecular_weight = (a_count * A_WEIGHT) + (t_count * T_WEIGHT) + (c_count *
                                                                                C_WEIGHT) + (g_count * G_WEIGHT) - ANHYDROUS_MOLECULAR_WEIGHT_CONSTANT
        return self._molecular_weight

    @molecular_weight.setter
    def molecular_weight(self, _):
        raise OPRBaseError(PRIMER_READ_ONLY_ATTRIBUTE_ERROR)

    @molecular_weight.deleter
    def molecular_weight(self, _):
        self._molecular_weight = None

    @property
    def gc_content(self):
        """
        Calculate gc content.

        :return: gc content
        """
        if self._gc_content is None:
            gc_count = self._sequence.count('G') + self._sequence.count('C')
            self._gc_content = gc_count / len(self._sequence)
        if self._gc_content < PRIMER_LOWEST_GC_RANGE or self._gc_content > PRIMER_HIGHEST_GC_RANGE:
            warn(PRIMER_SEQUENCE_VALID_GC_CONTENT_RANGE_WARNING, RuntimeWarning)
        return self._gc_content

    @gc_content.setter
    def gc_content(self, _):
        raise OPRBaseError(PRIMER_READ_ONLY_ATTRIBUTE_ERROR)

    @gc_content.deleter
    def gc_content(self, _):
        raise OPRBaseError(PRIMER_NOT_REMOVABLE_ATTRIBUTE_ERROR)

    @property
    def melting_temperature(self):
        """
        Calculate(if needed) the melting temperature.

        :return: approximated melting temperatures
        """
        a_count = self._sequence.count('A')
        t_count = self._sequence.count('T')
        c_count = self._sequence.count('C')
        g_count = self._sequence.count('G')
        warn(PRIMER_SUPPORTED_MELTING_TEMPERATURE_CALCULATIONS)
        if self._melting_temperature[MeltingTemperatureMode.BASIC] is None:
            if len(self) <= 13:
                # Tm= (wA+xT) * 2 + (yG+zC) * 4
                # where w,x,y,z are the number of the bases A,T,G,C in the sequence,
                # respectively (from Marmur,J., and Doty,P. (1962) J Mol Biol 5:109-118
                # [PubMed]).
                self._melting_temperature[MeltingTemperatureMode.BASIC] = (
                    a_count + t_count) * 2 + (g_count + c_count) * 4
            else:
                # Tm= 64.9 +41 * (yG+zC-16.4)/(wA+xT+yG+zC)
                # See Wallace,R.B., Shaffer,J., Murphy,R.F., Bonner,J., Hirose,T., and
                # Itakura,K. (1979) Nucleic Acids Res 6:3543-3557 (Abstract) and
                # Sambrook,J., and Russell,D.W. (2001) Molecular Cloning: A Laboratory
                # Manual. Cold Spring Harbor Laboratory Press; Cold Spring Harbor, NY.
                # (CHSL Press)
                self._melting_temperature[MeltingTemperatureMode.BASIC] = 64.9 + 41 * \
                    ((g_count + c_count - 16.4) / (a_count + t_count + g_count + c_count))
        return self._melting_temperature

    @melting_temperature.setter
    def melting_temperature(self, _):
        raise OPRBaseError(PRIMER_READ_ONLY_ATTRIBUTE_ERROR)

    @melting_temperature.deleter
    def melting_temperature(self, _):
        self.melting_temperature = {
            MeltingTemperatureMode.BASIC: None,
            MeltingTemperatureMode.SALT_ADJUSTED: None,
            MeltingTemperatureMode.NEAREST_NEIGHBOR: None,
        }
