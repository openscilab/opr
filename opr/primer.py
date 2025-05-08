# -*- coding: utf-8 -*-
"""OPR primer."""
import re
import itertools
from enum import Enum
from warnings import warn
from .errors import OPRBaseError
from .params import DEFAULT_PRIMER_NAME, VALID_BASES
from .params import PRIMER_SEQUENCE_TYPE_ERROR, PRIMER_SEQUENCE_LENGTH_WARNING, PRIMER_SEQUENCE_VALID_BASES_ERROR, PRIMER_SEQUENCE_VALID_GC_CONTENT_RANGE_WARNING
from .params import PRIMER_LOWER_LENGTH, PRIMER_HIGHEST_LENGTH, PRIMER_LOWEST_GC_RANGE, PRIMER_HIGHEST_GC_RANGE
from .params import DNA_COMPLEMENT_MAP
from .params import PRIMER_ADDITION_ERROR, PRIMER_MULTIPLICATION_ERROR
from .params import PRIMER_MELTING_TEMPERATURE_NOT_IMPLEMENTED_ERROR
from .params import PRIMER_ATTRIBUTE_NOT_COMPUTABLE_ERROR
from .functions import molecular_weight_calc, basic_melting_temperature_calc, salt_adjusted_melting_temperature_calc, gc_clamp_calc
from .functions import nearest_neighbor_melting_temperature_calc
from .functions import e260_ssnn_calc, protein_seq_calc


class MeltingTemperature(Enum):
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

    def __init__(self, sequence, name=DEFAULT_PRIMER_NAME, salt=50):
        """
        Initialize the Primer instance.

        :param sequence: primer nucleotides sequence
        :type sequence: str
        :param name: primer name
        :type name: str
        :param salt: Sodium ion concentration in millimoles (unit mM)
        :type salt: float
        :return: an instance of the Primer class
        """
        self._sequence = Primer.validate_primer(sequence)
        self._name = name
        self._molecular_weight = None
        self._gc_content = None
        self._gc_clamp = None
        self._single_runs = None
        self._double_runs = None
        self._E260 = None
        self._salt_level = salt
        self._melting_temperature = {
            MeltingTemperature.BASIC: None,
            MeltingTemperature.SALT_ADJUSTED: None,
            MeltingTemperature.NEAREST_NEIGHBOR: None,
        }

        # Track computed attributes
        self._computed = {
            "molecular_weight": False,
            "gc_content": False,
            "gc_clamp": False,
            "single_runs": False,
            "double_runs": False,
            "E260": False,
            "melting_temperature": {
                MeltingTemperature.BASIC: False,
                MeltingTemperature.SALT_ADJUSTED: False,
                MeltingTemperature.NEAREST_NEIGHBOR: False,
            },
        }

    def is_computed(self, attr):
        """
        Check whether the given attribute has been computed.

        :param attr: The attribute to check.
        :type attr: str
        :return: True if the attribute is previously computed, False otherwise.
        """
        if attr not in self._computed:
            raise OPRBaseError(PRIMER_ATTRIBUTE_NOT_COMPUTABLE_ERROR)
        return self._computed.get(attr, False)

    def __len__(self):
        """
        Return the length of the Primer sequence.

        :return: length of the Primer sequence
        """
        return len(self._sequence)

    def __eq__(self, other_primer):
        """
        Check primers equality.

        :param other_primer: another Primer
        :type other_primer: Primer
        :return: result as bool
        """
        if isinstance(other_primer, Primer):
            return self._sequence == other_primer._sequence
        return False

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
        raise OPRBaseError(PRIMER_MULTIPLICATION_ERROR)

    def __contains__(self, sequence):
        """
        Check if the Primer contains the given sequence.

        :param sequence: sequence
        :type sequence: str/Primer
        :return: bool
        """
        if isinstance(sequence, str):
            return sequence in self._sequence
        elif isinstance(sequence, Primer):
            return sequence._sequence in self._sequence
        return False

    def __str__(self):
        """
        Primer object string representation method.

        :return: primer sequence as str
        """
        return self._sequence

    def __iter__(self):
        """
        Iterate through Primer.

        :return: base as Generator[str]
        """
        yield from self.sequence

    def reverse(self, inplace=False):
        """
        Reverse sequence.

        :param inplace: inplace flag
        :type inplace: bool
        :return: new Primer object or None
        """
        new_sequence = self._sequence[::-1]
        if inplace:
            self._sequence = new_sequence
        else:
            return Primer(sequence=new_sequence)

    def complement(self, inplace=False):
        """
        Complement sequence.

        :param inplace: inplace flag
        :type inplace: bool
        :return: new Primer object or None
        """
        new_sequence = ""
        for item in self._sequence:
            new_sequence += DNA_COMPLEMENT_MAP[item]
        if inplace:
            self._sequence = new_sequence
        else:
            return Primer(sequence=new_sequence)

    def to_rna(self):
        """
        Convert DNA sequence to RNA.

        :return: str
        """
        return self._sequence.replace('T', 'U')

    def to_protein(self, frame=1):
        """
        Convert DNA sequence to protein.

        :param frame: reading frame (1, 2, or 3)
        :type frame: int

        :return: str
        """
        return protein_seq_calc(self.to_rna(), frame)

    @staticmethod
    def validate_primer(sequence):
        """
        Validate the given primer sequence.

        :param sequence: primer nucleotides sequence
        :type sequence: any
        :return: an uppercased primer sequence
        """
        if not isinstance(sequence, str):
            raise OPRBaseError(PRIMER_SEQUENCE_TYPE_ERROR)
        sequence = sequence.upper()

        if len(sequence) < PRIMER_LOWER_LENGTH or len(sequence) > PRIMER_HIGHEST_LENGTH:
            warn(PRIMER_SEQUENCE_LENGTH_WARNING, RuntimeWarning)

        if not all(base in VALID_BASES for base in sequence):
            raise OPRBaseError(PRIMER_SEQUENCE_VALID_BASES_ERROR)
        return sequence

    @property
    def sequence(self):
        """
        Return the primer sequence.

        :return: primer sequence
        """
        return self._sequence

    @property
    def name(self):
        """
        Return the primer name.

        :return: primer name
        """
        return self._name

    @property
    def molecular_weight(self):
        """
        Calculate(if needed) the molecular weight.

        :return: molecular weight
        """
        if not self._computed["molecular_weight"]:
            self._molecular_weight = molecular_weight_calc(self._sequence)
            self._computed["molecular_weight"] = True
        return self._molecular_weight

    @property
    def gc_content(self):
        """
        Calculate gc content.

        :return: gc content
        """
        if not self._computed["gc_content"]:
            gc_count = self._sequence.count('G') + self._sequence.count('C')
            self._gc_content = gc_count / len(self._sequence)
            self._computed["gc_content"] = True
        if self._gc_content < PRIMER_LOWEST_GC_RANGE or self._gc_content > PRIMER_HIGHEST_GC_RANGE:
            warn(PRIMER_SEQUENCE_VALID_GC_CONTENT_RANGE_WARNING, RuntimeWarning)
        return self._gc_content

    @property
    def gc_clamp(self):
        """
        Calculate GC clamp of the primer.

        :return: GC clamp of the primer
        """
        if not self._computed["gc_clamp"]:
            self._gc_clamp = gc_clamp_calc(self._sequence)
            self._computed["gc_clamp"] = True
        return self._gc_clamp

    @property
    def single_runs(self):
        """
        Calculate Single Runs of the primer.

        Run length refers to how many times a single base is repeated consecutively in the primer.

        :return: single runs of the primer
        """
        if not self._computed["single_runs"]:
            self._single_runs = {}
            for base in VALID_BASES:
                self._single_runs[base] = self.repeats(base, consecutive=True)
            self._computed["single_runs"] = True
        return self._single_runs

    @property
    def double_runs(self):
        """
        Calculate Double Runs of the primer.

        It refers to how many times each 2-base pairs occurs consecutively in the primer.

        :return: Dictionary of double runs (2-base pairs) and their counts in the primer
        """
        if not self._computed["double_runs"]:
            pairs = [''.join(pair) for pair in itertools.product(VALID_BASES, repeat=2) if pair[0] != pair[1]]
            counts = {}
            for pair in pairs:
                counts[pair] = self.repeats(pair, consecutive=True)
            self._double_runs = counts
            self._computed["double_runs"] = True
        return self._double_runs

    @property
    def E260(self):
        """
        Calculate the extinction coefficient at 260 nm.

        :return: extinction coefficient at 260 nm
        """
        if not self._computed["E260"]:
            self._E260 = e260_ssnn_calc(self._sequence)
            self._computed["E260"] = True
        return self._E260

    def repeats(self, sequence, consecutive=False):
        """
        Count occurrences of a subsequence in a given sequence.

        :param sequence: The sequence to search within.
        :type sequence: str
        :param consecutive: Whether to count only consecutive repeats.
        :type consecutive: bool
        :return: The count of occurrences.
        """
        if consecutive:
            pattern = f"(?:{re.escape(sequence)})+"
            matches = re.findall(f"({pattern})+", self.sequence)
            result = max((len(match) // len(sequence) for match in matches), default=0)
            return result
        else:
            return self.sequence.count(sequence)

    def melting_temperature(self, method=MeltingTemperature.BASIC):
        """
        Calculate(if needed) the melting temperature.

        :param method: requested calculation mode for melting temperature
        :type method: MeltingTemperature
        :return: approximated melting temperature
        """
        if method not in self._computed["melting_temperature"]:
            raise NotImplementedError(PRIMER_MELTING_TEMPERATURE_NOT_IMPLEMENTED_ERROR)
        if self._computed["melting_temperature"][method]:
            return self._melting_temperature[method]

        if method == MeltingTemperature.BASIC:
            self._melting_temperature[MeltingTemperature.BASIC] = basic_melting_temperature_calc(self._sequence)
        elif method == MeltingTemperature.SALT_ADJUSTED:
            self._melting_temperature[MeltingTemperature.SALT_ADJUSTED] = salt_adjusted_melting_temperature_calc(
                self._sequence, self._salt_level)
        elif method == MeltingTemperature.NEAREST_NEIGHBOR:
            self._melting_temperature[MeltingTemperature.NEAREST_NEIGHBOR] = nearest_neighbor_melting_temperature_calc(
                self._sequence, self._salt_level)
        else:
            raise NotImplementedError(PRIMER_MELTING_TEMPERATURE_NOT_IMPLEMENTED_ERROR)
        self._computed["melting_temperature"][method] = True
        return self._melting_temperature[method]
