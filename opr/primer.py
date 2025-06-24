# -*- coding: utf-8 -*-
"""OPR primer."""
from __future__ import annotations
from typing import Union, Generator, Optional
from typing import Dict
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
from .params import FRAME_ERROR
from .params import CODONS_TO_AMINO_ACIDS_LONG, CODONS_TO_AMINO_ACIDS_SHORT
from .functions import molecular_weight_calc, basic_melting_temperature_calc, salt_adjusted_melting_temperature_calc, gc_clamp_calc
from .functions import nearest_neighbor_melting_temperature_calc, calculate_thermodynamics_constants
from .functions import e260_ssnn_calc


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

    def __init__(self, sequence: str, name: str = DEFAULT_PRIMER_NAME, salt: float = 50) -> None:
        """
        Initialize the Primer instance.

        :param sequence: primer nucleotides sequence
        :param name: primer name
        :param salt: Sodium ion concentration in millimoles (unit mM)
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
        self._delta_s = None
        self._delta_h = None
        self._protein_seq = {"AA1": {}, "AA3": {}}

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
            "delta_s": False,
            "delta_h": False,
        }

    def is_computed(self, attr: str) -> bool:
        """
        Check whether the given attribute has been computed. Return true if it has been previously computed.

        :param attr: The attribute to check.
        """
        if attr not in self._computed:
            raise OPRBaseError(PRIMER_ATTRIBUTE_NOT_COMPUTABLE_ERROR)
        return self._computed.get(attr, False)

    def __len__(self):
        """Return the length of the Primer sequence."""
        return len(self._sequence)

    def __eq__(self, other_primer: Primer) -> bool:
        """
        Check primers equality. Return true if the sequences are equal.

        :param other_primer: another Primer
        """
        if isinstance(other_primer, Primer):
            return self._sequence == other_primer._sequence
        return False

    def __add__(self, other_primer: Primer) -> Primer:
        """
        Concatenate the sequences of the current Primer with another one and return a new Primer.

        :param other_primer: another Primer to concat its sequence to the current Primer
        """
        if isinstance(other_primer, Primer):
            return Primer(self._sequence + other_primer._sequence)
        raise OPRBaseError(PRIMER_ADDITION_ERROR)

    def __mul__(self, number: int) -> Primer:
        """
        Multiply the Primer sequence `number` times and return a new Primer.

        :param number: times to concat the Primer sequence to itself
        """
        if isinstance(number, int):
            return Primer(self._sequence * number)
        raise OPRBaseError(PRIMER_MULTIPLICATION_ERROR)

    def __contains__(self, sequence: Union[str, Primer]) -> bool:
        """
        Check if the Primer contains the given sequence.

        :param sequence: sequence
        """
        if isinstance(sequence, str):
            return sequence in self._sequence
        elif isinstance(sequence, Primer):
            return sequence._sequence in self._sequence
        return False

    def __str__(self) -> str:
        """Primer object string representation method. Returns the primer sequence."""
        return self._sequence

    def __iter__(self) -> Generator[str, None, None]:
        """Iterate through the primer sequence."""
        yield from self.sequence

    def reverse(self, inplace: bool = False) -> Optional[Primer]:
        """
        Reverse the sequence.

        :param inplace: inplace flag
        """
        new_sequence = self._sequence[::-1]
        if inplace:
            self._sequence = new_sequence
        else:
            return Primer(sequence=new_sequence)

    def complement(self, inplace: bool = False) -> Optional[Primer]:
        """
        Complement sequence.

        :param inplace: inplace flag
        """
        new_sequence = ""
        for item in self._sequence:
            new_sequence += DNA_COMPLEMENT_MAP[item]
        if inplace:
            self._sequence = new_sequence
        else:
            return Primer(sequence=new_sequence)

    def to_rna(self) -> str:
        """Convert DNA sequence to RNA."""
        return self._sequence.replace('T', 'U')

    def to_protein(self, frame: int = 1, multi_letter: bool = False) -> str:
        """
        Convert DNA sequence to protein and return the protein sequence.

        :param frame: reading frame (1, 2, or 3)
        :param multi_letter: whether to return amino acids in 1-letter codes (False) or 3-letter codes (True)
        """
        if frame not in [1, 2, 3]:
            raise OPRBaseError(FRAME_ERROR)

        key = "AA3" if multi_letter else "AA1"
        if frame in self._protein_seq[key]:
            return self._protein_seq[key][frame]

        rna_sequence = self.to_rna()
        start = frame - 1
        protein_aa1 = []
        protein_aa3 = []
        for i in range(start, len(rna_sequence) - 2, 3):
            codon = rna_sequence[i:i+3]
            protein_aa1.append(CODONS_TO_AMINO_ACIDS_SHORT[codon])
            protein_aa3.append(CODONS_TO_AMINO_ACIDS_LONG[codon])

        self._protein_seq["AA1"][frame] = ''.join(protein_aa1)
        self._protein_seq["AA3"][frame] = '-'.join(protein_aa3)

        result = self._protein_seq["AA3"][frame] if multi_letter else self._protein_seq["AA1"][frame]
        return result

    @staticmethod
    def validate_primer(sequence: str) -> str:
        """
        Validate the given primer sequence and return it in uppercase.

        :param sequence: primer nucleotides sequence
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
    def sequence(self) -> str:
        """Return the primer sequence."""
        return self._sequence

    @property
    def name(self) -> str:
        """Return the primer name."""
        return self._name

    @property
    def molecular_weight(self) -> float:
        """Calculate the molecular weight and return it."""
        if not self._computed["molecular_weight"]:
            self._molecular_weight = molecular_weight_calc(self._sequence)
            self._computed["molecular_weight"] = True
        return self._molecular_weight

    @property
    def gc_content(self) -> float:
        """Calculate gc content and return it."""
        if not self._computed["gc_content"]:
            gc_count = self._sequence.count('G') + self._sequence.count('C')
            self._gc_content = gc_count / len(self._sequence)
            self._computed["gc_content"] = True
        if self._gc_content < PRIMER_LOWEST_GC_RANGE or self._gc_content > PRIMER_HIGHEST_GC_RANGE:
            warn(PRIMER_SEQUENCE_VALID_GC_CONTENT_RANGE_WARNING, RuntimeWarning)
        return self._gc_content

    @property
    def gc_clamp(self) -> int:
        """Calculate GC clamp of the primer and return it."""
        if not self._computed["gc_clamp"]:
            self._gc_clamp = gc_clamp_calc(self._sequence)
            self._computed["gc_clamp"] = True
        return self._gc_clamp

    @property
    def single_runs(self) -> Dict[str, int]:
        """
        Calculate Single Runs of the primer and return them.

        Run length refers to how many times a single base is repeated consecutively in the primer.
        """
        if not self._computed["single_runs"]:
            self._single_runs = {}
            for base in VALID_BASES:
                self._single_runs[base] = self.repeats(base, consecutive=True)
            self._computed["single_runs"] = True
        return self._single_runs

    @property
    def double_runs(self) -> Dict[str, int]:
        """
        Calculate Double Run Counts (2-base pairs) of the primer and return them.

        It refers to how many times each 2-base pairs occurs consecutively in the primer.
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
    def E260(self) -> float:
        """Calculate the extinction coefficient at 260 nm and return it."""
        if not self._computed["E260"]:
            self._E260 = e260_ssnn_calc(self._sequence)
            self._computed["E260"] = True
        return self._E260

    @property
    def delta_s(self) -> float:
        """Calculate entropy change, ΔS (in kcal/mol·K), and return it."""
        if not self._computed["delta_s"]:
            self._delta_h , self._delta_s = calculate_thermodynamics_constants(self._sequence)
            self._computed["delta_s"] = True
            self._computed["delta_h"] = True
        return self._delta_s

    @property
    def delta_h(self) -> float:
        """Calculate enthalpy change, ΔH (in kcal/mol), and return it."""
        if not self._computed["delta_h"]:
            self._delta_h , self._delta_s = calculate_thermodynamics_constants(self._sequence)
            self._computed["delta_s"] = True
            self._computed["delta_h"] = True
        return self._delta_h

    def repeats(self, sequence: str, consecutive: bool = False) -> int:
        """
        Count occurrences of a subsequence in a given sequence and return it.

        :param sequence: The sequence to search within.
        :param consecutive: Whether to count only consecutive repeats.
        """
        if consecutive:
            pattern = f"(?:{re.escape(sequence)})+"
            matches = re.findall(f"({pattern})+", self.sequence)
            result = max((len(match) // len(sequence) for match in matches), default=0)
            return result
        else:
            return self.sequence.count(sequence)

    def melting_temperature(self, method: MeltingTemperature = MeltingTemperature.BASIC) -> float:
        """
        Calculate the approximated melting temperature and return it.

        :param method: requested calculation mode for melting temperature
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
        else:
            # the method is MeltingTemperature.NEAREST_NEIGHBOR
            self._melting_temperature[MeltingTemperature.NEAREST_NEIGHBOR] = nearest_neighbor_melting_temperature_calc(
                self._sequence,
                self._salt_level,
                (self.delta_h, self.delta_s)
                )
        self._computed["melting_temperature"][method] = True
        return self._melting_temperature[method]
