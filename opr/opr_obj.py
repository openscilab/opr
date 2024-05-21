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
