from opr import Primer

TEST_CASE_NAME = "Operations tests"


def test_reverse_1():  # Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGATCGATCGATCGAT")
    oprimer_reversed = oprimer.reverse(inplace=False)
    assert oprimer_reversed.sequence == "TAGCTAGCTAGCTAGCTA"


def test_reverse_2():  # Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGATCGATCGATCGAT")
    oprimer_reversed = oprimer.reverse()
    assert oprimer_reversed.sequence == "TAGCTAGCTAGCTAGCTA"


def test_reverse_3():  # Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGATCGATCGATCGAT")
    oprimer.reverse(inplace=True)
    assert oprimer.sequence == "TAGCTAGCTAGCTAGCTA"


def test_complement_1():  # Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGGCTAAATCGGCTAA")
    oprimer_complemented = oprimer.complement(inplace=False)
    assert oprimer_complemented.sequence == "TAGCCGATTTAGCCGATT"


def test_complement_2():  # Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGGCTAAATCGGCTAA")
    oprimer_complemented = oprimer.complement()
    assert oprimer_complemented.sequence == "TAGCCGATTTAGCCGATT"


def test_complement_3():  # Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGGCTAAATCGGCTAA")
    oprimer.complement(inplace=True)
    assert oprimer.sequence == "TAGCCGATTTAGCCGATT"


def test_to_rna_1():  # Reference: https://biomodel.uah.es/en/lab/cybertory/analysis/trans.htm
    oprimer = Primer("ATCGGCTAAATCGGCTAA")
    oprimer_rna = oprimer.to_rna()
    assert oprimer_rna == "AUCGGCUAAAUCGGCUAA"


def test_to_rna_2():  # Reference: https://biomodel.uah.es/en/lab/cybertory/analysis/trans.htm
    oprimer = Primer("ATCGATCGATCG")
    oprimer_rna = oprimer.to_rna()
    assert oprimer_rna == "AUCGAUCGAUCG"


def test_to_protein():  # Reference: https://biomodel.uah.es/en/lab/cybertory/analysis/trans.htm
    oprimer = Primer("ATCGATCG")
    oprimer_protein = oprimer.to_protein()
    assert oprimer_protein == "ID"


def test_to_protein_2():  # Reference: https://en.vectorbuilder.com/tool/dna-translation.html
    oprimer = Primer("ATCGATCG")
    oprimer_protein = oprimer.to_protein(frame=2)
    assert oprimer_protein == "SI"


def test_to_protein_3():  # Reference: https://en.vectorbuilder.com/tool/dna-translation.html
    oprimer = Primer("ATCGATCG")
    oprimer_protein = oprimer.to_protein(frame=3)
    assert oprimer_protein == "RS"


def test_to_protein_4():  # Reference: https://en.vectorbuilder.com/tool/dna-translation.html
    oprimer = Primer("ATCGATCGATCGTAA")
    oprimer_protein = oprimer.to_protein()
    assert oprimer_protein == "IDRS*"


def test_length():
    oprimer = Primer("ATCGGCTAAATCGGCTAA")
    assert len(oprimer) == 18


def test_addition():
    oprimer_1 = Primer("ATCG")
    oprimer_2 = Primer("GATC")
    oprimer_concat = oprimer_1 + oprimer_2
    assert oprimer_concat.sequence == "ATCGGATC"


def test_multiply():
    oprimer_1 = Primer("ATCG")
    oprimer_concat = oprimer_1 * 4
    assert oprimer_concat.sequence == "ATCGATCGATCGATCG"


def test_equality1():
    oprimer_1 = Primer("ATCG")
    oprimer_2 = Primer("ATCG")
    assert oprimer_1 == oprimer_2


def test_equality2():
    oprimer_1 = Primer("ATCG")
    oprimer_2 = Primer("ATCGC")
    assert oprimer_1 != oprimer_2


def test_equality3():
    oprimer = Primer("ATCG")
    assert oprimer != 2


def test_str():
    oprimer_1 = Primer("ATCG")
    oprimer_2 = Primer("ATCGC")
    oprimer_concat = oprimer_1 + oprimer_2
    assert str(oprimer_1) + str(oprimer_2) == oprimer_concat.sequence


def test_iter():
    oprimer_1 = Primer("ATCG")
    sequence = ""
    for base in oprimer_1:
        sequence += base
    assert oprimer_1.sequence == sequence


def test_contains1():
    oprimer_1 = Primer("ATCG")
    oprimer_2 = Primer("ATCGAT")
    assert oprimer_1 in oprimer_2


def test_contains2():
    sequence = "ATCG"
    oprimer = Primer("ATCGAT")
    assert sequence in oprimer


def test_contains3():
    sequence = "ATCG"
    oprimer = Primer("TCGAT")
    assert sequence not in oprimer


def test_contains4():
    oprimer = Primer("TCGAT")
    assert 2 not in oprimer
