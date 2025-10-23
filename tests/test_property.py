from opr import Primer

TEST_CASE_NAME = "Property tests"


def test_sequence():
    oprimer = Primer("ATCGATCGATCGATCGAT")
    assert oprimer.sequence == "ATCGATCGATCGATCGAT"


def test_name():
    oprimer = Primer("ATCGATCGATCGATCGAT", "primer1")
    assert oprimer.name == "primer1"


def test_molecular_formula1(): # Reference: https://atdbio.com/tools/oligo-calculator
    oprimer = Primer("ATCGGCTAAATCGGCTAA")
    assert oprimer.molecular_formula == "C176 H221 N70 O104 P17"


def test_molecular_formula2(): # Reference: https://atdbio.com/tools/oligo-calculator
    oprimer = Primer("AAAAAAAAAAAAAAAAAA")
    assert oprimer.molecular_formula == "C180 H217 N90 O88 P17"


def test_molecular_formula3(): # Reference: https://atdbio.com/tools/oligo-calculator
    oprimer = Primer("ATCG")
    assert oprimer.molecular_formula == "C39 H50 N15 O22 P3"


def test_molecular_formula4(): # Reference: https://atdbio.com/tools/oligo-calculator
    oprimer = Primer("T")
    assert oprimer.molecular_formula == "C10 H14 N2 O5"
