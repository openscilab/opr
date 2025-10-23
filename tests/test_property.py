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
    assert oprimer.molecular_formula == "C176H221N70O104P17"


def test_molecular_formula2(): # Reference: https://atdbio.com/tools/oligo-calculator
    oprimer = Primer("AAAAAAAAAAAAAAAAAA")
    assert oprimer.molecular_formula == "C180H217N90O88P17"


def test_molecular_formula3(): # Reference: https://atdbio.com/tools/oligo-calculator
    oprimer = Primer("ATCG")
    assert oprimer.molecular_formula == "C39H50N15O22P3"


def test_molecular_formula4(): # Reference: https://atdbio.com/tools/oligo-calculator
    oprimer = Primer("T")
    assert oprimer.molecular_formula == "C10H14N2O5"
