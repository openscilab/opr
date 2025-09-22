from opr import Primer

TEST_CASE_NAME = "Property tests"


def test_sequence():
    oprimer = Primer("ATCGATCGATCGATCGAT")
    assert oprimer.sequence == "ATCGATCGATCGATCGAT"


def test_name():
    oprimer = Primer("ATCGATCGATCGATCGAT", "primer1")
    assert oprimer.name == "primer1"


def test_molecular_formula(): #Reference: https://atdbio.com/tools/oligo-calculator
    oprimer = Primer("ATCGGCTAAATCGGCTAA")
    assert oprimer.molecular_formula == "C176H221N70O104P17"
