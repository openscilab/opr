from opr.opr_obj import Primer

TEST_CASE_NAME = "molecular weight calculation testcase"

def test_mwc():
    oprimer = Primer("ATCGATCGATCGATCGAT")
    assert round(oprimer.molecular_weight, 1) == 5498.7
