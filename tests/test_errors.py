import pytest
from opr import Primer, MeltingTemperature, OPRBaseError

TEST_CASE_NAME = "errors testcase"

def test_addition():
    oprimer = Primer("ATCGATCGATCGATCGAT")
    with pytest.raises(OPRBaseError, match=r"You can only add two Primer objects."):
        oprimer + 2

def test_multiply():
    oprimer = Primer("ATCGATCGATCGATCGAT")
    with pytest.raises(OPRBaseError, match=r"The primer sequence can only be multiplied by an integer."):
        oprimer * "2"

def test_validate_primer():
    with pytest.raises(OPRBaseError, match=r"Primer sequence should only contain the nucleotide bases A, T, C, and G."):
        oprimer = Primer("ATCGATCGATCGATCGAF")

def test_melting_temperature():
    oprimer = Primer("ATCGATCGATCGATCGAT")
    with pytest.raises(NotImplementedError, match=r"This method for calculating melting temperature has not been implemented."):
        oprimer.melting_temperature(method=MeltingTemperature.NEAREST_NEIGHBOR)




