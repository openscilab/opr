import pytest
from opr import Primer, OPRBaseError

TEST_CASE_NAME = "errors testcase"

def test_add_error():
    oprimer = Primer("ATCGATCGATCGATCGAT")
    with pytest.raises(OPRBaseError, match=r"You can only add two Primer objects."):
        oprimer + 2

