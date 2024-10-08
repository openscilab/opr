from opr import Primer
from opr.primer import MeltingTemperature

TEST_CASE_NAME = "Calculations tests"

def test_mwc():
    oprimer = Primer("ATCGATCGATCGATCGAT")
    assert round(oprimer.molecular_weight, 1) == 5498.7

def test_gc_content_1(): #Reference: https://jamiemcgowan.ie/bioinf/gc_content.html
    oprimer = Primer("ATCG")
    assert oprimer.gc_content == 0.5

def test_gc_content_2(): #Reference: https://jamiemcgowan.ie/bioinf/gc_content.html
    oprimer = Primer("ATTCG")
    assert oprimer.gc_content == 0.4

def test_gc_content_3(): #Reference: https://jamiemcgowan.ie/bioinf/gc_content.html
    oprimer = Primer("ATTTTTT")
    assert oprimer.gc_content == 0

def test_melt_temp_1(): #Reference: http://biotools.nubic.northwestern.edu/OligoCalc.html
    oprimer = Primer("ATCGATCGATCGATCGATCG")
    basic_melt_temp = oprimer.melting_temperature(MeltingTemperature.BASIC)
    assert round(basic_melt_temp,1) == 51.8

def test_melt_temp_2(): #Reference: http://biotools.nubic.northwestern.edu/OligoCalc.html
    oprimer = Primer("ATCG")
    basic_melt_temp = oprimer.melting_temperature(method=MeltingTemperature.BASIC)
    assert round(basic_melt_temp,1) == 12
