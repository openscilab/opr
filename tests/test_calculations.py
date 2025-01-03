from opr import Primer
from opr.primer import MeltingTemperature

TEST_CASE_NAME = "Calculations tests"

def test_mwc_1():
    oprimer = Primer("ATCGATCGATCGATCGAT")
    assert round(oprimer.molecular_weight, 1) == 5498.7

def test_mwc_2():
    oprimer = Primer("ATCGATCGATCGATCGAT")
    molecular_weight = oprimer.molecular_weight
    assert round(oprimer.molecular_weight, 1) == round(molecular_weight, 1)

def test_gc_content_1(): #Reference: https://jamiemcgowan.ie/bioinf/gc_content.html
    oprimer = Primer("ATCG")
    assert oprimer.gc_content == 0.5

def test_gc_content_2(): #Reference: https://jamiemcgowan.ie/bioinf/gc_content.html
    oprimer = Primer("ATTCG")
    assert oprimer.gc_content == 0.4

def test_gc_content_3(): #Reference: https://jamiemcgowan.ie/bioinf/gc_content.html
    oprimer = Primer("ATTTTTT")
    assert oprimer.gc_content == 0

def test_gc_content_4(): #Reference: https://jamiemcgowan.ie/bioinf/gc_content.html
    oprimer = Primer("ATTCG")
    gc_content = oprimer.gc_content
    assert oprimer.gc_content == gc_content

def test_gc_clamp_1(): #Reference: https://www.bioinformatics.org/sms2/pcr_primer_stats.html
    oprimer = Primer("ATCGATCGATCGATCGGTCG")
    assert oprimer.gc_clamp == 4

def test_gc_clamp_2(): #Reference: https://www.bioinformatics.org/sms2/pcr_primer_stats.html
    oprimer = Primer("ATCG")
    assert oprimer.gc_clamp == 0

def test_gc_clamp_3(): #Reference: https://www.bioinformatics.org/sms2/pcr_primer_stats.html
    oprimer = Primer("ACTTA")
    assert oprimer.gc_clamp == 1

def test_gc_clamp_4(): #Reference: https://www.bioinformatics.org/sms2/pcr_primer_stats.html
    oprimer = Primer("ATCGATCGATCGATCGGTCG")
    gc_clamp = oprimer.gc_clamp
    assert oprimer.gc_clamp == gc_clamp

def test_melt_temp_1(): #Reference: http://biotools.nubic.northwestern.edu/OligoCalc.html
    oprimer = Primer("ATCGATCGATCGATCGATCG")
    basic_melt_temp = oprimer.melting_temperature(MeltingTemperature.BASIC)
    assert round(basic_melt_temp,1) == 51.8

def test_melt_temp_2(): #Reference: http://biotools.nubic.northwestern.edu/OligoCalc.html
    oprimer = Primer("ATCG")
    basic_melt_temp = oprimer.melting_temperature(method=MeltingTemperature.BASIC)
    assert round(basic_melt_temp,1) == 12

def test_melt_temp_3(): #Reference: http://biotools.nubic.northwestern.edu/OligoCalc.html
    oprimer = Primer("ATCGATCGATCGATCGATCG")
    basic_melt_temp = oprimer.melting_temperature(MeltingTemperature.BASIC)
    assert round(oprimer.melting_temperature(MeltingTemperature.BASIC),1) == round(basic_melt_temp,1)

def test_single_runs_1(): #Reference: https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("ATCGATCG")
    runs = oprimer.single_runs
    assert runs['A'] == 0 and runs['T'] == 0 and runs['C'] == 0 and runs['G'] == 0

def test_single_runs_2(): #Reference: https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("ATTCGATCCCCG")
    runs = oprimer.single_runs
    assert runs['A'] == 0 and runs['T'] == 2 and runs['C'] == 4 and runs['G'] == 0

def test_single_runs_3(): #Reference: https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("AAAAATTCGGGGATCCCCG")
    runs = oprimer.single_runs
    assert runs['A'] == 5 and runs['T'] == 2 and runs['C'] == 4 and runs['G'] == 4

def test_single_runs_4(): #Reference: https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("AAAAATTCGGGGATCCCCG")
    runs = oprimer.single_runs
    assert oprimer.single_runs['A'] == runs['A'] and oprimer.single_runs['T'] == runs['T'] and oprimer.single_runs['C'] == runs['C'] and oprimer.single_runs['G'] == runs['G']
