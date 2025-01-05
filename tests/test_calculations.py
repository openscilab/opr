from opr import Primer, MeltingTemperature

TEST_CASE_NAME = "Calculations tests"


def test_mwc():
    oprimer = Primer("ATCGATCGATCGATCGAT")
    assert round(oprimer.molecular_weight, 1) == 5498.7


def test_gc_content_1():  # Reference: https://jamiemcgowan.ie/bioinf/gc_content.html
    oprimer = Primer("ATCG")
    assert oprimer.gc_content == 0.5


def test_gc_content_2():  # Reference: https://jamiemcgowan.ie/bioinf/gc_content.html
    oprimer = Primer("ATTCG")
    assert oprimer.gc_content == 0.4


def test_gc_content_3():  # Reference: https://jamiemcgowan.ie/bioinf/gc_content.html
    oprimer = Primer("ATTTTTT")
    assert oprimer.gc_content == 0


def test_gc_clamp_1():  # Reference: https://www.bioinformatics.org/sms2/pcr_primer_stats.html
    oprimer = Primer("ATCGATCGATCGATCGGTCG")
    assert oprimer.gc_clamp == 4


def test_gc_clamp_2():  # Reference: https://www.bioinformatics.org/sms2/pcr_primer_stats.html
    oprimer = Primer("ATCG")
    assert oprimer.gc_clamp == 0


def test_gc_clamp_3():  # Reference: https://www.bioinformatics.org/sms2/pcr_primer_stats.html
    oprimer = Primer("ACTTA")
    assert oprimer.gc_clamp == 1


def test_melt_temp_1():  # Reference: http://biotools.nubic.northwestern.edu/OligoCalc.html
    oprimer = Primer("ATCGATCGATCGATCGATCG")
    basic_melt_temp = oprimer.melting_temperature(MeltingTemperature.BASIC)
    assert round(basic_melt_temp, 1) == 51.8


def test_melt_temp_2():  # Reference: http://biotools.nubic.northwestern.edu/OligoCalc.html
    oprimer = Primer("ATCG")
    basic_melt_temp = oprimer.melting_temperature(method=MeltingTemperature.BASIC)
    assert round(basic_melt_temp, 1) == 12


def test_single_runs_1():  # Reference: https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("ATCGATCG")
    runs = oprimer.single_runs
    assert runs['A'] == 0 and runs['T'] == 0 and runs['C'] == 0 and runs['G'] == 0


def test_single_runs_2():  # Reference: https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("ATTCGATCCCCG")
    runs = oprimer.single_runs
    assert runs['A'] == 0 and runs['T'] == 2 and runs['C'] == 4 and runs['G'] == 0


def test_single_runs_3():  # Reference: https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("AAAAATTCGGGGATCCCCG")
    runs = oprimer.single_runs
    assert runs['A'] == 5 and runs['T'] == 2 and runs['C'] == 4 and runs['G'] == 4
