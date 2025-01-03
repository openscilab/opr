from opr import Primer, MeltingTemperature

TEST_CASE_NAME = "Cache tests"


def test_mwc():
    oprimer = Primer("ATCGATCGATCGATCGAT")
    molecular_weight = oprimer.molecular_weight
    assert round(oprimer.molecular_weight, 1) == round(molecular_weight, 1)


def test_gc_content():
    oprimer = Primer("ATTCG")
    gc_content = oprimer.gc_content
    assert oprimer.gc_content == gc_content


def test_gc_clamp():
    oprimer = Primer("ATCGATCGATCGATCGGTCG")
    gc_clamp = oprimer.gc_clamp
    assert oprimer.gc_clamp == gc_clamp


def test_melt_temp():
    oprimer = Primer("ATCGATCGATCGATCGATCG")
    basic_melt_temp = oprimer.melting_temperature(MeltingTemperature.BASIC)
    assert round(oprimer.melting_temperature(MeltingTemperature.BASIC), 1) == round(basic_melt_temp, 1)


def test_single_runs():
    oprimer = Primer("AAAAATTCGGGGATCCCCG")
    runs = oprimer.single_runs
    assert oprimer.single_runs['A'] == runs['A'] and oprimer.single_runs['T'] == runs[
        'T'] and oprimer.single_runs['C'] == runs['C'] and oprimer.single_runs['G'] == runs['G']
