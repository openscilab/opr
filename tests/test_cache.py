import itertools
from opr.params import VALID_BASES
from opr import Primer, MeltingTemperature

TEST_CASE_NAME = "Cache tests"


def test_mwc():
    oprimer = Primer("ATCGATCGATCGATCGAT")
    molecular_weight_first = oprimer.molecular_weight
    assert oprimer.is_computed("molecular_weight")

    molecular_weight_second = oprimer.molecular_weight
    assert round(molecular_weight_first, 1) == round(molecular_weight_second, 1)


def test_gc_content():
    oprimer = Primer("ATTCG")
    gc_content_first = oprimer.gc_content
    assert oprimer.is_computed("gc_content")

    gc_content_second = oprimer.gc_content
    assert gc_content_first == gc_content_second


def test_gc_clamp():
    oprimer = Primer("ATCGATCGATCGATCGGTCG")
    gc_clamp_first = oprimer.gc_clamp
    assert oprimer.is_computed("gc_clamp")

    gc_clamp_second = oprimer.gc_clamp
    assert gc_clamp_first == gc_clamp_second


def test_melt_temp():
    oprimer = Primer("ATCGATCGATCGATCGATCG")
    basic_melt_temp_first = oprimer.melting_temperature(MeltingTemperature.BASIC)
    assert oprimer.is_computed("melting_temperature")[MeltingTemperature.BASIC]

    basic_melt_temp_second = oprimer.melting_temperature(MeltingTemperature.BASIC)
    assert round(basic_melt_temp_first, 1) == round(basic_melt_temp_second, 1)


def test_single_runs():
    oprimer = Primer("AAAAATTCGGGGATCCCCG")
    runs_first = oprimer.single_runs
    assert oprimer.is_computed("single_runs")

    runs_second = oprimer.single_runs
    assert all(runs_first[pair] == runs_second[pair] for pair in runs_first)
    assert all(runs_second[pair] == runs_first[pair] for pair in runs_second)


def test_double_runs():
    p1 = Primer("ATATCGAACACACACACA")
    double_runs_first = p1.double_runs
    assert p1.is_computed("double_runs")

    double_runs_second = p1.double_runs
    assert len(double_runs_second) == len(double_runs_first)
    assert all(double_runs_first[pair] == double_runs_second[pair] for pair in double_runs_first)
    assert all(double_runs_second[pair] == double_runs_first[pair] for pair in double_runs_second)
