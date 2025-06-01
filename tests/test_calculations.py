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


def test_melt_temp_3():  # Reference: http://biotools.nubic.northwestern.edu/OligoCalc.html
    oprimer = Primer("CTGGAGGACGGAAGAGGAAGTAA")
    salt_adjusted_melt_temp = oprimer.melting_temperature(method=MeltingTemperature.SALT_ADJUSTED)
    assert round(salt_adjusted_melt_temp, 0) == 65.0


def test_melt_temp_4():  # Reference: http://biotools.nubic.northwestern.edu/OligoCalc.html
    oprimer = Primer("CTGGAGGACGGAAGAGGAAGTAAA", salt=65)
    salt_adjusted_melt_temp = oprimer.melting_temperature(method=MeltingTemperature.SALT_ADJUSTED)
    assert round(salt_adjusted_melt_temp, 0) == 67.0


def test_melt_temp_5():  # Reference: http://biotools.nubic.northwestern.edu/OligoCalc.html
    oprimer = Primer("CTGGAGG")
    salt_adjusted_melt_temp = oprimer.melting_temperature(method=MeltingTemperature.SALT_ADJUSTED)
    assert round(salt_adjusted_melt_temp, 0) == 24.0


def test_melt_temp_6():  # Reference: http://biotools.nubic.northwestern.edu/OligoCalc.html
    oprimer = Primer("CTGGAGG", salt=65)
    salt_adjusted_melt_temp = oprimer.melting_temperature(method=MeltingTemperature.SALT_ADJUSTED)
    assert round(salt_adjusted_melt_temp, 0) == 26.0


def test_melt_temp_7():  
    # References: 
    # https://www.sigmaaldrich.com/CA/en/technical-documents/protocol/genomics/pcr/oligos-melting-temp
    # https://www.sigmaaldrich.com/deepweb/assets/sigmaaldrich/marketing/global/documents/367/000/meltingtemp1.pdf
    # https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("AAAAACCCCCGGGGGTTTTT", salt=50)
    nearest_neighbor_melt_temp = oprimer.melting_temperature(method=MeltingTemperature.NEAREST_NEIGHBOR)
    assert round(nearest_neighbor_melt_temp, 1) == 69.6


def test_melt_temp_8():  
    # References: 
    # https://www.sigmaaldrich.com/CA/en/technical-documents/protocol/genomics/pcr/oligos-melting-temp
    # https://www.sigmaaldrich.com/deepweb/assets/sigmaaldrich/marketing/global/documents/367/000/meltingtemp1.pdf
    # https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("ACGTACGTACGTACGTACGT", salt=50)
    nearest_neighbor_melt_temp = oprimer.melting_temperature(method=MeltingTemperature.NEAREST_NEIGHBOR)
    assert round(nearest_neighbor_melt_temp, 1) == 57.1


def test_melt_temp_9():  
    # References: 
    # https://www.sigmaaldrich.com/CA/en/technical-documents/protocol/genomics/pcr/oligos-melting-temp
    # https://www.sigmaaldrich.com/deepweb/assets/sigmaaldrich/marketing/global/documents/367/000/meltingtemp1.pdf
    # https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("GATCGATCGATCGATCGATC", salt=50)
    nearest_neighbor_melt_temp = oprimer.melting_temperature(method=MeltingTemperature.NEAREST_NEIGHBOR)
    assert round(nearest_neighbor_melt_temp, 1) == 64.5


def test_melt_temp_10():
    # References: 
    # https://www.sigmaaldrich.com/CA/en/technical-documents/protocol/genomics/pcr/oligos-melting-temp
    # https://www.sigmaaldrich.com/deepweb/assets/sigmaaldrich/marketing/global/documents/367/000/meltingtemp1.pdf
    # https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("ATATATATATCGCGCGCGCG", salt=50)
    nearest_neighbor_melt_temp = oprimer.melting_temperature(method=MeltingTemperature.NEAREST_NEIGHBOR)
    assert round(nearest_neighbor_melt_temp, 1) == 66.3


def test_thermodynamic_constants_1():
    # References: 
    # https://www.sigmaaldrich.com/CA/en/technical-documents/protocol/genomics/pcr/oligos-melting-temp
    # https://www.sigmaaldrich.com/deepweb/assets/sigmaaldrich/marketing/global/documents/367/000/meltingtemp1.pdf
    # https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("AAAAACCCCCGGGGGTTTTT", salt=50)
    assert round(oprimer.delta_h, 1) == -185.7 and round(oprimer.delta_s, 3) == -0.467


def test_thermodynamic_constants_2():
    # References: 
    # https://www.sigmaaldrich.com/CA/en/technical-documents/protocol/genomics/pcr/oligos-melting-temp
    # https://www.sigmaaldrich.com/deepweb/assets/sigmaaldrich/marketing/global/documents/367/000/meltingtemp1.pdf
    # https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("ACGTACGTACGTACGTACGT", salt=50)
    assert round(oprimer.delta_h, 1) == -148.5 and round(oprimer.delta_s, 3) == -0.38


def test_thermodynamic_constants_3():
    # References: 
    # https://www.sigmaaldrich.com/CA/en/technical-documents/protocol/genomics/pcr/oligos-melting-temp
    # https://www.sigmaaldrich.com/deepweb/assets/sigmaaldrich/marketing/global/documents/367/000/meltingtemp1.pdf
    # https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("GATCGATCGATCGATCGATC", salt=50)
    assert round(oprimer.delta_h, 1) == -146.6 and round(oprimer.delta_s, 3) == -0.366


def test_thermodynamic_constants_4():
    # References: 
    # https://www.sigmaaldrich.com/CA/en/technical-documents/protocol/genomics/pcr/oligos-melting-temp
    # https://www.sigmaaldrich.com/deepweb/assets/sigmaaldrich/marketing/global/documents/367/000/meltingtemp1.pdf
    # https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("ATATATATATCGCGCGCGCG", salt=50)
    assert round(oprimer.delta_h, 1) == -176.5 and round(oprimer.delta_s, 3) == -0.446


def test_single_runs_1():  # Reference: https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("ATCGATCG")
    runs = oprimer.single_runs
    assert runs['A'] == 1 and runs['T'] == 1 and runs['C'] == 1 and runs['G'] == 1


def test_single_runs_2():  # Reference: https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("ATTCGATCCCCG")
    runs = oprimer.single_runs
    assert runs['A'] == 1 and runs['T'] == 2 and runs['C'] == 4 and runs['G'] == 1


def test_single_runs_3():  # Reference: https://www.oligoevaluator.com/OligoCalcServlet
    oprimer = Primer("AAAAATTCGGGGATCCCCG")
    runs = oprimer.single_runs
    assert runs['A'] == 5 and runs['T'] == 2 and runs['C'] == 4 and runs['G'] == 4


def test_double_runs():
    p1 = Primer("ATATCGAACACACACACA")
    double_runs = p1.double_runs
    true_answer = {
        'AT': 2,
        'AC': 5,
        'AG': 0,
        'TA': 1,
        'TC': 1,
        'TG': 0,
        'CA': 5,
        'CT': 0,
        'CG': 1,
        'GA': 1,
        'GT': 0,
        'GC': 0
    }
    assert len(true_answer) == len(double_runs) and all(double_runs[pair] == true_answer[pair] for pair in double_runs)


def test_repeats_1():
    p = Primer("ATCG")
    assert (
        p.repeats(sequence="A", consecutive=False) == 1 and
        p.repeats(sequence="AT", consecutive=False) == 1 and
        p.repeats(sequence="AC", consecutive=False) == 0 and
        p.repeats(sequence="A", consecutive=True) == 1 and
        p.repeats(sequence="AT", consecutive=True) == 1
    )


def test_repeats_2():
    p = Primer("AAAATCGTGT")
    assert (
        p.repeats(sequence="AA", consecutive=True) == 2 and
        p.repeats(sequence="GT", consecutive=True) == 2
    )


def test_repeats_3():
    p = Primer("ATCGATCGATCG")
    assert p.repeats(sequence="ATCG", consecutive=True) == 3


def test_e260_1():  # https://atdbio.com/tools/oligo-calculator
    oprimer = Primer("ATCGATCGATCGATCGAT")
    assert round(oprimer.E260, 1) == 179.6


def test_e260_2():  # https://atdbio.com/tools/oligo-calculator
    oprimer = Primer("ACGT")
    assert round(oprimer.E260, 1) == 40.3
