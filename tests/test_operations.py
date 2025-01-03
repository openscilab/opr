from opr import Primer

TEST_CASE_NAME = "Operations tests"


def test_reverse_1():  # Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGATCGATCGATCGAT")
    oprimer_reversed = oprimer.reverse(inplace=False)
    assert oprimer_reversed.sequence == "TAGCTAGCTAGCTAGCTA"


def test_reverse_2():  # Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGATCGATCGATCGAT")
    oprimer_reversed = oprimer.reverse()
    assert oprimer_reversed.sequence == "TAGCTAGCTAGCTAGCTA"


def test_reverse_3():  # Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGATCGATCGATCGAT")
    oprimer.reverse(inplace=True)
    assert oprimer.sequence == "TAGCTAGCTAGCTAGCTA"


def test_complement_1():  # Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGGCTAAATCGGCTAA")
    oprimer_complemented = oprimer.complement(inplace=False)
    assert oprimer_complemented.sequence == "TAGCCGATTTAGCCGATT"


def test_complement_2():  # Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGGCTAAATCGGCTAA")
    oprimer_complemented = oprimer.complement()
    assert oprimer_complemented.sequence == "TAGCCGATTTAGCCGATT"


def test_complement_3():  # Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGGCTAAATCGGCTAA")
    oprimer.complement(inplace=True)
    assert oprimer.sequence == "TAGCCGATTTAGCCGATT"


def test_length():
    oprimer = Primer("ATCGGCTAAATCGGCTAA")
    assert len(oprimer) == 18


def test_addition():
    oprimer_1 = Primer("ATCG")
    oprimer_2 = Primer("GATC")
    oprimer_concat = oprimer_1 + oprimer_2
    assert oprimer_concat.sequence == "ATCGGATC"


def test_multiply():
    oprimer_1 = Primer("ATCG")
    oprimer_concat = oprimer_1 * 4
    assert oprimer_concat.sequence == "ATCGATCGATCGATCG"


def test_equality1():
    oprimer_1 = Primer("ATCG")
    oprimer_2 = Primer("ATCG")
    assert oprimer_1 == oprimer_2


def test_equality2():
    oprimer_1 = Primer("ATCG")
    oprimer_2 = Primer("ATCGC")
    assert oprimer_1 != oprimer_2


def test_str():
    oprimer_1 = Primer("ATCG")
    oprimer_2 = Primer("ATCGC")
    oprimer_concat = oprimer_1 + oprimer_2
    assert str(oprimer_1) + str(oprimer_2) == oprimer_concat.sequence
