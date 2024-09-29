from opr import Primer

TEST_CASE_NAME = "operations testcase"

def test_reverse_1(): #Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGATCGATCGATCGAT")
    oprimer_reversed = oprimer.reverse(inplace=False)
    assert oprimer_reversed.sequence == "TAGCTAGCTAGCTAGCTA"

def test_reverse_2(): #Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGATCGATCGATCGAT")
    oprimer_reversed = oprimer.reverse()
    assert oprimer_reversed.sequence == "TAGCTAGCTAGCTAGCTA"

def test_reverse_3(): #Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGATCGATCGATCGAT")
    oprimer.reverse(inplace=True)
    assert oprimer.sequence == "TAGCTAGCTAGCTAGCTA"


def test_complement_1(): #Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGGCTAAATCGGCTAA")
    oprimer_complemented = oprimer.complement(inplace=False)
    assert oprimer_complemented.sequence == "TAGCCGATTTAGCCGATT"

def test_complement_2(): #Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
    oprimer = Primer("ATCGGCTAAATCGGCTAA")
    oprimer_complemented = oprimer.complement()
    assert oprimer_complemented.sequence == "TAGCCGATTTAGCCGATT"

def test_complement_3(): #Reference: https://www.qiagen.com/us/applications/enzymes/tools-and-calculators/reverse-complement-converter
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
