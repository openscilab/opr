from opr import Primer

def test_reverse_1():
    oprimer = Primer("ATCGATCGATCGATCGAT")
    oprimer_reversed = oprimer.reverse(inplace=False)
    assert oprimer_reversed.sequence == "TAGCTAGCTAGCTAGCTA"

def test_reverse_2():
    oprimer = Primer("ATCGATCGATCGATCGAT")
    oprimer_reversed = oprimer.reverse()
    assert oprimer_reversed.sequence == "TAGCTAGCTAGCTAGCTA"

def test_reverse_3():
    oprimer = Primer("ATCGATCGATCGATCGAT")
    oprimer.reverse(inplace=True)
    assert oprimer.sequence == "TAGCTAGCTAGCTAGCTA"


def test_complement_1():
    oprimer = Primer("ATCGGCTAAATCGGCTAA")
    oprimer_complemented = oprimer.complement(inplace=False)
    assert oprimer_complemented.sequence == "TAGCCGATTTAGCCGATT"

def test_complement_2():
    oprimer = Primer("ATCGGCTAAATCGGCTAA")
    oprimer_complemented = oprimer.complement()
    assert oprimer_complemented.sequence == "TAGCCGATTTAGCCGATT"

def test_complement_3():
    oprimer = Primer("ATCGGCTAAATCGGCTAA")
    oprimer.complement(inplace=True)
    assert oprimer.sequence == "TAGCCGATTTAGCCGATT"

