import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from orfannotate.kozak import score_kozak

def test_strong_kozak():
    seq   = "AAAAGGATGGGGAA"
    start = 9
    strength, _ = score_kozak(seq, start)
    assert strength == "strong"

def test_moderate_kozak_minus3():
    seq   = "TGCCTGGAAGATGGAGC"
    start = 9
    strength, _ = score_kozak(seq, start)
    assert strength == "moderate"

def test_weak_kozak():
    seq   = "TTTTTTTTTTATGCCC"
    start = 11
    strength, _ = score_kozak(seq, start)
    assert strength == "weak"

def test_short_sequence_returns_na():
    seq   = "ATG"
    start = 2
    strength, context = score_kozak(seq, start)
    assert strength == "weak"
    assert context == "NA"