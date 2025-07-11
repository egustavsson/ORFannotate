from orfannotate.nmd import predict_nmd

def test_nmd_positive_strand_true():
    assert predict_nmd(orf_end=100, junctions=[(160, 200)], strand="+") == "TRUE"

def test_nmd_positive_strand_false():
    assert predict_nmd(orf_end=180, junctions=[(200, 250)], strand="+") == "FALSE"

def test_nmd_negative_strand_true():
    assert predict_nmd(orf_end=300, junctions=[(100, 240)], strand="-") == "TRUE"

def test_nmd_negative_strand_false():
    assert predict_nmd(orf_end=230, junctions=[(100, 240)], strand="-") == "FALSE"

def test_nmd_no_junctions():
    assert predict_nmd(orf_end=100, junctions=[], strand="+") == "FALSE"
