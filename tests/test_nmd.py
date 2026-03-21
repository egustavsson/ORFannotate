from orfannotate.nmd import predict_nmd

def test_nmd_positive_strand_true():
    # EJ at end of exon1 = 160, stop at 100 → 160-100=60 > 50 → TRUE
    assert predict_nmd(orf_end=100, junctions=[(50, 160), (210, 300)], strand="+") == "TRUE"

def test_nmd_positive_strand_false():
    # EJ at end of exon1 = 210, stop at 180 → 210-180=30 < 50 → FALSE
    assert predict_nmd(orf_end=180, junctions=[(100, 210), (260, 300)], strand="+") == "FALSE"

def test_nmd_negative_strand_true():
    # EJ at start of exon2 = 240, stop at 300 → 300-240=60 > 50 → TRUE
    assert predict_nmd(orf_end=300, junctions=[(100, 200), (240, 400)], strand="-") == "TRUE"

def test_nmd_negative_strand_false():
    # EJ at start of exon2 = 200, stop at 230 → 230-200=30 < 50 → FALSE
    assert predict_nmd(orf_end=230, junctions=[(100, 180), (200, 400)], strand="-") == "FALSE"

def test_nmd_no_junctions():
    # Empty list → FALSE
    assert predict_nmd(orf_end=100, junctions=[], strand="+") == "FALSE"

def test_nmd_mono_exonic():
    # Single exon → len < 2 → FALSE regardless of coordinates
    assert predict_nmd(orf_end=100, junctions=[(1, 300)], strand="+") == "FALSE"

def test_nmd_junction_only_in_cds():
    # EJ at 80, stop at 100 → 80 not > 100 → no downstream EJ → FALSE
    assert predict_nmd(orf_end=100, junctions=[(1, 80), (130, 300)], strand="+") == "FALSE"

def test_nmd_multiple_utr3_junctions_last_triggers():
    # Two downstream EJs at 130 and 200, stop at 100
    # last downstream EJ = 200 → 200-100=100 > 50 → TRUE
    assert predict_nmd(orf_end=100, junctions=[(1, 80), (90, 130), (150, 200), (250, 300)], strand="+") == "TRUE"

def test_nmd_multiple_utr3_junctions_last_does_not_trigger():
    # Two downstream EJs at 130 and 120, stop at 100
    # last downstream EJ = 130 → 130-100=30 < 50 → FALSE
    assert predict_nmd(orf_end=100, junctions=[(1, 80), (90, 120), (140, 300)], strand="+") == "FALSE"