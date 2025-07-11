import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def test_import_nmd():
    from orfannotate import nmd
    assert hasattr(nmd, "predict_nmd")

def test_import_kozak():
    from orfannotate import kozak
    assert hasattr(kozak, "score_kozak")

def test_import_orf_filter():
    from orfannotate import orf_filter
    assert hasattr(orf_filter, "get_best_orfs_by_cpat")
    assert hasattr(orf_filter, "build_cds_features")

def test_import_gtf_annotation():
    from orfannotate import gtf_annotation
    assert hasattr(gtf_annotation, "annotate_gtf_with_cds")

def test_import_transcript_extraction():
    from orfannotate import transcript_extraction
    assert hasattr(transcript_extraction, "extract_transcripts_from_gtf")

def test_import_orf_prediction():
    from orfannotate import orf_prediction
    assert hasattr(orf_prediction, "run_cpat")

def test_import_summarise():
    from orfannotate import summarise
    assert hasattr(summarise, "generate_summary")
