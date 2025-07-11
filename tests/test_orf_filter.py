import sys, os, io, tempfile
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from orfannotate.orf_filter import get_best_orfs_by_cpat

def test_get_best_orfs_by_cpat_selects_highest_scoring_orf():
    mock_cpat_tsv = """ID\tORF\tORF_start\tORF_end\tORF_strand\tORF_frame\tCoding_prob
TX001_ORF_1\t150\t1\t150\t+\t0\t0.12
TX001_ORF_2\t300\t10\t309\t+\t0\t0.85
TX002_ORF_1\t180\t5\t184\t+\t0\t0.67
TX002_ORF_2\t190\t4\t193\t+\t0\t0.35
"""

    # Write mock TSV to a temp file so we can pass a real path
    with tempfile.NamedTemporaryFile("w+", delete=False) as tmp:
        tmp.write(mock_cpat_tsv)
        tmp_path = tmp.name

    # Run the function under test
    best_orfs = get_best_orfs_by_cpat(cpat_best_path=tmp_path)

    # TX001 should pick ORF_2 (coding_prob 0.85, start 10)
    assert best_orfs["TX001"]["start"] == 10
    assert best_orfs["TX001"]["coding_prob"] == 0.85

    # TX002 should pick ORF_1 (coding_prob 0.67, start 5)
    assert best_orfs["TX002"]["start"] == 5
    assert best_orfs["TX002"]["coding_prob"] == 0.67
