import json
import shutil
import subprocess
import pandas as pd
from pathlib import Path
import pytest
import threading
import time
import sys
import gffutils
import tempfile


DEFAULT_CONFIG = Path(__file__).parent / "gtf_configs" / "toy.json"


def load_config(config_path):
    with open(config_path) as f:
        return json.load(f)


def _stream_progress(stop_event, label="Running ORFannotate"):
    """Prints a dot every 30 seconds until stop_event is set."""
    start = time.time()
    interval = 30
    next_print = interval
    sys.stdout.write(f"\n  {label} ")
    sys.stdout.flush()
    while not stop_event.is_set():
        elapsed = time.time() - start
        if elapsed >= next_print:
            sys.stdout.write(".")
            sys.stdout.flush()
            next_print += interval
        time.sleep(1)
    elapsed = time.time() - start
    mins, secs = divmod(int(elapsed), 60)
    sys.stdout.write(f" done ({mins}m{secs:02d}s)\n\n")
    sys.stdout.flush()


@pytest.fixture(scope="module")
def integration_config(request):
    config_path = request.config.getoption("--run-config", default=None)
    path = Path(config_path) if config_path else DEFAULT_CONFIG
    if not path.exists():
        pytest.skip(f"Integration config not found: {path}")
    return load_config(path)


@pytest.fixture(scope="module")
def run_output(integration_config, tmp_path_factory):
    """Returns (summary_df, outdir) so tests can access both."""
    gtf           = integration_config["gtf"]
    fasta         = integration_config["fasta"]
    species       = integration_config["species"]
    coding_cutoff = integration_config["coding_cutoff"]

    if not Path(gtf).exists():
        pytest.fail(f"GTF not found: {gtf}")
    if not Path(fasta).exists():
        pytest.fail(f"FASTA not found: {fasta}")

    # Clear stale tmp folder to avoid cross-run contamination
    tmp_dir = Path("./tmp")
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
        sys.stdout.write("  Cleared stale ./tmp folder\n")
        sys.stdout.flush()
        
    outdir = tmp_path_factory.mktemp("orfannotate_integration")

    stop_event = threading.Event()
    progress_thread = threading.Thread(
        target=_stream_progress,
        args=(stop_event, f"Running ORFannotate on {Path(gtf).name}"),
        daemon=True
    )
    progress_thread.start()

    result = subprocess.run(
        ["python", "ORFannotate.py",
         "--gtf",           gtf,
         "--fa",            fasta,
         "--species",       species,
         "--coding-cutoff", str(coding_cutoff),
         "--outdir",        str(outdir)],
        capture_output=False, text=True
    )

    stop_event.set()
    progress_thread.join()

    assert result.returncode == 0, "ORFannotate failed — see output above"

    summary = pd.read_csv(
        outdir / "ORFannotate_summary.tsv",
        sep="\t",
        dtype={"chrom": str}
    )
    return summary, outdir


@pytest.fixture(scope="module")
def run_summary(run_output):
    return run_output[0]


@pytest.fixture(scope="module")
def run_outdir(run_output):
    return run_output[1]


# ---------------------------------------------------------------------------
# Config-driven tests
# ---------------------------------------------------------------------------

def test_known_coding_transcripts(run_summary, integration_config):
    entries = [e for e in integration_config.get("expected_fields", [])
               if e["fields"].get("coding_class") == "coding"]
    if not entries:
        pytest.skip("No coding transcripts in expected_fields")

    failures = []
    for entry in entries:
        tid = entry["tid"]
        matches = run_summary[run_summary.transcript_id == tid]
        if len(matches) != 1:
            failures.append(f"  [{tid}] not found in summary")
            continue
        row = matches.iloc[0]
        if row.coding_class != "coding":
            failures.append(f"  [{tid}] coding_class: expected coding, got {row.coding_class}")
        expected_aa = entry["fields"].get("orf_aa_len")
        if expected_aa is not None and row.orf_aa_len != expected_aa:
            failures.append(f"  [{tid}] orf_aa_len: expected {expected_aa}, got {row.orf_aa_len}")

    if failures:
        pytest.fail("Coding transcript mismatches:\n" + "\n".join(failures))


def test_known_noncoding_transcripts(run_summary, integration_config):
    entries = [e for e in integration_config.get("expected_fields", [])
               if e["fields"].get("coding_class") == "noncoding"]
    if not entries:
        pytest.skip("No noncoding transcripts in expected_fields")

    failures = []
    for entry in entries:
        tid = entry["tid"]
        matches = run_summary[run_summary.transcript_id == tid]
        if len(matches) != 1:
            failures.append(f"  [{tid}] not found in summary")
            continue
        row = matches.iloc[0]
        if row.coding_class != "noncoding":
            failures.append(f"  [{tid}] coding_class: expected noncoding, got {row.coding_class}")

    if failures:
        pytest.fail("Noncoding transcript mismatches:\n" + "\n".join(failures))


def test_expected_fields(run_summary, integration_config):
    FLOAT_TOLERANCES = {
        "coding_prob":           1e-6,
        "coding_prob_best_uORF": 1e-4,
    }

    failures = []
    for entry in integration_config.get("expected_fields", []):
        tid = entry["tid"]
        matches = run_summary[run_summary.transcript_id == tid]
        if len(matches) != 1:
            failures.append(f"  [{tid}] not found in summary")
            continue
        row = matches.iloc[0]

        for field, expected in entry["fields"].items():
            actual = row[field]
            if expected is None:
                if not pd.isna(actual):
                    failures.append(f"  [{tid}] {field}: expected NaN, got {actual!r}")
            elif isinstance(expected, float):
                tol = FLOAT_TOLERANCES.get(field, 1e-6)
                if pd.isna(actual) or abs(actual - expected) >= tol:
                    failures.append(f"  [{tid}] {field}: expected {expected}, got {actual}")
            else:
                if actual != expected:
                    failures.append(f"  [{tid}] {field}: expected {expected!r}, got {actual!r}")

    if failures:
        pytest.fail("Field mismatches:\n" + "\n".join(failures))


def test_below_cutoff_transcripts(run_summary, integration_config):
    """
    Transcripts with has_orf=True but coding_prob below the species cutoff
    should be classified as noncoding.
    """
    cutoff = integration_config.get("coding_cutoff")
    if cutoff is None:
        pytest.skip("No coding_cutoff defined in integration config")

    below_cutoff = [
        entry for entry in integration_config.get("expected_fields", [])
        if entry["fields"].get("coding_class") == "noncoding"
        and entry["fields"].get("has_orf") is True
        and entry["fields"].get("coding_prob") is not None
        and entry["fields"]["coding_prob"] < cutoff
    ]

    if not below_cutoff:
        pytest.skip("No below-cutoff transcripts found in expected_fields")

    failures = []
    for entry in below_cutoff:
        tid = entry["tid"]
        matches = run_summary[run_summary.transcript_id == tid]
        if len(matches) != 1:
            failures.append(f"  [{tid}] not found in summary")
            continue
        row = matches.iloc[0]

        if not (row.coding_prob < cutoff):
            failures.append(
                f"  [{tid}] coding_prob {row.coding_prob} is NOT below cutoff {cutoff}"
            )
        if row.coding_class != "noncoding":
            failures.append(
                f"  [{tid}] coding_class: expected noncoding, got {row.coding_class} "
                f"(coding_prob={row.coding_prob}, cutoff={cutoff})"
            )
        expected_prob = entry["fields"]["coding_prob"]
        if abs(row.coding_prob - expected_prob) >= 1e-4:
            failures.append(
                f"  [{tid}] coding_prob: expected {expected_prob}, got {row.coding_prob}"
            )

    if failures:
        pytest.fail("Below-cutoff transcript mismatches:\n" + "\n".join(failures))


# ---------------------------------------------------------------------------
# Structural integrity rules (apply to all coding transcripts in summary)
# ---------------------------------------------------------------------------

def test_nmd_requires_stop_to_last_ej_gt50(run_summary):
    """Rule: NMD_sensitive=TRUE only if stop_to_last_EJ > 50."""
    coding = run_summary[run_summary.coding_class == "coding"]
    failures = []
    for _, row in coding.iterrows():
        if row.NMD_sensitive == "TRUE" or row.NMD_sensitive is True:
            if pd.isna(row.stop_to_last_EJ) or row.stop_to_last_EJ <= 50:
                failures.append(
                    f"  [{row.transcript_id}] NMD_sensitive=TRUE but "
                    f"stop_to_last_EJ={row.stop_to_last_EJ}"
                )
    if failures:
        pytest.fail("NMD rule violated:\n" + "\n".join(failures))


def test_stop_to_last_ej_nonzero_requires_utr3_junctions(run_summary):
    """Rule: stop_to_last_EJ != 0 only if utr3_junctions > 0."""
    coding = run_summary[run_summary.coding_class == "coding"]
    failures = []
    for _, row in coding.iterrows():
        if pd.isna(row.stop_to_last_EJ):
            continue
        if row.stop_to_last_EJ != 0:
            if pd.isna(row.utr3_junctions) or row.utr3_junctions == 0:
                failures.append(
                    f"  [{row.transcript_id}] stop_to_last_EJ={row.stop_to_last_EJ} "
                    f"but utr3_junctions={row.utr3_junctions}"
                )
    if failures:
        pytest.fail("stop_to_last_EJ rule violated:\n" + "\n".join(failures))


def test_stop_to_last_ej_na_for_monoexonic(run_summary):
    """Rule: stop_to_last_EJ is NA or 0 for mono-exonic transcripts (total_junctions == 0)."""
    coding = run_summary[run_summary.coding_class == "coding"]
    failures = []
    for _, row in coding.iterrows():
        if row.total_junctions == 0:
            is_null = pd.isna(row.stop_to_last_EJ)
            is_zero = row.stop_to_last_EJ == 0
            if not (is_null or is_zero):
                failures.append(
                    f"  [{row.transcript_id}] mono-exonic but "
                    f"stop_to_last_EJ={row.stop_to_last_EJ} (expected NA or 0)"
                )
    if failures:
        pytest.fail("stop_to_last_EJ mono-exonic rule violated:\n" + "\n".join(failures))


def test_orf_nt_len_matches_annotated_gtf_cds(run_summary, run_outdir):
    """
    Rule: orf_nt_len must equal the total length of all CDS features
    in ORFannotate_annotated.gtf for that transcript.
    """
    annotated_gtf = run_outdir / "ORFannotate_annotated.gtf"
    if not annotated_gtf.exists():
        pytest.skip("ORFannotate_annotated.gtf not found in output directory")

    # Parse CDS lengths from annotated GTF directly (avoid gffutils overhead)
    cds_lengths = {}
    with open(annotated_gtf) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            # Extract transcript_id from attributes
            attrs = parts[8]
            tid = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("transcript_id"):
                    tid = attr.split('"')[1]
                    break
            if tid is None:
                continue
            cds_len = int(parts[4]) - int(parts[3]) + 1
            cds_lengths[tid] = cds_lengths.get(tid, 0) + cds_len

    coding = run_summary[run_summary.coding_class == "coding"]
    failures = []
    for _, row in coding.iterrows():
        tid = row.transcript_id
        if pd.isna(row.orf_nt_len):
            continue
        expected_len = cds_lengths.get(tid)
        if expected_len is None:
            failures.append(f"  [{tid}] no CDS features found in annotated GTF")
            continue
        if int(row.orf_nt_len) != expected_len:
            failures.append(
                f"  [{tid}] orf_nt_len={int(row.orf_nt_len)} but "
                f"sum of CDS lengths in GTF={expected_len}"
            )

    if failures:
        pytest.fail("orf_nt_len vs annotated GTF CDS mismatch:\n" + "\n".join(failures))
        
# import json
# import subprocess
# import pandas as pd
# from pathlib import Path
# import pytest
# import threading
# import time
# import sys


# DEFAULT_CONFIG = Path(__file__).parent / "gtf_configs" / "toy.json"


# def load_config(config_path):
#     with open(config_path) as f:
#         return json.load(f)


# def _stream_progress(stop_event, label="Running ORFannotate"):
#     """Prints a dot every 30 seconds until stop_event is set."""
#     start = time.time()
#     interval = 30
#     next_print = interval
#     sys.stdout.write(f"\n  {label} ")
#     sys.stdout.flush()
#     while not stop_event.is_set():
#         elapsed = time.time() - start
#         if elapsed >= next_print:
#             sys.stdout.write(".")
#             sys.stdout.flush()
#             next_print += interval
#         time.sleep(1)
#     elapsed = time.time() - start
#     mins, secs = divmod(int(elapsed), 60)
#     sys.stdout.write(f" done ({mins}m{secs:02d}s)\n\n")
#     sys.stdout.flush()


# @pytest.fixture(scope="module")
# def integration_config(request):
#     config_path = request.config.getoption("--run-config", default=None)
#     path = Path(config_path) if config_path else DEFAULT_CONFIG
#     if not path.exists():
#         pytest.skip(f"Integration config not found: {path}")
#     return load_config(path)


# @pytest.fixture(scope="module")
# def run_summary(integration_config, tmp_path_factory):
#     gtf         = integration_config["gtf"]
#     fasta       = integration_config["fasta"]
#     species     = integration_config["species"]
#     coding_cutoff = integration_config["coding_cutoff"]

#     if not Path(gtf).exists():
#         pytest.skip(f"GTF not found: {gtf}")
#     if not Path(fasta).exists():
#         pytest.skip(f"FASTA not found: {fasta}")

#     outdir = tmp_path_factory.mktemp("orfannotate_integration")

#     stop_event = threading.Event()
#     progress_thread = threading.Thread(
#         target=_stream_progress,
#         args=(stop_event, f"Running ORFannotate on {Path(gtf).name}"),
#         daemon=True
#     )
#     progress_thread.start()

#     result = subprocess.run(
#         ["python", "ORFannotate.py",
#          "--gtf",    gtf,
#          "--fa",     fasta,
#          "--species", species,
#          "--coding-cutoff", coding_cutoff,
#          "--outdir", str(outdir)],
#         capture_output=True, text=True
#     )

#     stop_event.set()
#     progress_thread.join()

#     assert result.returncode == 0, f"ORFannotate failed:\n{result.stderr}"
#     return pd.read_csv(outdir / "ORFannotate_summary.tsv", sep="\t", dtype={"chrom": str})


# def test_known_coding_transcripts(run_summary, integration_config):
#     entries = [e for e in integration_config.get("expected_fields", [])
#                if e["fields"].get("coding_class") == "coding"]
#     if not entries:
#         pytest.skip("No coding transcripts in expected_fields")

#     failures = []
#     for entry in entries:
#         tid = entry["tid"]
#         matches = run_summary[run_summary.transcript_id == tid]
#         if len(matches) != 1:
#             failures.append(f"  [{tid}] not found in summary")
#             continue
#         row = matches.iloc[0]
#         if row.coding_class != "coding":
#             failures.append(f"  [{tid}] coding_class: expected coding, got {row.coding_class}")
#         expected_aa = entry["fields"].get("orf_aa_len")
#         if expected_aa is not None and row.orf_aa_len != expected_aa:
#             failures.append(f"  [{tid}] orf_aa_len: expected {expected_aa}, got {row.orf_aa_len}")

#     if failures:
#         pytest.fail("Coding transcript mismatches:\n" + "\n".join(failures))


# def test_known_noncoding_transcripts(run_summary, integration_config):
#     entries = [e for e in integration_config.get("expected_fields", [])
#                if e["fields"].get("coding_class") == "noncoding"]
#     if not entries:
#         pytest.skip("No noncoding transcripts in expected_fields")

#     failures = []
#     for entry in entries:
#         tid = entry["tid"]
#         matches = run_summary[run_summary.transcript_id == tid]
#         if len(matches) != 1:
#             failures.append(f"  [{tid}] not found in summary")
#             continue
#         row = matches.iloc[0]
#         if row.coding_class != "noncoding":
#             failures.append(f"  [{tid}] coding_class: expected noncoding, got {row.coding_class}")

#     if failures:
#         pytest.fail("Noncoding transcript mismatches:\n" + "\n".join(failures))


# def test_expected_fields(run_summary, integration_config):
#     FLOAT_TOLERANCES = {
#         "coding_prob":           1e-6,
#         "coding_prob_best_uORF": 1e-4,
#     }

#     failures = []
#     for entry in integration_config.get("expected_fields", []):
#         tid = entry["tid"]
#         matches = run_summary[run_summary.transcript_id == tid]
#         if len(matches) != 1:
#             failures.append(f"  [{tid}] not found in summary")
#             continue
#         row = matches.iloc[0]

#         for field, expected in entry["fields"].items():
#             actual = row[field]
#             if expected is None:
#                 if not pd.isna(actual):
#                     failures.append(f"  [{tid}] {field}: expected NaN, got {actual!r}")
#             elif isinstance(expected, float):
#                 tol = FLOAT_TOLERANCES.get(field, 1e-6)
#                 if pd.isna(actual) or abs(actual - expected) >= tol:
#                     failures.append(f"  [{tid}] {field}: expected {expected}, got {actual}")
#             else:
#                 if actual != expected:
#                     failures.append(f"  [{tid}] {field}: expected {expected!r}, got {actual!r}")

#     if failures:
#         pytest.fail("Field mismatches:\n" + "\n".join(failures))


# def test_below_cutoff_transcripts(run_summary, integration_config):
#     """
#     Transcripts with has_orf=True but coding_prob below the species cutoff
#     should be classified as noncoding.
#     """
#     cutoff = integration_config.get("coding_cutoff")
#     if cutoff is None:
#         pytest.skip("No coding_cutoff defined in integration config")

#     below_cutoff = [
#         entry for entry in integration_config.get("expected_fields", [])
#         if entry["fields"].get("coding_class") == "noncoding"
#         and entry["fields"].get("has_orf") is True
#         and entry["fields"].get("coding_prob") is not None
#         and entry["fields"]["coding_prob"] < cutoff
#     ]

#     if not below_cutoff:
#         pytest.skip("No below-cutoff transcripts found in expected_fields")

#     failures = []
#     for entry in below_cutoff:
#         tid = entry["tid"]
#         matches = run_summary[run_summary.transcript_id == tid]
#         if len(matches) != 1:
#             failures.append(f"  [{tid}] not found in summary")
#             continue
#         row = matches.iloc[0]

#         if not (row.coding_prob < cutoff):
#             failures.append(
#                 f"  [{tid}] coding_prob {row.coding_prob} is NOT below cutoff {cutoff}"
#             )
#         if row.coding_class != "noncoding":
#             failures.append(
#                 f"  [{tid}] coding_class: expected noncoding, got {row.coding_class} "
#                 f"(coding_prob={row.coding_prob}, cutoff={cutoff})"
#             )
#         expected_prob = entry["fields"]["coding_prob"]
#         if abs(row.coding_prob - expected_prob) >= 1e-4:
#             failures.append(
#                 f"  [{tid}] coding_prob: expected {expected_prob}, got {row.coding_prob}"
#             )

#     if failures:
#         pytest.fail("Below-cutoff transcript mismatches:\n" + "\n".join(failures))

