import json
import subprocess
import pandas as pd
from pathlib import Path
import pytest
import threading
import time
import sys


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
def run_summary(integration_config, tmp_path_factory):
    gtf   = integration_config["gtf"]
    fasta = integration_config["fasta"]

    if not Path(gtf).exists():
        pytest.skip(f"GTF not found: {gtf}")
    if not Path(fasta).exists():
        pytest.skip(f"FASTA not found: {fasta}")

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
         "--gtf",    gtf,
         "--fa",     fasta,
         "--outdir", str(outdir)],
        capture_output=True, text=True
    )

    stop_event.set()
    progress_thread.join()

    assert result.returncode == 0, f"ORFannotate failed:\n{result.stderr}"
    return pd.read_csv(outdir / "ORFannotate_summary.tsv", sep="\t", dtype={"chrom": str})


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
