from tests.getrm_harness import (
    SequencingDataUnavailable,
    bcftools_available,
    fetch_1000genomes_region,
    load_consensus_table,
    run_concordance_check,
    summarize,
)


def test_consensus_table_loads_and_has_expected_shape():
    samples = load_consensus_table()
    assert len(samples) >= 10
    for entry in samples:
        assert entry["gene"] == "CYP2D6"
        assert entry["coriell_id"].startswith("NA")
        assert "/" in entry["consensus_diplotype"]


def test_fetch_without_bcftools_raises_actionable_error():
    if bcftools_available():
        import pytest
        pytest.skip("bcftools is installed in this environment; this test covers its absence")
    try:
        fetch_1000genomes_region("22", 42522500, 42526900, "NA17244")
        assert False, "expected SequencingDataUnavailable"
    except SequencingDataUnavailable as exc:
        assert "bcftools" in str(exc)


def test_unconfigured_chromosome_raises_clear_error():
    if not bcftools_available():
        import pytest
        pytest.skip("bcftools absence is covered by test_fetch_without_bcftools_raises_actionable_error")
    try:
        fetch_1000genomes_region("99", 1, 100, "NA17244")
        assert False, "expected SequencingDataUnavailable"
    except SequencingDataUnavailable as exc:
        assert "No indexed 1000 Genomes URL" in str(exc)


def test_run_concordance_check_never_fabricates_a_match_without_data():
    outcomes = run_concordance_check()
    assert len(outcomes) >= 10
    # Every outcome must be honestly reported -- since end-to-end calling
    # isn't wired up yet (see getrm_harness.py docstring), nothing may be
    # reported as "match"/"mismatch" here; everything is "skipped_no_data".
    assert all(o.status == "skipped_no_data" for o in outcomes)


def test_summarize_reports_skipped_not_fabricated_concordance():
    outcomes = run_concordance_check()
    summary = summarize(outcomes)
    assert summary["matches"] == 0
    assert summary["concordance_pct_of_attempted"] is None
    assert summary["skipped_no_data"] == summary["total"]
