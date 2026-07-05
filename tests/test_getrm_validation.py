import pytest

from tests.getrm_harness import (
    SequencingDataUnavailable,
    bcftools_available,
    fetch_1000genomes_panel_samples,
    fetch_1000genomes_region,
    load_consensus_table,
    run_concordance_check,
    summarize,
)

# Samples confirmed (via the real 1000 Genomes panel file, cross-checked
# against the actual VCF header) to NOT be in the 1000 Genomes phase 3
# panel -- see validation/getrm_concordance_report.md for the full table.
KNOWN_NOT_IN_1000G = {"NA17222", "NA17287", "NA18524", "NA19790", "NA17244"}
KNOWN_IN_1000G = {"NA18563", "NA18564", "NA18565", "NA19785", "NA19908", "NA20289"}


def test_consensus_table_loads_and_has_expected_shape():
    samples = load_consensus_table()
    assert len(samples) >= 10
    for entry in samples:
        assert entry["gene"] == "CYP2D6"
        assert entry["coriell_id"].startswith("NA")
        assert "/" in entry["consensus_diplotype"]


def test_fetch_without_bcftools_raises_actionable_error():
    if bcftools_available():
        pytest.skip("bcftools is installed in this environment; this test covers its absence")
    try:
        fetch_1000genomes_region("22", 42500000, 42550000, "NA18563")
        assert False, "expected SequencingDataUnavailable"
    except SequencingDataUnavailable as exc:
        assert "bcftools" in str(exc)


def test_unconfigured_chromosome_raises_clear_error():
    if not bcftools_available():
        pytest.skip("bcftools absence is covered by test_fetch_without_bcftools_raises_actionable_error")
    try:
        fetch_1000genomes_region("99", 1, 100, "NA18563")
        assert False, "expected SequencingDataUnavailable"
    except SequencingDataUnavailable as exc:
        assert "No indexed 1000 Genomes URL" in str(exc)


def _panel_reachable() -> bool:
    try:
        fetch_1000genomes_panel_samples()
        return True
    except SequencingDataUnavailable:
        return False


def test_run_concordance_check_never_fabricates_a_result():
    """Regardless of network/tooling availability, every outcome must be
    an honestly-earned match/mismatch, or an explicitly reasoned skip --
    never a match/mismatch that wasn't actually computed from real data.
    """
    outcomes = run_concordance_check()
    assert len(outcomes) >= 10
    assert all(o.status in ("match", "mismatch", "skipped_no_data") for o in outcomes)
    for o in outcomes:
        if o.status in ("match", "mismatch"):
            assert o.called_diplotype is not None
            assert o.n_variants_considered > 0
        else:
            assert o.detail  # a skip must always carry a reason


def test_samples_not_in_1000g_panel_are_skipped_with_correct_reason():
    if not _panel_reachable():
        pytest.skip("1000 Genomes panel file not reachable in this environment")
    outcomes = {o.coriell_id: o for o in run_concordance_check()}
    for sample_id in KNOWN_NOT_IN_1000G:
        assert outcomes[sample_id].status == "skipped_no_data"
        assert "not present in the 1000 Genomes" in outcomes[sample_id].detail


def test_samples_in_1000g_panel_get_a_real_diplotype_call():
    if not _panel_reachable() or not bcftools_available():
        pytest.skip("1000 Genomes panel/bcftools not reachable in this environment")
    outcomes = {o.coriell_id: o for o in run_concordance_check()}
    for sample_id in KNOWN_IN_1000G:
        outcome = outcomes[sample_id]
        # A real fetch+call was attempted -- result is match or mismatch,
        # never skipped, since these are confirmed panel members.
        assert outcome.status in ("match", "mismatch"), (sample_id, outcome.detail)


def test_summarize_reflects_real_or_honest_skip_counts():
    outcomes = run_concordance_check()
    summary = summarize(outcomes)
    assert summary["total"] == len(outcomes)
    assert summary["matches"] + summary["mismatches"] + summary["skipped_no_data"] == summary["total"]
    if summary["matches"] + summary["mismatches"] == 0:
        assert summary["concordance_pct_of_attempted"] is None
