"""Tests for rule-based QC decisions."""

from __future__ import annotations

from autobiopipe.config import AutoBioPipeConfig, QCConfig
from autobiopipe.decision import DecisionEngine
from autobiopipe.biology import analyze_biology
from autobiopipe.parser import FastaRecord, FastqRecord
from autobiopipe.qc import run_qc


def _config() -> AutoBioPipeConfig:
    return AutoBioPipeConfig(
        qc=QCConfig(
            min_avg_quality=20.0,
            gc_content_min=35.0,
            gc_content_max=65.0,
            min_read_length=4,
            length_cv_threshold=0.20,
            n_fraction_warning_threshold=0.05,
            n_fraction_critical_threshold=0.20,
            small_dataset_threshold=100,
            enable_small_dataset_warning=True,
        )
    )


def test_pass_when_metrics_are_within_thresholds() -> None:
    records = [FastqRecord("r1", "", "ATGC", "IIII")] * 120
    qc_result = run_qc(records, "reads.fastq", "fastq", quality_threshold=20.0)
    decision = DecisionEngine(_config()).evaluate(qc_result)
    assert decision.overall_status == "PASS"


def test_fail_for_low_quality() -> None:
    records = [FastqRecord("r1", "", "ATGC", "!!!!")] * 120
    qc_result = run_qc(records, "reads.fastq", "fastq", quality_threshold=20.0)
    decision = DecisionEngine(_config()).evaluate(qc_result)
    assert decision.overall_status == "FAIL"
    assert any(f.rule_id == "QC001" for f in decision.findings)


def test_warning_for_gc_anomaly() -> None:
    records = [FastaRecord("s1", "", "AAAA")] * 120
    qc_result = run_qc(records, "reads.fasta", "fasta")
    decision = DecisionEngine(_config()).evaluate(qc_result)
    assert any(f.rule_id == "QC002" for f in decision.findings)


def test_warning_for_length_variability() -> None:
    records = [FastaRecord("s1", "", "A"), FastaRecord("s2", "", "ATGCATGC")] * 60
    qc_result = run_qc(records, "reads.fasta", "fasta")
    decision = DecisionEngine(_config()).evaluate(qc_result)
    assert any(f.rule_id == "QC003" for f in decision.findings)


def test_fail_for_high_n_fraction() -> None:
    records = [FastaRecord("s1", "", "NNNN")] * 120
    qc_result = run_qc(records, "reads.fasta", "fasta")
    decision = DecisionEngine(_config()).evaluate(qc_result)
    assert any(f.rule_id == "QC004" for f in decision.findings)
    assert decision.overall_status == "FAIL"


def test_info_for_small_dataset() -> None:
    records = [FastaRecord("s1", "", "ATGC")] * 5
    qc_result = run_qc(records, "reads.fasta", "fasta")
    decision = DecisionEngine(_config()).evaluate(qc_result)
    assert any(f.rule_id == "QC005" for f in decision.findings)


def test_small_dataset_rule_can_be_disabled() -> None:
    config = _config()
    config.qc.enable_small_dataset_warning = False
    records = [FastaRecord("s1", "", "ATGC")] * 5
    qc_result = run_qc(records, "reads.fasta", "fasta")
    decision = DecisionEngine(config).evaluate(qc_result)
    assert all(f.rule_id != "QC005" for f in decision.findings)


def test_rule_qc014_triggers_for_duplicate_sequences() -> None:
    config = _config()
    records = [FastaRecord("s1", "", "ATGAAACCCGGGTAA")] * 20
    qc_result = run_qc(records, "reads.fasta", "fasta")
    qc_result.biology = analyze_biology(records, config.biology)
    qc_result.sequence_scale = qc_result.biology.sequence_scale
    decision = DecisionEngine(config).evaluate(qc_result)
    assert any(f.rule_id == "QC014" for f in decision.findings)


def test_rule_qc015_triggers_for_gc_skew_imbalance() -> None:
    config = _config()
    records = [FastaRecord("s1", "", "GGGGGGGGGGAAAAAA")] * 20
    qc_result = run_qc(records, "reads.fasta", "fasta")
    qc_result.biology = analyze_biology(records, config.biology)
    qc_result.sequence_scale = qc_result.biology.sequence_scale
    decision = DecisionEngine(config).evaluate(qc_result)
    assert any(f.rule_id == "QC015" for f in decision.findings)
