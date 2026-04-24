"""Tests for biological interpretation helpers."""

from __future__ import annotations

from autobiopipe.biology import (
    GCReferenceDatabase,
    analyze_biology,
    classify_sequence_scale,
    estimate_coding_potential,
    evaluate_gc_context,
    infer_organism_type,
)
from autobiopipe.config import BiologyConfig
from autobiopipe.parser import FastaRecord


def _config() -> BiologyConfig:
    return BiologyConfig(
        gene_max_length=10_000,
        genome_min_length=1_000_000,
        coding_score_coding_threshold=0.65,
        coding_score_noncoding_threshold=0.35,
        low_complexity_entropy_threshold=0.55,
        homopolymer_run_threshold=8,
    )


def test_gc_reference_database_contains_major_groups() -> None:
    database = GCReferenceDatabase()
    assert database.get_range("bacteria").min_fraction == 0.30
    assert database.get_range("human").max_fraction == 0.44
    assert database.get_range("virus").description


def test_classify_sequence_scale_gene() -> None:
    assert classify_sequence_scale(5000, _config()) == "gene/transcript"


def test_classify_sequence_scale_contig() -> None:
    assert classify_sequence_scale(500_000, _config()) == "contig"


def test_classify_sequence_scale_genome() -> None:
    assert classify_sequence_scale(2_000_000, _config()) == "genome"


def test_infer_organism_type_small_sequence() -> None:
    assert infer_organism_type(45.0, 5000, _config()) == "likely_gene_or_viral"


def test_evaluate_gc_context_normal_for_human_like_gc() -> None:
    evaluation = evaluate_gc_context(41.0, "likely_human_genome")
    assert evaluation.status == "normal"


def test_estimate_coding_potential_for_simple_orf() -> None:
    result = estimate_coding_potential("ATGAAACCCGGGTAA", _config())
    assert result.score > 0.65
    assert result.classification == "coding"


def test_estimate_coding_potential_for_non_coding_repeat() -> None:
    result = estimate_coding_potential("AAAAAAAAAAAAAA", _config())
    assert result.classification == "non-coding"


def test_analyze_biology_summarizes_dataset() -> None:
    records = [
        FastaRecord("s1", "", "ATGAAACCCGGGTAA"),
        FastaRecord("s2", "", "ATGAAACCCGGGTAA"),
        FastaRecord("s3", "", "AAAAAAAAAAAAAAA"),
    ]
    summary = analyze_biology(records, _config())
    assert summary.sequence_scale == "gene/transcript"
    assert summary.inferred_organism_type
    assert summary.duplicate_fraction > 0.0
