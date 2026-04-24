"""Quality control metrics for sequence datasets."""

from __future__ import annotations

import statistics
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Dict, Iterable, Optional

from autobiopipe.parser import AnyRecord, FastqRecord

if TYPE_CHECKING:
    from autobiopipe.biology import BiologicalSummary


@dataclass
class GCStats:
    """Summary statistics for GC content."""

    mean: float
    std_dev: float
    min_value: float
    max_value: float
    histogram: Dict[int, int]
    warnings: list[str] = field(default_factory=list)


@dataclass
class LengthStats:
    """Summary statistics for read or sequence lengths."""

    mean: float
    std_dev: float
    min_value: int
    max_value: int
    median: float
    coefficient_of_variation: float
    total_bases: int
    histogram: Dict[int, int]
    warnings: list[str] = field(default_factory=list)


@dataclass
class QualityStats:
    """Summary statistics for FASTQ quality scores."""

    mean_per_read: float
    std_dev_per_read: float
    min_avg_quality: float
    max_avg_quality: float
    reads_below_threshold: int
    reads_below_threshold_pct: float
    per_position_means: list[float]
    threshold_used: float
    warnings: list[str] = field(default_factory=list)


@dataclass
class QCResult:
    """Complete QC output for one sequence file."""

    file_name: str
    file_type: str
    total_records: int
    gc: GCStats
    length: LengthStats
    quality: Optional[QualityStats]
    n_fraction: float
    sequence_scale: str = "unknown"
    biology: Optional["BiologicalSummary"] = None
    warnings: list[str] = field(default_factory=list)


def _safe_std(values: Iterable[float]) -> float:
    """Return the sample standard deviation or 0.0 for short inputs.

    Args:
        values: Numeric values.

    Returns:
        Sample standard deviation, or 0.0 when fewer than two values exist.
    """
    materialized = list(values)
    return statistics.stdev(materialized) if len(materialized) > 1 else 0.0


def _build_histogram(values: Iterable[float], step: int) -> dict[int, int]:
    """Bucket numeric values into fixed-width integer bins.

    Args:
        values: Numeric values to bin.
        step: Bin width.

    Returns:
        Histogram keyed by bin start.
    """
    histogram: dict[int, int] = {}
    for value in values:
        bucket = int(value // step) * step
        histogram[bucket] = histogram.get(bucket, 0) + 1
    return dict(sorted(histogram.items()))


def compute_gc_stats(records: list[AnyRecord]) -> GCStats:
    """Compute GC-content summary statistics.

    Args:
        records: Parsed FASTA or FASTQ records.

    Returns:
        GCStats including a 5% histogram.
    """
    if not records:
        return GCStats(0.0, 0.0, 0.0, 0.0, {}, ["No records available for GC analysis."])

    gc_values = [record.gc_content() for record in records]
    warnings: list[str] = []
    mean_gc = statistics.mean(gc_values)
    if mean_gc < 20.0 or mean_gc > 80.0:
        warnings.append("Mean GC content is highly atypical.")
    if _safe_std(gc_values) > 20.0:
        warnings.append("GC-content variability is unusually high across records.")

    return GCStats(
        mean=mean_gc,
        std_dev=_safe_std(gc_values),
        min_value=min(gc_values),
        max_value=max(gc_values),
        histogram=_build_histogram(gc_values, step=5),
        warnings=warnings,
    )


def compute_length_stats(records: list[AnyRecord]) -> LengthStats:
    """Compute read or sequence length statistics.

    Args:
        records: Parsed FASTA or FASTQ records.

    Returns:
        LengthStats including coefficient of variation and histogram.
    """
    if not records:
        return LengthStats(
            mean=0.0,
            std_dev=0.0,
            min_value=0,
            max_value=0,
            median=0.0,
            coefficient_of_variation=0.0,
            total_bases=0,
            histogram={},
            warnings=["No records available for length analysis."],
        )

    lengths = [len(record) for record in records]
    mean_length = statistics.mean(lengths)
    std_dev = _safe_std([float(length) for length in lengths])
    warnings: list[str] = []
    if mean_length and (std_dev / mean_length) > 0.5:
        warnings.append("Read lengths are highly variable.")
    if min(lengths) < 10:
        warnings.append("Very short records are present.")

    return LengthStats(
        mean=mean_length,
        std_dev=std_dev,
        min_value=min(lengths),
        max_value=max(lengths),
        median=statistics.median(lengths),
        coefficient_of_variation=(std_dev / mean_length) if mean_length else 0.0,
        total_bases=sum(lengths),
        histogram=_build_histogram(lengths, step=max(1, max(lengths) // 10 or 1)),
        warnings=warnings,
    )


def compute_quality_stats(
    records: list[FastqRecord],
    quality_threshold: float = 20.0,
) -> QualityStats:
    """Compute FASTQ quality-score statistics.

    Args:
        records: Parsed FASTQ records.
        quality_threshold: Minimum acceptable average quality threshold.

    Returns:
        QualityStats including per-position means.
    """
    if not records:
        return QualityStats(
            mean_per_read=0.0,
            std_dev_per_read=0.0,
            min_avg_quality=0.0,
            max_avg_quality=0.0,
            reads_below_threshold=0,
            reads_below_threshold_pct=0.0,
            per_position_means=[],
            threshold_used=quality_threshold,
            warnings=["No FASTQ records available for quality analysis."],
        )

    per_read_means = [record.avg_quality() for record in records]
    max_length = max(len(record) for record in records)
    position_sums = [0.0] * max_length
    position_counts = [0] * max_length

    for record in records:
        for index, score in enumerate(record.quality_scores):
            position_sums[index] += score
            position_counts[index] += 1

    per_position_means = [
        position_sums[index] / position_counts[index]
        for index in range(max_length)
        if position_counts[index] > 0
    ]
    reads_below_threshold = sum(1 for score in per_read_means if score < quality_threshold)
    reads_below_threshold_pct = (reads_below_threshold / len(records)) * 100.0
    warnings: list[str] = []
    if statistics.mean(per_read_means) < quality_threshold:
        warnings.append("Average read quality is below the configured threshold.")
    if reads_below_threshold_pct > 10.0:
        warnings.append("More than 10% of reads fall below the quality threshold.")

    return QualityStats(
        mean_per_read=statistics.mean(per_read_means),
        std_dev_per_read=_safe_std(per_read_means),
        min_avg_quality=min(per_read_means),
        max_avg_quality=max(per_read_means),
        reads_below_threshold=reads_below_threshold,
        reads_below_threshold_pct=reads_below_threshold_pct,
        per_position_means=per_position_means,
        threshold_used=quality_threshold,
        warnings=warnings,
    )


def compute_n_fraction(records: list[AnyRecord]) -> float:
    """Compute the fraction of bases represented by ``N``.

    Args:
        records: Parsed FASTA or FASTQ records.

    Returns:
        Fraction of ambiguous ``N`` bases in the interval 0.0 to 1.0.
    """
    if not records:
        return 0.0
    total_bases = sum(len(record) for record in records)
    if total_bases == 0:
        return 0.0
    ambiguous_bases = sum(record.sequence.upper().count("N") for record in records)
    return ambiguous_bases / total_bases


def run_qc(
    records: list[AnyRecord],
    file_name: str,
    file_type: str,
    quality_threshold: float = 20.0,
) -> QCResult:
    """Run the complete QC metric suite.

    Args:
        records: Parsed FASTA or FASTQ records.
        file_name: Original input filename.
        file_type: Detected file type.
        quality_threshold: Minimum acceptable quality threshold for FASTQ.

    Returns:
        QCResult containing all computed metrics and warnings.
    """
    warnings: list[str] = []
    gc_stats = compute_gc_stats(records)
    length_stats = compute_length_stats(records)
    n_fraction = compute_n_fraction(records)
    warnings.extend(gc_stats.warnings)
    warnings.extend(length_stats.warnings)

    quality_stats: Optional[QualityStats] = None
    if file_type.lower() == "fastq":
        fastq_records = [record for record in records if isinstance(record, FastqRecord)]
        quality_stats = compute_quality_stats(fastq_records, quality_threshold=quality_threshold)
        warnings.extend(quality_stats.warnings)

    if n_fraction > 0.05:
        warnings.append(f"Ambiguous base fraction is elevated ({n_fraction:.2%}).")
    if not records:
        warnings.append("Input dataset contains no parseable records.")

    return QCResult(
        file_name=file_name,
        file_type=file_type,
        total_records=len(records),
        gc=gc_stats,
        length=length_stats,
        quality=quality_stats,
        n_fraction=n_fraction,
        warnings=warnings,
    )
