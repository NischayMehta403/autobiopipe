"""Biological interpretation helpers for AutoBioPipe."""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, field

from autobiopipe.config import BiologyConfig
from autobiopipe.parser import AnyRecord

STOP_CODONS = {"TAA", "TAG", "TGA"}


@dataclass(frozen=True)
class GCReferenceRange:
    """Expected GC-content range for an organism group."""

    name: str
    min_fraction: float
    max_fraction: float
    description: str


class GCReferenceDatabase:
    """Reference GC-content ranges for broad organism groups."""

    def __init__(self) -> None:
        """Initialize the built-in reference ranges."""
        self._ranges: dict[str, GCReferenceRange] = {
            "bacteria": GCReferenceRange(
                name="bacteria",
                min_fraction=0.30,
                max_fraction=0.70,
                description="Broad bacterial genomes often span moderate to high GC content.",
            ),
            "human": GCReferenceRange(
                name="human",
                min_fraction=0.38,
                max_fraction=0.44,
                description="Human genomic sequence clusters near 41% GC on average.",
            ),
            "virus": GCReferenceRange(
                name="virus",
                min_fraction=0.20,
                max_fraction=0.80,
                description="Viruses show broad GC variability across families.",
            ),
            "generic_fragment": GCReferenceRange(
                name="generic_fragment",
                min_fraction=0.25,
                max_fraction=0.75,
                description="Generic fragments and contigs tolerate wider GC variability.",
            ),
        }

    def get_range(self, key: str) -> GCReferenceRange:
        """Return the reference range for a key."""
        return self._ranges[key]


@dataclass
class GCContextEvaluation:
    """Interpretation of GC content relative to an inferred source."""

    status: str
    explanation: str
    reference_group: str


@dataclass
class CodingPotentialResult:
    """Heuristic coding potential summary for a sequence."""

    score: float
    classification: str
    longest_orf_nt: int
    orf_count: int
    codon_usage_bias: float


@dataclass
class BiologicalSummary:
    """Aggregated biological interpretation for a dataset."""

    sequence_scale: str
    scale_distribution: dict[str, int]
    inferred_organism_type: str
    gc_context: str
    gc_explanation: str
    coding_potential_score: float
    coding_potential_classification: str
    orf_presence_fraction: float
    mean_longest_orf_nt: float
    mean_codon_usage_bias: float
    low_complexity_fraction: float
    homopolymer_fraction: float
    duplicate_fraction: float
    gc_skew_mean_abs: float
    short_sequence_fraction: float
    warnings: list[str] = field(default_factory=list)


def classify_sequence_scale(length: int, config: BiologyConfig) -> str:
    """Classify the biological scale implied by sequence length.

    Args:
        length: Sequence length in bases.
        config: Biology-related thresholds.

    Returns:
        One of ``gene/transcript``, ``contig``, or ``genome``.
    """
    if length < config.gene_max_length:
        return "gene/transcript"
    if length > config.genome_min_length:
        return "genome"
    return "contig"


def infer_organism_type(gc_content: float, length: int, config: BiologyConfig) -> str:
    """Infer a broad organism/context label from GC and length.

    Args:
        gc_content: GC percentage in the range 0-100.
        length: Sequence length in bases.
        config: Biology-related thresholds.

    Returns:
        A coarse organism/context classification label.
    """
    gc_fraction = gc_content / 100.0
    if length < config.gene_max_length:
        return "likely_gene_or_viral"
    if length > config.genome_min_length:
        if 0.38 <= gc_fraction <= 0.44:
            return "likely_human_genome"
        if 0.30 <= gc_fraction <= 0.70:
            return "likely_bacterial_genome"
        return "likely_extreme_gc_genome"
    return "contig_or_fragment"


def evaluate_gc_context(gc_value: float, inferred_type: str) -> GCContextEvaluation:
    """Interpret observed GC content in biological context.

    Args:
        gc_value: Observed GC content in percent.
        inferred_type: Output from ``infer_organism_type``.

    Returns:
        GCContextEvaluation with status and explanation.
    """
    database = GCReferenceDatabase()
    mapping = {
        "likely_gene_or_viral": "virus",
        "likely_human_genome": "human",
        "likely_bacterial_genome": "bacteria",
        "likely_extreme_gc_genome": "generic_fragment",
        "contig_or_fragment": "generic_fragment",
    }
    reference_key = mapping.get(inferred_type, "generic_fragment")
    reference = database.get_range(reference_key)
    gc_fraction = gc_value / 100.0

    if reference.min_fraction <= gc_fraction <= reference.max_fraction:
        status = "normal"
    elif gc_fraction < reference.min_fraction - 0.10 or gc_fraction > reference.max_fraction + 0.10:
        status = "extreme"
    else:
        status = "atypical"

    explanation = (
        f"Observed GC is {gc_value:.2f}%. Relative to the {reference.name} reference "
        f"range of {reference.min_fraction:.2%}-{reference.max_fraction:.2%}, "
        f"this is considered {status}. {reference.description}"
    )
    return GCContextEvaluation(
        status=status,
        explanation=explanation,
        reference_group=reference.name,
    )


def _find_orfs(sequence: str) -> list[int]:
    """Find simple forward-strand ORF lengths.

    Args:
        sequence: DNA sequence.

    Returns:
        List of ORF lengths in nucleotides.
    """
    clean = sequence.upper()
    orfs: list[int] = []
    for frame in range(3):
        for start in range(frame, len(clean) - 2, 3):
            if clean[start : start + 3] != "ATG":
                continue
            for stop in range(start + 3, len(clean) - 2, 3):
                codon = clean[stop : stop + 3]
                if codon in STOP_CODONS:
                    orfs.append(stop + 3 - start)
                    break
    return orfs


def _estimate_codon_usage_bias(sequence: str) -> float:
    """Estimate a simple codon-usage bias score.

    Args:
        sequence: DNA sequence, ideally coding.

    Returns:
        A 0-1 score describing codon-pattern enrichment.
    """
    codons = [
        sequence[index : index + 3]
        for index in range(0, len(sequence) - 2, 3)
        if len(sequence[index : index + 3]) == 3 and "N" not in sequence[index : index + 3]
    ]
    if not codons:
        return 0.0
    counts = Counter(codons)
    dominant_fraction = counts.most_common(1)[0][1] / len(codons)
    diversity_fraction = len(counts) / len(codons)
    return min(1.0, (dominant_fraction * 0.6) + (diversity_fraction * 0.4))


def estimate_coding_potential(
    sequence: str,
    config: BiologyConfig,
) -> CodingPotentialResult:
    """Estimate coding potential using simple ORF and codon heuristics.

    Args:
        sequence: DNA sequence.
        config: Biology-related thresholds.

    Returns:
        CodingPotentialResult with a 0-1 score and coarse classification.
    """
    clean = sequence.upper()
    if not clean:
        return CodingPotentialResult(0.0, "non-coding", 0, 0, 0.0)

    orf_lengths = _find_orfs(clean)
    longest_orf_nt = max(orf_lengths, default=0)
    longest_fraction = longest_orf_nt / len(clean)
    longest_orf_sequence = clean[:longest_orf_nt] if longest_orf_nt else ""
    codon_usage_bias = _estimate_codon_usage_bias(longest_orf_sequence)
    start_stop_bonus = 1.0 if orf_lengths else 0.0
    score = min(
        1.0,
        (0.6 * longest_fraction) + (0.25 * start_stop_bonus) + (0.15 * codon_usage_bias),
    )

    if score >= config.coding_score_coding_threshold:
        classification = "coding"
    elif score <= config.coding_score_noncoding_threshold:
        classification = "non-coding"
    else:
        classification = "uncertain"

    return CodingPotentialResult(
        score=score,
        classification=classification,
        longest_orf_nt=longest_orf_nt,
        orf_count=len(orf_lengths),
        codon_usage_bias=codon_usage_bias,
    )


def _normalized_entropy(sequence: str) -> float:
    """Compute normalized mononucleotide Shannon entropy.

    Args:
        sequence: DNA sequence.

    Returns:
        Entropy normalized to 0-1.
    """
    clean = [base for base in sequence.upper() if base in {"A", "C", "G", "T"}]
    if not clean:
        return 0.0
    counts = Counter(clean)
    total = len(clean)
    entropy = 0.0
    for count in counts.values():
        fraction = count / total
        entropy -= fraction * math.log2(fraction)
    return entropy / 2.0


def _max_homopolymer_run(sequence: str) -> int:
    """Return the longest homopolymer run in a sequence."""
    max_run = 0
    current_run = 0
    previous = ""
    for base in sequence.upper():
        if base == previous:
            current_run += 1
        else:
            previous = base
            current_run = 1
        max_run = max(max_run, current_run)
    return max_run


def _gc_skew(sequence: str) -> float:
    """Compute GC skew for one sequence."""
    clean = sequence.upper()
    g_count = clean.count("G")
    c_count = clean.count("C")
    denominator = g_count + c_count
    if denominator == 0:
        return 0.0
    return (g_count - c_count) / denominator


def analyze_biology(records: list[AnyRecord], config: BiologyConfig) -> BiologicalSummary:
    """Aggregate biological interpretation across records.

    Args:
        records: Parsed FASTA or FASTQ records.
        config: Biology-related thresholds.

    Returns:
        BiologicalSummary for the dataset.
    """
    if not records:
        return BiologicalSummary(
            sequence_scale="unknown",
            scale_distribution={},
            inferred_organism_type="unknown",
            gc_context="atypical",
            gc_explanation="No records available for biological interpretation.",
            coding_potential_score=0.0,
            coding_potential_classification="uncertain",
            orf_presence_fraction=0.0,
            mean_longest_orf_nt=0.0,
            mean_codon_usage_bias=0.0,
            low_complexity_fraction=0.0,
            homopolymer_fraction=0.0,
            duplicate_fraction=0.0,
            gc_skew_mean_abs=0.0,
            short_sequence_fraction=0.0,
            warnings=["No biological analysis was possible because no records were parsed."],
        )

    lengths = [len(record) for record in records]
    scale_distribution = Counter(
        classify_sequence_scale(length, config) for length in lengths
    )
    sequence_scale = classify_sequence_scale(int(sum(lengths) / len(lengths)), config)
    gc_values = [record.gc_content() for record in records]
    mean_gc = sum(gc_values) / len(gc_values)
    inferred_type = infer_organism_type(mean_gc, int(sum(lengths) / len(lengths)), config)
    gc_context = evaluate_gc_context(mean_gc, inferred_type)

    coding_results = [estimate_coding_potential(record.sequence, config) for record in records]
    coding_score = sum(result.score for result in coding_results) / len(coding_results)
    if coding_score >= config.coding_score_coding_threshold:
        coding_classification = "coding"
    elif coding_score <= config.coding_score_noncoding_threshold:
        coding_classification = "non-coding"
    else:
        coding_classification = "uncertain"

    low_complexity_fraction = (
        sum(
            1
            for record in records
            if _normalized_entropy(record.sequence) < config.low_complexity_entropy_threshold
        )
        / len(records)
    )
    homopolymer_fraction = (
        sum(
            1
            for record in records
            if _max_homopolymer_run(record.sequence) >= config.homopolymer_run_threshold
        )
        / len(records)
    )
    duplicate_fraction = 1.0 - (
        len({record.sequence.upper() for record in records}) / len(records)
    )
    gc_skew_mean_abs = (
        sum(abs(_gc_skew(record.sequence)) for record in records) / len(records)
    )
    orf_presence_fraction = (
        sum(1 for result in coding_results if result.longest_orf_nt > 0) / len(coding_results)
    )
    short_sequence_fraction = (
        sum(
            1
            for record in records
            if len(record) < config.short_sequence_length_threshold
        )
        / len(records)
    )
    warnings: list[str] = []
    if low_complexity_fraction > config.low_complexity_fraction_threshold:
        warnings.append("A substantial fraction of records appear low complexity.")
    if duplicate_fraction > config.duplicate_fraction_warning_threshold:
        warnings.append("Sequence duplication is elevated across records.")

    return BiologicalSummary(
        sequence_scale=sequence_scale,
        scale_distribution=dict(scale_distribution),
        inferred_organism_type=inferred_type,
        gc_context=gc_context.status,
        gc_explanation=gc_context.explanation,
        coding_potential_score=coding_score,
        coding_potential_classification=coding_classification,
        orf_presence_fraction=orf_presence_fraction,
        mean_longest_orf_nt=(
            sum(result.longest_orf_nt for result in coding_results) / len(coding_results)
        ),
        mean_codon_usage_bias=(
            sum(result.codon_usage_bias for result in coding_results) / len(coding_results)
        ),
        low_complexity_fraction=low_complexity_fraction,
        homopolymer_fraction=homopolymer_fraction,
        duplicate_fraction=duplicate_fraction,
        gc_skew_mean_abs=gc_skew_mean_abs,
        short_sequence_fraction=short_sequence_fraction,
        warnings=warnings,
    )
