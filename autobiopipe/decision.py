"""Rule-based decision engine for AutoBioPipe."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Callable, Optional

from autobiopipe.config import AutoBioPipeConfig
from autobiopipe.qc import QCResult

logger = logging.getLogger(__name__)


class Severity(str, Enum):
    """Severity assigned to a QC finding."""

    INFO = "INFO"
    WARNING = "WARNING"
    CRITICAL = "CRITICAL"

    @property
    def weight(self) -> int:
        """Return the weight associated with a severity level."""
        return {"INFO": 1, "WARNING": 2, "CRITICAL": 3}[self.value]


@dataclass
class Finding:
    """A single rule evaluation result."""

    rule_id: str
    severity: Severity
    metric: str
    message: str
    recommendation: str
    value: Optional[float] = None
    threshold: Optional[float] = None
    possible_cause: str = ""
    next_step: str = ""
    validation_method: str = ""

    def as_dict(self) -> dict[str, object]:
        """Serialize the finding to a JSON-friendly dictionary."""
        return {
            "rule_id": self.rule_id,
            "severity": self.severity.value,
            "metric": self.metric,
            "message": self.message,
            "recommendation": self.recommendation,
            "value": self.value,
            "threshold": self.threshold,
            "possible_cause": self.possible_cause,
            "next_step": self.next_step,
            "validation_method": self.validation_method,
        }


@dataclass
class DecisionReport:
    """Aggregated findings and final status for a QC run."""

    findings: list[Finding] = field(default_factory=list)
    overall_status: str = "PASS"
    summary: str = ""
    total_score: int = 0

    def as_dict(self) -> dict[str, object]:
        """Serialize the report to a JSON-friendly dictionary."""
        return {
            "overall_status": self.overall_status,
            "summary": self.summary,
            "total_score": self.total_score,
            "findings": [finding.as_dict() for finding in self.findings],
        }


RuleFn = Callable[[QCResult, AutoBioPipeConfig], list[Finding]]
RULE_REGISTRY: list[tuple[str, RuleFn]] = []


def register_rule(rule_id: str) -> Callable[[RuleFn], RuleFn]:
    """Register a decision rule function by QC identifier."""

    def decorator(rule: RuleFn) -> RuleFn:
        RULE_REGISTRY.append((rule_id, rule))
        return rule

    return decorator


def _biology_required(qc_result: QCResult) -> bool:
    """Return whether biological annotations are available."""
    return qc_result.biology is not None


@register_rule("QC001")
def qc001_low_quality(qc_result: QCResult, config: AutoBioPipeConfig) -> list[Finding]:
    """QC001: Detect low average FASTQ quality."""
    if qc_result.quality is None:
        return []
    threshold = config.qc.min_avg_quality
    observed = qc_result.quality.mean_per_read
    if observed < threshold:
        return [
            Finding(
                rule_id="QC001",
                severity=Severity.CRITICAL,
                metric="mean_quality",
                message=(
                    f"Average read quality is {observed:.2f}, below the configured minimum "
                    f"of {threshold:.2f}."
                ),
                recommendation="Low quality may indicate sequencing deterioration or adapter contamination. Trim or filter low-quality reads before downstream analysis.",
                value=observed,
                threshold=threshold,
                possible_cause="Sequencing chemistry decay, poor basecalling, or contamination.",
                next_step="Run adapter and quality trimming, then rerun QC.",
                validation_method="Inspect FastQC-style quality profiles or compare before/after trimming.",
            )
        ]
    return []


@register_rule("QC002")
def qc002_gc_content(qc_result: QCResult, config: AutoBioPipeConfig) -> list[Finding]:
    """QC002: Detect GC-content anomalies."""
    gc_mean = qc_result.gc.mean
    if gc_mean < config.qc.gc_content_min:
        return [
            Finding(
                rule_id="QC002",
                severity=Severity.WARNING,
                metric="gc_content",
                message=(
                    f"Mean GC content is {gc_mean:.2f}%, below the configured lower bound "
                    f"of {config.qc.gc_content_min:.2f}%."
                ),
                recommendation="Low GC content may indicate contamination or sequencing bias. Compare against the expected reference organism or assembly target.",
                value=gc_mean,
                threshold=config.qc.gc_content_min,
                possible_cause="AT-rich contamination, amplification bias, or unexpected organism origin.",
                next_step="Compare the sample against an expected reference or taxonomic classification.",
                validation_method="Review alignment or k-mer taxonomy against a known reference set.",
            )
        ]
    if gc_mean > config.qc.gc_content_max:
        return [
            Finding(
                rule_id="QC002",
                severity=Severity.WARNING,
                metric="gc_content",
                message=(
                    f"Mean GC content is {gc_mean:.2f}%, above the configured upper bound "
                    f"of {config.qc.gc_content_max:.2f}%."
                ),
                recommendation="Elevated GC content may reflect contamination, amplification bias, or a high-GC organism. Validate against the expected sample type.",
                value=gc_mean,
                threshold=config.qc.gc_content_max,
                possible_cause="High-GC organism background, PCR bias, or library-prep skew.",
                next_step="Cross-check expected GC distribution and inspect organism identity.",
                validation_method="Compare with a reference genome or taxonomic profiling output.",
            )
        ]
    return []


@register_rule("QC003")
def qc003_length_consistency(qc_result: QCResult, config: AutoBioPipeConfig) -> list[Finding]:
    """QC003: Detect inconsistent or too-short reads."""
    findings: list[Finding] = []
    cv = qc_result.length.coefficient_of_variation
    if cv > config.qc.length_cv_threshold:
        findings.append(
            Finding(
                rule_id="QC003",
                severity=Severity.WARNING,
                metric="length_cv",
                message=(
                    f"Read-length coefficient of variation is {cv:.3f}, above the configured "
                    f"threshold of {config.qc.length_cv_threshold:.3f}."
                ),
                recommendation="Length inconsistency may indicate mixed library types or trimming artifacts. Review the distribution before downstream assembly or mapping.",
                value=cv,
                threshold=config.qc.length_cv_threshold,
                possible_cause="Mixed inserts, trimming artifacts, or corrupted input.",
                next_step="Inspect the length histogram and consider length-based filtering.",
                validation_method="Review the generated length distribution plot.",
            )
        )
    if qc_result.length.min_value < config.qc.min_read_length:
        findings.append(
            Finding(
                rule_id="QC003",
                severity=Severity.WARNING,
                metric="min_read_length",
                message=(
                    f"Shortest record length is {qc_result.length.min_value}, below the configured "
                    f"minimum of {config.qc.min_read_length}."
                ),
                recommendation="Very short sequences may map nonspecifically. Filter short reads prior to downstream analysis.",
                value=float(qc_result.length.min_value),
                threshold=float(config.qc.min_read_length),
                possible_cause="Aggressive trimming, partial records, or fragmented assemblies.",
                next_step="Apply a minimum-length filter and rerun QC.",
                validation_method="Check read-length summaries before and after filtering.",
            )
        )
    return findings


@register_rule("QC004")
def qc004_n_fraction(qc_result: QCResult, config: AutoBioPipeConfig) -> list[Finding]:
    """QC004: Detect elevated ambiguous-base fraction."""
    n_fraction = qc_result.n_fraction
    if n_fraction >= config.qc.n_fraction_critical_threshold:
        return [
            Finding(
                rule_id="QC004",
                severity=Severity.CRITICAL,
                metric="n_fraction",
                message=(
                    f"Ambiguous base fraction is {n_fraction:.2%}, exceeding the configured "
                    f"critical threshold of {config.qc.n_fraction_critical_threshold:.2%}."
                ),
                recommendation="High N content can invalidate downstream interpretation. Investigate sequence generation quality or masked assembly regions.",
                value=n_fraction,
                threshold=config.qc.n_fraction_critical_threshold,
                possible_cause="Low-confidence basecalling, assembly gaps, or poor coverage.",
                next_step="Filter low-confidence sequences or re-evaluate the source data.",
                validation_method="Inspect base-level quality or assembly coverage tracks.",
            )
        ]
    if n_fraction >= config.qc.n_fraction_warning_threshold:
        return [
            Finding(
                rule_id="QC004",
                severity=Severity.WARNING,
                metric="n_fraction",
                message=(
                    f"Ambiguous base fraction is {n_fraction:.2%}, exceeding the configured "
                    f"warning threshold of {config.qc.n_fraction_warning_threshold:.2%}."
                ),
                recommendation="Moderate N content warrants review before analysis-sensitive tasks such as variant calling.",
                value=n_fraction,
                threshold=config.qc.n_fraction_warning_threshold,
                possible_cause="Coverage gaps, low-confidence positions, or masked reference segments.",
                next_step="Assess whether ambiguous positions can be filtered or imputed.",
                validation_method="Inspect ambiguous base locations in a viewer or downstream summary.",
            )
        ]
    return []


@register_rule("QC005")
def qc005_small_dataset(qc_result: QCResult, config: AutoBioPipeConfig) -> list[Finding]:
    """QC005: Warn when the dataset is very small."""
    if not config.qc.enable_small_dataset_warning:
        return []
    if qc_result.total_records < config.qc.small_dataset_threshold:
        return [
            Finding(
                rule_id="QC005",
                severity=Severity.WARNING,
                metric="total_records",
                message=(
                    f"Dataset contains {qc_result.total_records} records, below the configured "
                    f"small-dataset threshold of {config.qc.small_dataset_threshold}."
                ),
                recommendation="Small datasets can give unstable summary statistics. Interpret the run cautiously or gather more data.",
                value=float(qc_result.total_records),
                threshold=float(config.qc.small_dataset_threshold),
                possible_cause="Limited sequencing depth or a test/demo dataset.",
                next_step="If this is real data, confirm the expected sequencing yield.",
                validation_method="Compare the record count to project expectations or sequencing logs.",
            )
        ]
    return []


@register_rule("QC006")
def qc006_extreme_gc_content(qc_result: QCResult, config: AutoBioPipeConfig) -> list[Finding]:
    """QC006: Flag extreme GC content."""
    gc_mean = qc_result.gc.mean
    if gc_mean < config.biology.gc_extreme_low or gc_mean > config.biology.gc_extreme_high:
        threshold = (
            config.biology.gc_extreme_low
            if gc_mean < config.biology.gc_extreme_low
            else config.biology.gc_extreme_high
        )
        return [
            Finding(
                rule_id="QC006",
                severity=Severity.CRITICAL,
                metric="extreme_gc_content",
                message=(
                    f"GC content of {gc_mean:.2f}% falls outside the extreme threshold band "
                    f"({config.biology.gc_extreme_low:.2f}-{config.biology.gc_extreme_high:.2f}%)."
                ),
                recommendation="Extreme GC content can distort amplification and mapping. Confirm organism identity and compare with expected GC content.",
                value=gc_mean,
                threshold=threshold,
                possible_cause="Strong taxonomic mismatch, contamination, or biased enrichment.",
                next_step="Validate taxonomy or reference identity before proceeding.",
                validation_method="Compare GC distribution against a trusted reference or control sample.",
            )
        ]
    return []


@register_rule("QC007")
def qc007_low_complexity_sequence(qc_result: QCResult, config: AutoBioPipeConfig) -> list[Finding]:
    """QC007: Flag low-complexity sequence content."""
    if not _biology_required(qc_result):
        return []
    observed = qc_result.biology.low_complexity_fraction
    if observed > config.biology.low_complexity_fraction_threshold:
        return [
            Finding(
                rule_id="QC007",
                severity=Severity.WARNING,
                metric="low_complexity_fraction",
                message=(
                    f"Low-complexity records represent {observed:.2%} of the dataset, above the "
                    f"configured threshold of {config.biology.low_complexity_fraction_threshold:.2%}."
                ),
                recommendation="Low-complexity sequence can drive spurious alignments. Consider masking, filtering, or validating against expected repeat content.",
                value=observed,
                threshold=config.biology.low_complexity_fraction_threshold,
                possible_cause="Simple repeats, microsatellites, or adapter-derived artifacts.",
                next_step="Inspect repetitive content and mask or filter problematic sequences.",
                validation_method="Review sequence complexity metrics or repeat-content analysis.",
            )
        ]
    return []


@register_rule("QC008")
def qc008_homopolymer_runs(qc_result: QCResult, config: AutoBioPipeConfig) -> list[Finding]:
    """QC008: Flag excessive homopolymer runs."""
    if not _biology_required(qc_result):
        return []
    observed = qc_result.biology.homopolymer_fraction
    if observed > config.biology.homopolymer_fraction_threshold:
        return [
            Finding(
                rule_id="QC008",
                severity=Severity.WARNING,
                metric="homopolymer_fraction",
                message=(
                    f"Homopolymer-rich records represent {observed:.2%} of the dataset, above the "
                    f"configured threshold of {config.biology.homopolymer_fraction_threshold:.2%}."
                ),
                recommendation="Homopolymers can be error-prone for some platforms. Validate whether they are expected biological signal or platform artifacts.",
                value=observed,
                threshold=config.biology.homopolymer_fraction_threshold,
                possible_cause="Platform-specific indel error, poly-A/T tails, or repetitive content.",
                next_step="Inspect representative sequences and trim platform artifacts if necessary.",
                validation_method="Compare affected regions against reference or known transcript structure.",
            )
        ]
    return []


@register_rule("QC009")
def qc009_abnormal_length_distribution(
    qc_result: QCResult,
    config: AutoBioPipeConfig,
) -> list[Finding]:
    """QC009: Flag abnormally wide length ranges."""
    if qc_result.length.min_value == 0 or qc_result.total_records < 2:
        return []
    observed = qc_result.length.max_value / qc_result.length.min_value
    if observed > config.biology.abnormal_length_ratio_threshold:
        return [
            Finding(
                rule_id="QC009",
                severity=Severity.WARNING,
                metric="length_range_ratio",
                message=(
                    f"Length range ratio is {observed:.2f}, above the configured threshold of "
                    f"{config.biology.abnormal_length_ratio_threshold:.2f}."
                ),
                recommendation="A very wide length spread can indicate mixed inputs or partial fragments. Review before assembly or quantification.",
                value=observed,
                threshold=config.biology.abnormal_length_ratio_threshold,
                possible_cause="Mixed library populations, partial transcripts, or fragmented assemblies.",
                next_step="Inspect the length histogram for multimodality or strong tails.",
                validation_method="Review the saved length distribution plot.",
            )
        ]
    return []


@register_rule("QC010")
def qc010_high_variance_in_read_length(
    qc_result: QCResult,
    config: AutoBioPipeConfig,
) -> list[Finding]:
    """QC010: Flag large standard deviation in sequence length."""
    observed = qc_result.length.std_dev
    if observed > config.biology.length_std_dev_threshold:
        return [
            Finding(
                rule_id="QC010",
                severity=Severity.WARNING,
                metric="length_std_dev",
                message=(
                    f"Length standard deviation is {observed:.2f}, above the configured threshold "
                    f"of {config.biology.length_std_dev_threshold:.2f}."
                ),
                recommendation="High length variance can complicate mapping and abundance estimation. Consider stratifying or filtering by size.",
                value=observed,
                threshold=config.biology.length_std_dev_threshold,
                possible_cause="Mixed fragments, inconsistent trimming, or biological heterogeneity.",
                next_step="Check whether the dataset should be partitioned by expected insert size.",
                validation_method="Inspect the spread in the length distribution plot.",
            )
        ]
    return []


@register_rule("QC011")
def qc011_coding_potential_anomaly(
    qc_result: QCResult,
    config: AutoBioPipeConfig,
) -> list[Finding]:
    """QC011: Flag unexpectedly low coding potential."""
    if not _biology_required(qc_result):
        return []
    observed = qc_result.biology.coding_potential_score
    if observed < config.biology.coding_potential_warning_threshold:
        return [
            Finding(
                rule_id="QC011",
                severity=Severity.WARNING,
                metric="coding_potential_score",
                message=(
                    f"Mean coding potential score is {observed:.2f}, below the configured threshold "
                    f"of {config.biology.coding_potential_warning_threshold:.2f}."
                ),
                recommendation="Low coding potential may indicate non-coding sequence, contamination, or fragmented ORFs. Validate against the expected assay type.",
                value=observed,
                threshold=config.biology.coding_potential_warning_threshold,
                possible_cause="Non-coding sequence enrichment, assembly fragmentation, or off-target material.",
                next_step="Compare against expected transcript or coding reference content.",
                validation_method="Run ORF prediction or annotation against a trusted gene set.",
            )
        ]
    return []


@register_rule("QC012")
def qc012_orf_absence_warning(qc_result: QCResult, config: AutoBioPipeConfig) -> list[Finding]:
    """QC012: Flag absence of detectable ORFs."""
    if not _biology_required(qc_result):
        return []
    observed = qc_result.biology.orf_presence_fraction
    if observed < config.biology.orf_presence_warning_threshold:
        return [
            Finding(
                rule_id="QC012",
                severity=Severity.WARNING,
                metric="orf_presence_fraction",
                message=(
                    f"Only {observed:.2%} of records contain detectable ORFs, below the configured "
                    f"threshold of {config.biology.orf_presence_warning_threshold:.2%}."
                ),
                recommendation="Sparse ORF signal may be appropriate for non-coding data, but it is suspicious for expected coding sequence.",
                value=observed,
                threshold=config.biology.orf_presence_warning_threshold,
                possible_cause="Non-coding molecules, assembly fragmentation, or poor sequence completeness.",
                next_step="Confirm whether the experiment targets coding or non-coding molecules.",
                validation_method="Compare against annotation or dedicated ORF-calling tools.",
            )
        ]
    return []


@register_rule("QC013")
def qc013_suspiciously_short_sequences(
    qc_result: QCResult,
    config: AutoBioPipeConfig,
) -> list[Finding]:
    """QC013: Flag a large fraction of very short sequences."""
    if not _biology_required(qc_result):
        return []
    observed = qc_result.biology.short_sequence_fraction
    if observed > config.biology.short_sequence_fraction_threshold:
        return [
            Finding(
                rule_id="QC013",
                severity=Severity.WARNING,
                metric="short_sequence_fraction",
                message=(
                    f"{observed:.2%} of records are shorter than "
                    f"{config.biology.short_sequence_length_threshold} bp, above the configured "
                    f"threshold of {config.biology.short_sequence_fraction_threshold:.2%}."
                ),
                recommendation="An excess of short sequences can reduce mapping specificity. Filter short records if they are not biologically expected.",
                value=observed,
                threshold=config.biology.short_sequence_fraction_threshold,
                possible_cause="Aggressive trimming, degraded input, or fragmented assemblies.",
                next_step="Inspect whether short sequences are expected for the assay type.",
                validation_method="Review length statistics and compare against protocol expectations.",
            )
        ]
    return []


@register_rule("QC014")
def qc014_sequence_duplication_detection(
    qc_result: QCResult,
    config: AutoBioPipeConfig,
) -> list[Finding]:
    """QC014: Flag excessive duplicate sequences."""
    if not _biology_required(qc_result):
        return []
    observed = qc_result.biology.duplicate_fraction
    if observed > config.biology.duplicate_fraction_warning_threshold:
        return [
            Finding(
                rule_id="QC014",
                severity=Severity.WARNING,
                metric="duplicate_fraction",
                message=(
                    f"Duplicate sequence fraction is {observed:.2%}, above the configured threshold "
                    f"of {config.biology.duplicate_fraction_warning_threshold:.2%}."
                ),
                recommendation="High duplication can indicate PCR bottlenecks or overrepresented molecules. Consider deduplication or reviewing library complexity.",
                value=observed,
                threshold=config.biology.duplicate_fraction_warning_threshold,
                possible_cause="PCR overamplification, low library complexity, or small targeted panels.",
                next_step="Quantify whether duplicates are expected for the experiment design.",
                validation_method="Compare duplicate rates before and after deduplication or UMI collapse.",
            )
        ]
    return []


@register_rule("QC015")
def qc015_gc_skew_imbalance(qc_result: QCResult, config: AutoBioPipeConfig) -> list[Finding]:
    """QC015: Flag GC skew imbalance."""
    if not _biology_required(qc_result):
        return []
    observed = qc_result.biology.gc_skew_mean_abs
    if observed >= config.biology.gc_skew_critical_threshold:
        severity = Severity.CRITICAL
        threshold = config.biology.gc_skew_critical_threshold
    elif observed >= config.biology.gc_skew_warning_threshold:
        severity = Severity.WARNING
        threshold = config.biology.gc_skew_warning_threshold
    else:
        return []
    return [
        Finding(
            rule_id="QC015",
            severity=severity,
            metric="gc_skew_mean_abs",
            message=(
                f"Mean absolute GC skew is {observed:.2f}, above the configured threshold "
                f"of {threshold:.2f}."
            ),
            recommendation="Strong GC skew may reflect strand bias or unusual composition. Validate whether this is expected for the organism or assembly context.",
            value=observed,
            threshold=threshold,
            possible_cause="Replication-associated strand bias, compositional imbalance, or assembly issues.",
            next_step="Compare against known GC skew profiles for the expected organism.",
            validation_method="Inspect strand-aware composition plots or reference annotations.",
        )
    ]


class DecisionEngine:
    """Evaluate registered QC rules and derive a final status."""

    def __init__(self, config: AutoBioPipeConfig) -> None:
        """Initialize the engine."""
        self.config = config

    def evaluate(self, qc_result: QCResult) -> DecisionReport:
        """Evaluate all registered rules."""
        findings: list[Finding] = []
        for rule_id, rule in RULE_REGISTRY:
            logger.debug("Running decision rule %s (%s)", rule.__name__, rule_id)
            findings.extend(rule(qc_result, self.config))

        total_score = sum(finding.severity.weight for finding in findings)
        if any(finding.severity == Severity.CRITICAL for finding in findings) or total_score >= 8:
            overall_status = "FAIL"
            summary = f"Critical QC issues were detected (score={total_score})."
        elif any(finding.severity == Severity.WARNING for finding in findings) or total_score >= 3:
            overall_status = "WARNING"
            summary = f"QC completed with warnings that should be reviewed (score={total_score})."
        else:
            overall_status = "PASS"
            summary = f"QC metrics are within configured thresholds (score={total_score})."

        return DecisionReport(
            findings=findings,
            overall_status=overall_status,
            summary=summary,
            total_score=total_score,
        )
