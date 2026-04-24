"""Configuration loading and validation for AutoBioPipe."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Optional

try:
    import tomllib
except ImportError:  # pragma: no cover
    import tomli as tomllib  # type: ignore


def _validate_range(name: str, value: float, low: float, high: float) -> float:
    """Validate that a numeric value falls inside an inclusive range.

    Args:
        name: Human-readable field name.
        value: Numeric value to validate.
        low: Lower inclusive bound.
        high: Upper inclusive bound.

    Returns:
        The validated value.

    Raises:
        ValueError: If the value falls outside the permitted range.
    """
    if not low <= value <= high:
        raise ValueError(f"{name} must be between {low} and {high}; got {value}")
    return value


def _validate_non_negative(name: str, value: float) -> float:
    """Validate that a value is non-negative.

    Args:
        name: Human-readable field name.
        value: Numeric value to validate.

    Returns:
        The validated value.

    Raises:
        ValueError: If the value is negative.
    """
    if value < 0:
        raise ValueError(f"{name} must be non-negative; got {value}")
    return value


@dataclass
class QCConfig:
    """Thresholds that drive QC metrics and rule evaluation."""

    min_avg_quality: float = 20.0
    gc_content_min: float = 35.0
    gc_content_max: float = 65.0
    min_read_length: int = 50
    length_cv_threshold: float = 0.25
    n_fraction_warning_threshold: float = 0.05
    n_fraction_critical_threshold: float = 0.20
    small_dataset_threshold: int = 100
    enable_small_dataset_warning: bool = True

    def __post_init__(self) -> None:
        """Validate QC threshold relationships."""
        _validate_non_negative("min_avg_quality", self.min_avg_quality)
        _validate_range("gc_content_min", self.gc_content_min, 0.0, 100.0)
        _validate_range("gc_content_max", self.gc_content_max, 0.0, 100.0)
        if self.gc_content_min > self.gc_content_max:
            raise ValueError("gc_content_min cannot exceed gc_content_max")
        if self.min_read_length < 1:
            raise ValueError("min_read_length must be at least 1")
        _validate_non_negative("length_cv_threshold", self.length_cv_threshold)
        _validate_range(
            "n_fraction_warning_threshold",
            self.n_fraction_warning_threshold,
            0.0,
            1.0,
        )
        _validate_range(
            "n_fraction_critical_threshold",
            self.n_fraction_critical_threshold,
            0.0,
            1.0,
        )
        if self.n_fraction_warning_threshold > self.n_fraction_critical_threshold:
            raise ValueError(
                "n_fraction_warning_threshold cannot exceed "
                "n_fraction_critical_threshold"
            )
        if self.small_dataset_threshold < 1:
            raise ValueError("small_dataset_threshold must be at least 1")


@dataclass
class BiologyConfig:
    """Thresholds for biological interpretation and advanced QC rules."""

    gene_max_length: int = 10_000
    genome_min_length: int = 1_000_000
    coding_score_coding_threshold: float = 0.65
    coding_score_noncoding_threshold: float = 0.35
    low_complexity_entropy_threshold: float = 0.55
    low_complexity_fraction_threshold: float = 0.30
    homopolymer_run_threshold: int = 8
    homopolymer_fraction_threshold: float = 0.15
    duplicate_fraction_warning_threshold: float = 0.20
    gc_extreme_low: float = 20.0
    gc_extreme_high: float = 80.0
    gc_skew_warning_threshold: float = 0.20
    gc_skew_critical_threshold: float = 0.35
    abnormal_length_ratio_threshold: float = 10.0
    length_std_dev_threshold: float = 100.0
    orf_presence_warning_threshold: float = 0.30
    coding_potential_warning_threshold: float = 0.35
    short_sequence_length_threshold: int = 75
    short_sequence_fraction_threshold: float = 0.20

    def __post_init__(self) -> None:
        """Validate biology thresholds."""
        if self.gene_max_length < 1:
            raise ValueError("gene_max_length must be at least 1")
        if self.genome_min_length <= self.gene_max_length:
            raise ValueError("genome_min_length must exceed gene_max_length")
        _validate_range(
            "coding_score_coding_threshold",
            self.coding_score_coding_threshold,
            0.0,
            1.0,
        )
        _validate_range(
            "coding_score_noncoding_threshold",
            self.coding_score_noncoding_threshold,
            0.0,
            1.0,
        )
        if self.coding_score_noncoding_threshold > self.coding_score_coding_threshold:
            raise ValueError(
                "coding_score_noncoding_threshold cannot exceed "
                "coding_score_coding_threshold"
            )
        _validate_range(
            "low_complexity_entropy_threshold",
            self.low_complexity_entropy_threshold,
            0.0,
            1.0,
        )
        _validate_range(
            "low_complexity_fraction_threshold",
            self.low_complexity_fraction_threshold,
            0.0,
            1.0,
        )
        if self.homopolymer_run_threshold < 2:
            raise ValueError("homopolymer_run_threshold must be at least 2")
        _validate_range(
            "homopolymer_fraction_threshold",
            self.homopolymer_fraction_threshold,
            0.0,
            1.0,
        )
        _validate_range(
            "duplicate_fraction_warning_threshold",
            self.duplicate_fraction_warning_threshold,
            0.0,
            1.0,
        )
        _validate_range("gc_extreme_low", self.gc_extreme_low, 0.0, 100.0)
        _validate_range("gc_extreme_high", self.gc_extreme_high, 0.0, 100.0)
        if self.gc_extreme_low > self.gc_extreme_high:
            raise ValueError("gc_extreme_low cannot exceed gc_extreme_high")
        _validate_range(
            "gc_skew_warning_threshold",
            self.gc_skew_warning_threshold,
            0.0,
            1.0,
        )
        _validate_range(
            "gc_skew_critical_threshold",
            self.gc_skew_critical_threshold,
            0.0,
            1.0,
        )
        if self.gc_skew_warning_threshold > self.gc_skew_critical_threshold:
            raise ValueError(
                "gc_skew_warning_threshold cannot exceed gc_skew_critical_threshold"
            )
        _validate_non_negative(
            "abnormal_length_ratio_threshold",
            self.abnormal_length_ratio_threshold,
        )
        _validate_non_negative("length_std_dev_threshold", self.length_std_dev_threshold)
        _validate_range(
            "orf_presence_warning_threshold",
            self.orf_presence_warning_threshold,
            0.0,
            1.0,
        )
        _validate_range(
            "coding_potential_warning_threshold",
            self.coding_potential_warning_threshold,
            0.0,
            1.0,
        )
        if self.short_sequence_length_threshold < 1:
            raise ValueError("short_sequence_length_threshold must be at least 1")
        _validate_range(
            "short_sequence_fraction_threshold",
            self.short_sequence_fraction_threshold,
            0.0,
            1.0,
        )


@dataclass
class VisualizationConfig:
    """Configuration for generated figures."""

    enabled: bool = True
    dpi: int = 150
    histogram_bins: int = 20

    def __post_init__(self) -> None:
        """Validate visualization options."""
        if self.dpi < 72:
            raise ValueError("dpi must be at least 72")
        if self.histogram_bins < 5:
            raise ValueError("histogram_bins must be at least 5")


@dataclass
class PipelineConfig:
    """Operational settings for pipeline execution."""

    verbose: bool = False
    output_dir: Path = Path("autobiopipe_output")
    max_records: Optional[int] = None

    def __post_init__(self) -> None:
        """Validate pipeline settings."""
        if self.max_records is not None and self.max_records < 1:
            raise ValueError("max_records must be at least 1 when set")
        self.output_dir = Path(self.output_dir)


@dataclass
class AutoBioPipeConfig:
    """Top-level application configuration."""

    qc: QCConfig = field(default_factory=QCConfig)
    biology: BiologyConfig = field(default_factory=BiologyConfig)
    visualization: VisualizationConfig = field(default_factory=VisualizationConfig)
    pipeline: PipelineConfig = field(default_factory=PipelineConfig)


def _expect_mapping(section_name: str, data: Any) -> Mapping[str, Any]:
    """Validate that a TOML section is a mapping.

    Args:
        section_name: TOML section name for error messages.
        data: Parsed TOML value.

    Returns:
        The same value narrowed to a mapping.

    Raises:
        ValueError: If the section is not a mapping.
    """
    if not isinstance(data, Mapping):
        raise ValueError(f"Section '{section_name}' must be a table")
    return data


def _merge_qc_config(values: Mapping[str, Any]) -> QCConfig:
    """Build a validated QCConfig from TOML values.

    Args:
        values: TOML section values.

    Returns:
        QCConfig populated from values and defaults.
    """
    defaults = QCConfig()
    return QCConfig(
        min_avg_quality=float(values.get("min_avg_quality", defaults.min_avg_quality)),
        gc_content_min=float(values.get("gc_content_min", defaults.gc_content_min)),
        gc_content_max=float(values.get("gc_content_max", defaults.gc_content_max)),
        min_read_length=int(values.get("min_read_length", defaults.min_read_length)),
        length_cv_threshold=float(
            values.get("length_cv_threshold", defaults.length_cv_threshold)
        ),
        n_fraction_warning_threshold=float(
            values.get(
                "n_fraction_warning_threshold",
                defaults.n_fraction_warning_threshold,
            )
        ),
        n_fraction_critical_threshold=float(
            values.get(
                "n_fraction_critical_threshold",
                defaults.n_fraction_critical_threshold,
            )
        ),
        small_dataset_threshold=int(
            values.get("small_dataset_threshold", defaults.small_dataset_threshold)
        ),
        enable_small_dataset_warning=bool(
            values.get(
                "enable_small_dataset_warning",
                defaults.enable_small_dataset_warning,
            )
        ),
    )


def _merge_pipeline_config(values: Mapping[str, Any]) -> PipelineConfig:
    """Build a validated PipelineConfig from TOML values.

    Args:
        values: TOML section values.

    Returns:
        PipelineConfig populated from values and defaults.
    """
    defaults = PipelineConfig()
    raw_max_records = values.get("max_records", defaults.max_records)
    return PipelineConfig(
        verbose=bool(values.get("verbose", defaults.verbose)),
        output_dir=Path(values.get("output_dir", defaults.output_dir)),
        max_records=int(raw_max_records) if raw_max_records is not None else None,
    )


def _merge_biology_config(values: Mapping[str, Any]) -> BiologyConfig:
    """Build a validated BiologyConfig from TOML values.

    Args:
        values: TOML section values.

    Returns:
        BiologyConfig populated from values and defaults.
    """
    defaults = BiologyConfig()
    return BiologyConfig(
        gene_max_length=int(values.get("gene_max_length", defaults.gene_max_length)),
        genome_min_length=int(
            values.get("genome_min_length", defaults.genome_min_length)
        ),
        coding_score_coding_threshold=float(
            values.get(
                "coding_score_coding_threshold",
                defaults.coding_score_coding_threshold,
            )
        ),
        coding_score_noncoding_threshold=float(
            values.get(
                "coding_score_noncoding_threshold",
                defaults.coding_score_noncoding_threshold,
            )
        ),
        low_complexity_entropy_threshold=float(
            values.get(
                "low_complexity_entropy_threshold",
                defaults.low_complexity_entropy_threshold,
            )
        ),
        low_complexity_fraction_threshold=float(
            values.get(
                "low_complexity_fraction_threshold",
                defaults.low_complexity_fraction_threshold,
            )
        ),
        homopolymer_run_threshold=int(
            values.get("homopolymer_run_threshold", defaults.homopolymer_run_threshold)
        ),
        homopolymer_fraction_threshold=float(
            values.get(
                "homopolymer_fraction_threshold",
                defaults.homopolymer_fraction_threshold,
            )
        ),
        duplicate_fraction_warning_threshold=float(
            values.get(
                "duplicate_fraction_warning_threshold",
                defaults.duplicate_fraction_warning_threshold,
            )
        ),
        gc_extreme_low=float(values.get("gc_extreme_low", defaults.gc_extreme_low)),
        gc_extreme_high=float(values.get("gc_extreme_high", defaults.gc_extreme_high)),
        gc_skew_warning_threshold=float(
            values.get(
                "gc_skew_warning_threshold",
                defaults.gc_skew_warning_threshold,
            )
        ),
        gc_skew_critical_threshold=float(
            values.get(
                "gc_skew_critical_threshold",
                defaults.gc_skew_critical_threshold,
            )
        ),
        abnormal_length_ratio_threshold=float(
            values.get(
                "abnormal_length_ratio_threshold",
                defaults.abnormal_length_ratio_threshold,
            )
        ),
        length_std_dev_threshold=float(
            values.get("length_std_dev_threshold", defaults.length_std_dev_threshold)
        ),
        orf_presence_warning_threshold=float(
            values.get(
                "orf_presence_warning_threshold",
                defaults.orf_presence_warning_threshold,
            )
        ),
        coding_potential_warning_threshold=float(
            values.get(
                "coding_potential_warning_threshold",
                defaults.coding_potential_warning_threshold,
            )
        ),
        short_sequence_length_threshold=int(
            values.get(
                "short_sequence_length_threshold",
                defaults.short_sequence_length_threshold,
            )
        ),
        short_sequence_fraction_threshold=float(
            values.get(
                "short_sequence_fraction_threshold",
                defaults.short_sequence_fraction_threshold,
            )
        ),
    )


def _merge_visualization_config(values: Mapping[str, Any]) -> VisualizationConfig:
    """Build a validated VisualizationConfig from TOML values.

    Args:
        values: TOML section values.

    Returns:
        VisualizationConfig populated from values and defaults.
    """
    defaults = VisualizationConfig()
    return VisualizationConfig(
        enabled=bool(values.get("enabled", defaults.enabled)),
        dpi=int(values.get("dpi", defaults.dpi)),
        histogram_bins=int(values.get("histogram_bins", defaults.histogram_bins)),
    )


def load_config(config_path: Optional[Path] = None) -> AutoBioPipeConfig:
    """Load configuration from TOML or return defaults.

    Args:
        config_path: Optional path to a TOML configuration file.

    Returns:
        A validated AutoBioPipeConfig instance.

    Raises:
        FileNotFoundError: If the supplied config file does not exist.
        ValueError: If the config cannot be parsed or validated.
    """
    if config_path is None:
        return AutoBioPipeConfig()

    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    try:
        with path.open("rb") as handle:
            data = tomllib.load(handle)
    except tomllib.TOMLDecodeError as exc:
        raise ValueError(f"Invalid TOML in config file '{path}': {exc}") from exc

    if not isinstance(data, Mapping):
        raise ValueError(f"Config file '{path}' must contain a TOML table")

    qc_data = _expect_mapping("qc", data.get("qc", {}))
    biology_data = _expect_mapping("biology", data.get("biology", {}))
    visualization_data = _expect_mapping(
        "visualization",
        data.get("visualization", {}),
    )
    pipeline_data = _expect_mapping("pipeline", data.get("pipeline", {}))
    return AutoBioPipeConfig(
        qc=_merge_qc_config(qc_data),
        biology=_merge_biology_config(biology_data),
        visualization=_merge_visualization_config(visualization_data),
        pipeline=_merge_pipeline_config(pipeline_data),
    )
