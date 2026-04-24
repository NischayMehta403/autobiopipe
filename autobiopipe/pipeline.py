"""End-to-end orchestration for the AutoBioPipe QC pipeline."""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from autobiopipe.biology import analyze_biology
from autobiopipe.config import AutoBioPipeConfig, load_config
from autobiopipe.decision import DecisionEngine, DecisionReport
from autobiopipe.detect import DetectionResult, detect_file
from autobiopipe.parser import AnyRecord, parse_file
from autobiopipe.qc import QCResult, run_qc
from autobiopipe.report import (
    print_terminal_report,
    write_csv_report,
    write_json_report,
    write_pdf_report,
)
from autobiopipe.visualization import VisualizationArtifacts, generate_visualizations

logger = logging.getLogger(__name__)


@dataclass
class PipelineResult:
    """Result bundle for a pipeline execution."""

    input_file: Path
    output_dir: Path
    detection: DetectionResult
    qc_result: QCResult
    decision: DecisionReport
    visualizations: VisualizationArtifacts
    output_json: Path
    output_csv: Path
    output_pdf: Path
    stage_timings: dict[str, float] = field(default_factory=dict)
    elapsed_seconds: float = 0.0


def run_pipeline(
    input_file: Path,
    output_dir: Optional[Path] = None,
    config: Optional[AutoBioPipeConfig] = None,
) -> PipelineResult:
    """Run detection, parsing, QC, biology, visualization, decisions, and reporting.

    Args:
        input_file: FASTA or FASTQ file to analyze.
        output_dir: Optional output directory override.
        config: Optional application configuration.

    Returns:
        PipelineResult containing all stage outputs and report paths.

    Raises:
        ValueError: If file type detection fails.
    """
    resolved_config = config or load_config()
    resolved_input = Path(input_file).expanduser().resolve()
    resolved_output_dir = Path(output_dir or resolved_config.pipeline.output_dir)
    resolved_output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Starting pipeline for %s", resolved_input)
    overall_start = time.perf_counter()
    timings: dict[str, float] = {}

    detect_start = time.perf_counter()
    detection = detect_file(resolved_input)
    timings["detect"] = time.perf_counter() - detect_start
    logger.info("Detection completed in %.3fs: %s", timings["detect"], detection.summary())
    if detection.file_type == "unknown":
        raise ValueError(f"Unable to determine file type for '{resolved_input.name}'")

    parse_start = time.perf_counter()
    records: list[AnyRecord] = list(
        parse_file(
            resolved_input,
            detection.file_type,
            max_records=resolved_config.pipeline.max_records,
        )
    )
    timings["parse"] = time.perf_counter() - parse_start
    logger.info("Parsed %d records in %.3fs", len(records), timings["parse"])

    qc_start = time.perf_counter()
    qc_result = run_qc(
        records,
        file_name=resolved_input.name,
        file_type=detection.file_type,
        quality_threshold=resolved_config.qc.min_avg_quality,
    )
    timings["qc"] = time.perf_counter() - qc_start
    logger.info("QC completed in %.3fs", timings["qc"])

    biology_start = time.perf_counter()
    qc_result.biology = analyze_biology(records, resolved_config.biology)
    qc_result.sequence_scale = qc_result.biology.sequence_scale
    qc_result.warnings.extend(qc_result.biology.warnings)
    timings["biology"] = time.perf_counter() - biology_start
    logger.info(
        "Biological interpretation completed in %.3fs (%s, %s)",
        timings["biology"],
        qc_result.sequence_scale,
        qc_result.biology.inferred_organism_type,
    )

    visualization_start = time.perf_counter()
    visualizations = generate_visualizations(
        records,
        qc_result,
        resolved_output_dir,
        resolved_config.visualization,
    )
    timings["visualization"] = time.perf_counter() - visualization_start
    logger.info("Visualization completed in %.3fs", timings["visualization"])

    decision_start = time.perf_counter()
    decision = DecisionEngine(config=resolved_config).evaluate(qc_result)
    timings["decision"] = time.perf_counter() - decision_start
    logger.info(
        "Decision stage completed in %.3fs with status %s",
        timings["decision"],
        decision.overall_status,
    )

    report_start = time.perf_counter()
    output_stem = resolved_input.name
    output_json = resolved_output_dir / f"{output_stem}_report.json"
    output_csv = resolved_output_dir / f"{output_stem}_report.csv"
    output_pdf = resolved_output_dir / f"{output_stem}_report.pdf"
    write_json_report(detection, qc_result, decision, visualizations, output_json)
    write_csv_report(detection, qc_result, decision, visualizations, output_csv)
    write_pdf_report(detection, qc_result, decision, visualizations, output_pdf)
    print_terminal_report(detection, qc_result, decision, visualizations)
    timings["report"] = time.perf_counter() - report_start
    logger.info("Reporting completed in %.3fs", timings["report"])

    elapsed_seconds = time.perf_counter() - overall_start
    logger.info("Pipeline finished in %.3fs", elapsed_seconds)
    return PipelineResult(
        input_file=resolved_input,
        output_dir=resolved_output_dir,
        detection=detection,
        qc_result=qc_result,
        decision=decision,
        visualizations=visualizations,
        output_json=output_json,
        output_csv=output_csv,
        output_pdf=output_pdf,
        stage_timings=timings,
        elapsed_seconds=elapsed_seconds,
    )
