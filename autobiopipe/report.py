"""Terminal and file report generation for AutoBioPipe."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Mapping, Optional

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from autobiopipe._version import __version__
from autobiopipe.decision import DecisionReport, Severity
from autobiopipe.detect import DetectionResult
from autobiopipe.qc import QCResult
from autobiopipe.visualization import VisualizationArtifacts

try:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import mm
    from reportlab.platypus import (
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table as PdfTable,
        TableStyle,
    )

    REPORTLAB_AVAILABLE = True
except ImportError:  # pragma: no cover
    REPORTLAB_AVAILABLE = False


STATUS_COLORS = {"PASS": "green", "WARNING": "yellow", "FAIL": "red"}
SEVERITY_COLORS = {
    Severity.INFO: "cyan",
    Severity.WARNING: "yellow",
    Severity.CRITICAL: "red",
}


def _json_ready(value: Any) -> Any:
    """Recursively convert objects into JSON-serializable values."""
    if is_dataclass(value):
        return _json_ready(asdict(value))
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    return value


def _status_color(status: str) -> str:
    """Return the Rich color for a pipeline status."""
    return STATUS_COLORS.get(status, "white")


def _pdf_escape(text: str) -> str:  # pragma: no cover
    """Escape text for insertion into a basic PDF content stream."""
    return text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def _write_basic_pdf(output_path: Path, lines: list[str]) -> Path:  # pragma: no cover
    """Write a minimal, standards-compliant PDF without external dependencies."""
    content_lines = ["BT", "/F1 11 Tf", "50 790 Td"]
    for index, line in enumerate(lines):
        if index:
            content_lines.append("0 -14 Td")
        content_lines.append(f"({_pdf_escape(line)}) Tj")
    content_lines.append("ET")
    content = "\n".join(content_lines).encode("latin-1", errors="replace")

    objects = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Count 1 /Kids [3 0 R] >>",
        (
            b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 595 842] "
            b"/Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>"
        ),
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
        f"<< /Length {len(content)} >>\nstream\n".encode("ascii") + content + b"\nendstream",
    ]

    pdf = bytearray(b"%PDF-1.4\n")
    offsets = [0]
    for index, obj in enumerate(objects, start=1):
        offsets.append(len(pdf))
        pdf.extend(f"{index} 0 obj\n".encode("ascii"))
        pdf.extend(obj)
        pdf.extend(b"\nendobj\n")

    xref_offset = len(pdf)
    pdf.extend(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
    pdf.extend(b"0000000000 65535 f \n")
    for offset in offsets[1:]:
        pdf.extend(f"{offset:010d} 00000 n \n".encode("ascii"))
    pdf.extend(
        (
            f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\n"
            f"startxref\n{xref_offset}\n%%EOF\n"
        ).encode("ascii")
    )
    output_path.write_bytes(pdf)
    return output_path


def print_terminal_report(
    detection: DetectionResult,
    qc_result: QCResult,
    decision: DecisionReport,
    visualizations: VisualizationArtifacts,
    console: Optional[Console] = None,
) -> None:
    """Render a Rich terminal report."""
    rich_console = console or Console()
    status_color = _status_color(decision.overall_status)

    rich_console.print(
        Panel(
            f"[bold]File:[/bold] {detection.file_path.name}\n"
            f"[bold]Status:[/bold] [{status_color}]{decision.overall_status}[/{status_color}]\n"
            f"[bold]Summary:[/bold] {decision.summary}",
            title="AutoBioPipe QC Report",
            border_style=status_color,
        )
    )

    detection_table = Table(title="Detection")
    detection_table.add_column("Field", style="bold cyan")
    detection_table.add_column("Value")
    detection_table.add_row("File type", detection.file_type)
    detection_table.add_row("Compressed", "yes" if detection.is_compressed else "no")
    detection_table.add_row("Sequencing", detection.sequencing_type)
    detection_table.add_row("Confidence", detection.confidence)
    detection_table.add_row("Method", detection.detection_method)
    rich_console.print(detection_table)

    qc_table = Table(title="QC Metrics")
    qc_table.add_column("Metric", style="bold cyan")
    qc_table.add_column("Value")
    qc_table.add_row("Total records", str(qc_result.total_records))
    qc_table.add_row("GC mean", f"{qc_result.gc.mean:.2f}%")
    qc_table.add_row("Length mean", f"{qc_result.length.mean:.2f}")
    qc_table.add_row("Length CV", f"{qc_result.length.coefficient_of_variation:.3f}")
    qc_table.add_row("N fraction", f"{qc_result.n_fraction:.2%}")
    qc_table.add_row("Sequence scale", qc_result.sequence_scale)
    if qc_result.quality is not None:
        qc_table.add_row("Mean quality", f"{qc_result.quality.mean_per_read:.2f}")
        qc_table.add_row(
            "Reads below threshold",
            (
                f"{qc_result.quality.reads_below_threshold} "
                f"({qc_result.quality.reads_below_threshold_pct:.1f}%)"
            ),
        )
    rich_console.print(qc_table)

    if qc_result.biology is not None:
        biology_table = Table(title="Biological Interpretation")
        biology_table.add_column("Field", style="bold cyan")
        biology_table.add_column("Value")
        biology_table.add_row("Inferred type", qc_result.biology.inferred_organism_type)
        biology_table.add_row("GC context", qc_result.biology.gc_context)
        biology_table.add_row(
            "Coding potential",
            (
                f"{qc_result.biology.coding_potential_classification} "
                f"({qc_result.biology.coding_potential_score:.2f})"
            ),
        )
        biology_table.add_row(
            "ORF presence",
            f"{qc_result.biology.orf_presence_fraction:.2%}",
        )
        biology_table.add_row(
            "Duplicate fraction",
            f"{qc_result.biology.duplicate_fraction:.2%}",
        )
        biology_table.add_row(
            "GC skew mean abs",
            f"{qc_result.biology.gc_skew_mean_abs:.2f}",
        )
        rich_console.print(biology_table)
        rich_console.print(Panel(qc_result.biology.gc_explanation, title="GC Context"))

    findings_table = Table(title="Findings")
    findings_table.add_column("Rule", style="bold")
    findings_table.add_column("Severity")
    findings_table.add_column("Metric")
    findings_table.add_column("Message")
    findings_table.add_column("Next step")
    if decision.findings:
        for finding in decision.findings:
            findings_table.add_row(
                finding.rule_id,
                f"[{SEVERITY_COLORS[finding.severity]}]{finding.severity.value}[/{SEVERITY_COLORS[finding.severity]}]",
                finding.metric,
                finding.message,
                finding.next_step or finding.recommendation,
            )
    else:
        findings_table.add_row("-", "INFO", "-", "No rule-based findings.", "No action required.")
    rich_console.print(findings_table)

    plot_table = Table(title="Generated Plots")
    plot_table.add_column("Plot", style="bold cyan")
    plot_table.add_column("Path")
    plot_table.add_row("GC distribution", str(visualizations.gc_distribution))
    plot_table.add_row("Length distribution", str(visualizations.length_distribution))
    plot_table.add_row("Quality profile", str(visualizations.quality_profile))
    rich_console.print(plot_table)


def _build_report_payload(
    detection: DetectionResult,
    qc_result: QCResult,
    decision: DecisionReport,
    visualizations: VisualizationArtifacts,
) -> dict[str, Any]:
    """Build a structured payload used by all report writers."""
    return {
        "autobiopipe_version": __version__,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "detection": _json_ready(detection),
        "qc": _json_ready(qc_result),
        "decision": _json_ready(decision.as_dict()),
        "visualizations": _json_ready(visualizations),
    }


def write_json_report(
    detection: DetectionResult,
    qc_result: QCResult,
    decision: DecisionReport,
    visualizations: VisualizationArtifacts,
    output_path: Path,
) -> Path:
    """Write a JSON report to disk."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(
            _build_report_payload(detection, qc_result, decision, visualizations),
            handle,
            indent=2,
        )
    return output_path


def write_csv_report(
    detection: DetectionResult,
    qc_result: QCResult,
    decision: DecisionReport,
    visualizations: VisualizationArtifacts,
    output_path: Path,
) -> Path:
    """Write a CSV report with metric, biology, visualization, and finding sections."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["section", "field", "value"])
        writer.writerow(["detection", "file_type", detection.file_type])
        writer.writerow(["detection", "is_compressed", detection.is_compressed])
        writer.writerow(["detection", "sequencing_type", detection.sequencing_type])
        writer.writerow(["detection", "confidence", detection.confidence])
        writer.writerow(["qc", "total_records", qc_result.total_records])
        writer.writerow(["qc", "gc_mean", f"{qc_result.gc.mean:.4f}"])
        writer.writerow(["qc", "length_mean", f"{qc_result.length.mean:.4f}"])
        writer.writerow(["qc", "length_cv", f"{qc_result.length.coefficient_of_variation:.4f}"])
        writer.writerow(["qc", "n_fraction", f"{qc_result.n_fraction:.6f}"])
        writer.writerow(["qc", "sequence_scale", qc_result.sequence_scale])
        if qc_result.quality is not None:
            writer.writerow(["qc", "mean_quality", f"{qc_result.quality.mean_per_read:.4f}"])
            writer.writerow(
                ["qc", "reads_below_threshold_pct", f"{qc_result.quality.reads_below_threshold_pct:.4f}"]
            )
        if qc_result.biology is not None:
            writer.writerow(["biology", "inferred_organism_type", qc_result.biology.inferred_organism_type])
            writer.writerow(["biology", "gc_context", qc_result.biology.gc_context])
            writer.writerow(["biology", "coding_potential_score", f"{qc_result.biology.coding_potential_score:.4f}"])
            writer.writerow(["biology", "coding_potential_classification", qc_result.biology.coding_potential_classification])
            writer.writerow(["biology", "duplicate_fraction", f"{qc_result.biology.duplicate_fraction:.4f}"])
            writer.writerow(["biology", "gc_skew_mean_abs", f"{qc_result.biology.gc_skew_mean_abs:.4f}"])
        writer.writerow(["visualization", "gc_distribution", str(visualizations.gc_distribution)])
        writer.writerow(["visualization", "length_distribution", str(visualizations.length_distribution)])
        writer.writerow(["visualization", "quality_profile", str(visualizations.quality_profile)])
        writer.writerow([])
        writer.writerow(
            [
                "findings",
                "rule_id",
                "severity",
                "metric",
                "message",
                "recommendation",
                "possible_cause",
                "next_step",
                "validation_method",
                "value",
                "threshold",
            ]
        )
        for finding in decision.findings:
            writer.writerow(
                [
                    "finding",
                    finding.rule_id,
                    finding.severity.value,
                    finding.metric,
                    finding.message,
                    finding.recommendation,
                    finding.possible_cause,
                    finding.next_step,
                    finding.validation_method,
                    finding.value,
                    finding.threshold,
                ]
            )

    return output_path


def _pdf_table(rows: list[list[str]]) -> PdfTable:
    """Create a consistently styled PDF table."""
    table = PdfTable(rows)
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1F4E79")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#C7D5E0")),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F4F8FB")]),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    return table


def write_pdf_report(
    detection: DetectionResult,
    qc_result: QCResult,
    decision: DecisionReport,
    visualizations: VisualizationArtifacts,
    output_path: Path,
) -> Path:
    """Write a PDF summary report."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not REPORTLAB_AVAILABLE:  # pragma: no cover
        lines = [
            "AutoBioPipe QC Report",
            f"File: {detection.file_path.name}",
            f"Status: {decision.overall_status}",
            f"Summary: {decision.summary}",
            f"File type: {detection.file_type}",
            f"Sequencing: {detection.sequencing_type}",
            f"GC mean: {qc_result.gc.mean:.2f}%",
            f"Length mean: {qc_result.length.mean:.2f}",
            f"N fraction: {qc_result.n_fraction:.2%}",
            f"Sequence scale: {qc_result.sequence_scale}",
        ]
        if qc_result.biology is not None:
            lines.append(f"Inferred type: {qc_result.biology.inferred_organism_type}")
            lines.append(f"GC context: {qc_result.biology.gc_context}")
            lines.append(f"Coding potential: {qc_result.biology.coding_potential_score:.2f}")
        lines.append(f"GC plot: {visualizations.gc_distribution}")
        lines.append(f"Length plot: {visualizations.length_distribution}")
        lines.append(f"Quality plot: {visualizations.quality_profile}")
        for finding in decision.findings:
            lines.append(f"{finding.rule_id} [{finding.severity.value}] {finding.message}")
        return _write_basic_pdf(output_path, lines)

    document = SimpleDocTemplate(
        str(output_path),
        pagesize=A4,
        leftMargin=18 * mm,
        rightMargin=18 * mm,
        topMargin=18 * mm,
        bottomMargin=18 * mm,
    )
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        "AutoBioPipeTitle",
        parent=styles["Heading1"],
        textColor=colors.HexColor("#1F4E79"),
    )
    elements: list[Any] = [
        Paragraph("AutoBioPipe QC Report", title_style),
        Paragraph(f"File: {detection.file_path.name}", styles["Normal"]),
        Paragraph(f"Status: {decision.overall_status}", styles["Normal"]),
        Paragraph(f"Summary: {decision.summary}", styles["Normal"]),
        Spacer(1, 6),
        _pdf_table(
            [
                ["Detection field", "Value"],
                ["File type", detection.file_type],
                ["Compressed", "yes" if detection.is_compressed else "no"],
                ["Sequencing", detection.sequencing_type],
                ["Confidence", detection.confidence],
            ]
        ),
        Spacer(1, 8),
        _pdf_table(
            [
                ["QC metric", "Value"],
                ["Total records", str(qc_result.total_records)],
                ["GC mean", f"{qc_result.gc.mean:.2f}%"],
                ["Length mean", f"{qc_result.length.mean:.2f}"],
                ["Length CV", f"{qc_result.length.coefficient_of_variation:.3f}"],
                ["N fraction", f"{qc_result.n_fraction:.2%}"],
                ["Sequence scale", qc_result.sequence_scale],
                [
                    "Mean quality",
                    f"{qc_result.quality.mean_per_read:.2f}" if qc_result.quality else "n/a",
                ],
            ]
        ),
        Spacer(1, 8),
    ]
    if qc_result.biology is not None:
        elements.extend(
            [
                _pdf_table(
                    [
                        ["Biological field", "Value"],
                        ["Inferred type", qc_result.biology.inferred_organism_type],
                        ["GC context", qc_result.biology.gc_context],
                        ["Coding potential", f"{qc_result.biology.coding_potential_score:.2f}"],
                        ["ORF presence", f"{qc_result.biology.orf_presence_fraction:.2%}"],
                        ["Duplicate fraction", f"{qc_result.biology.duplicate_fraction:.2%}"],
                        ["GC skew mean abs", f"{qc_result.biology.gc_skew_mean_abs:.2f}"],
                    ]
                ),
                Spacer(1, 4),
                Paragraph(qc_result.biology.gc_explanation, styles["BodyText"]),
                Spacer(1, 8),
            ]
        )

    finding_rows = [["Rule", "Severity", "Metric", "Message"]]
    for finding in decision.findings:
        finding_rows.append(
            [
                finding.rule_id,
                finding.severity.value,
                finding.metric,
                finding.message,
            ]
        )
    if len(finding_rows) == 1:
        finding_rows.append(["-", "INFO", "-", "No rule-based findings."])
    elements.extend(
        [
            _pdf_table(finding_rows),
            Spacer(1, 8),
            _pdf_table(
                [
                    ["Plot", "Path"],
                    ["GC distribution", str(visualizations.gc_distribution)],
                    ["Length distribution", str(visualizations.length_distribution)],
                    ["Quality profile", str(visualizations.quality_profile)],
                ]
            ),
        ]
    )
    document.build(elements)
    return output_path


def print_saved_report(json_file: Path, console: Optional[Console] = None) -> None:
    """Print a previously generated JSON report to the terminal."""
    rich_console = console or Console()
    with Path(json_file).open("r", encoding="utf-8") as handle:
        payload: Mapping[str, Any] = json.load(handle)

    status = str(payload["decision"]["overall_status"])
    panel_color = _status_color(status)
    rich_console.print(
        Panel(
            f"[bold]File:[/bold] {payload['detection']['file_path']}\n"
            f"[bold]Status:[/bold] [{panel_color}]{status}[/{panel_color}]\n"
            f"[bold]Generated:[/bold] {payload['generated_at']}",
            title="Saved AutoBioPipe Report",
            border_style=panel_color,
        )
    )

    if "visualizations" in payload:
        rich_console.print(
            f"Plots: {payload['visualizations']['gc_distribution']}, "
            f"{payload['visualizations']['length_distribution']}, "
            f"{payload['visualizations']['quality_profile']}"
        )

    findings_table = Table(title="Findings")
    findings_table.add_column("Rule")
    findings_table.add_column("Severity")
    findings_table.add_column("Metric")
    findings_table.add_column("Message")
    for finding in payload["decision"]["findings"]:
        findings_table.add_row(
            str(finding["rule_id"]),
            str(finding["severity"]),
            str(finding["metric"]),
            str(finding["message"]),
        )
    rich_console.print(findings_table)
