# AutoBioPipe

AutoBioPipe is a modular Python CLI for automated quality control of biological sequence files in FASTA and FASTQ format. It now combines core QC metrics with biological interpretation, scientific visualizations, and a configurable multi-rule decision engine to produce researcher-friendly reports.

## Features

- Detects FASTA or FASTQ by content and extension, including `.gz` inputs
- Parses records lazily with custom `ParseError` exceptions and line numbers
- Computes GC content, length statistics, quality-score summaries, and N-base fraction
- Infers sequence scale, organism context, GC context, ORF presence, coding potential, duplication, and GC skew
- Generates `gc_distribution.png`, `length_distribution.png`, and `quality_profile.png`
- Applies configurable QC rules `QC001` to `QC015` with weighted scoring
- Writes terminal, JSON, CSV, and PDF reports with biological interpretation and plot references
- Ships with a Typer CLI, examples, tests, and GitHub Actions

## Installation

```bash
git clone https://github.com/example/autobiopipe.git
cd autobiopipe
pip install -e .
```

For development:

```bash
pip install -e .[dev]
```

## Quick Start

Run the full pipeline on a FASTQ file:

```bash
autobio run examples/sample.fastq --output-dir results
```

Inspect file detection only:

```bash
autobio detect examples/sample.fasta
```

Render a saved JSON report:

```bash
autobio report results/sample.fastq_report.json
```

Show the installed version:

```bash
autobio version
```

## Example Outputs

AutoBioPipe writes:

- `*_report.json`
- `*_report.csv`
- `*_report.pdf`
- `gc_distribution.png`
- `length_distribution.png`
- `quality_profile.png`

The JSON report is the canonical machine-readable output. CSV is useful for spreadsheets, and PDF is intended for sharing or archival summaries.

## Configuration

AutoBioPipe reads TOML configuration files. A default example is provided at [examples/config.toml](/home/Fred/autobiopipe/examples/config.toml) and a stricter variant at [examples/config_strict.toml](/home/Fred/autobiopipe/examples/config_strict.toml).

Main config sections:

- `[qc]`: core QC thresholds
- `[biology]`: biological heuristics and advanced rule thresholds
- `[visualization]`: plot generation settings
- `[pipeline]`: runtime and output controls

Run with a custom config:

```bash
autobio run examples/sample.fastq --config examples/config.toml
```

## CLI Reference

`autobio run <input_file> [--output-dir DIR] [--config FILE] [--verbose]`

- Runs detect -> parse -> qc -> biology -> visualization -> decision -> report
- Exit code `0`: success with `PASS` or `WARNING`
- Exit code `1`: runtime or input error
- Exit code `2`: pipeline completed but overall status is `FAIL`

`autobio detect <input_file>`

- Prints the detected format, compression, and inferred sequencing layout

`autobio report <json_file>`

- Prints a previously generated JSON report in the terminal

`autobio version`

- Prints package and Python version information

## Development

Run tests with coverage:

```bash
pytest
```

Main project docs:

- [docs/CONTRIBUTING.md](/home/Fred/autobiopipe/docs/CONTRIBUTING.md)
- [docs/ARCHITECTURE.md](/home/Fred/autobiopipe/docs/ARCHITECTURE.md)
- [CLAUDE.md](/home/Fred/autobiopipe/CLAUDE.md)

## Notes

- Empty files, malformed FASTQ records, multiline FASTA sequences, gzip compression, and paired-end filename hints are handled explicitly.
- The bundled example reads are intentionally very short to keep tests and demos lightweight.
- The tiny bundled `examples/sample.fastq` is intentionally harsh under the advanced rule set and may return `FAIL`; use it as a regression fixture rather than a "good biological sample" reference.
