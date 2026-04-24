# AutoBioPipe V 0.1.0 

[![Tests](https://github.com/NischayMehta403/autobiopipe/actions/workflows/test.yml/badge.svg)](https://github.com/NischayMehta403/autobiopipe/actions/workflows/test.yml)
[![Python](https://img.shields.io/badge/python-3.9%2B-1f6feb.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-0f766e.svg)](LICENSE)

AutoBioPipe is a Python CLI for sequence-file quality control that goes beyond basic FASTA/FASTQ validation.

It detects file structure, computes core QC metrics, adds lightweight biological interpretation, renders plots, and produces terminal, JSON, CSV, and PDF reports in one run.

AutoBioPipe is aimed at fast first-pass QC for researchers, student projects, prototypes, and command-line workflows where you want something more informative than "file parsed successfully" but lighter than a full downstream pipeline stack.

## Why AutoBioPipe

Most quick QC scripts stop at counts and averages. AutoBioPipe is designed to give you a fuller first-pass read on a dataset:

- format detection for FASTA, FASTQ, and `.gz` inputs
- streaming parsers with explicit error handling
- GC, length, quality, and N-content QC metrics
- lightweight biological heuristics for coding potential, GC context, duplication, low complexity, homopolymers, and GC skew
- generated plots for GC distribution, read lengths, and quality profiles
- a configurable decision engine with rule-based findings and severity scoring
- machine-readable and human-readable reports from the same run

## Pipeline At A Glance

AutoBioPipe runs the following stages:

`detect -> parse -> qc -> biology -> visualization -> decision -> report`

That means a single command can tell you:

- what kind of file you gave it
- whether the records parse cleanly
- whether basic QC metrics look suspicious
- whether the sequences show biologically unusual patterns
- what concrete findings were triggered

## Supported Inputs

AutoBioPipe currently supports:

- FASTA
- FASTQ
- gzipped FASTA and FASTQ
- single-end files
- paired-end filename hints inferred from names such as `R1` and `R2`

Important detail:

- paired-end detection is heuristic and filename-based
- it is helpful metadata, not a guarantee about the experiment design
- AutoBioPipe currently analyzes one input file at a time rather than jointly processing read pairs

## Installation

```bash
git clone https://github.com/NischayMehta403/autobiopipe.git
cd autobiopipe
pip install -e .
```

For development:

```bash
pip install -e .[dev]
```

## Quick Start

Run the full pipeline on the bundled sample FASTQ:

```bash
autobio run examples/sample.fastq --output-dir results
```

Inspect detection only:

```bash
autobio detect examples/sample.fasta
```

Render a saved JSON report back to the terminal:

```bash
autobio report results/sample.fastq_report.json
```

Show the installed version:

```bash
autobio version
```

Run a FASTA example:

```bash
autobio run examples/sample.fasta --output-dir results_fasta
```

Run with more verbose logging:

```bash
autobio run examples/sample.fastq --output-dir results --verbose
```

## What You Get

A typical run writes:

- `*_report.json`
- `*_report.csv`
- `*_report.pdf`
- `gc_distribution.png`
- `length_distribution.png`
- `quality_profile.png`

The JSON report is the canonical structured output. CSV is useful for spreadsheets and downstream tabular review. PDF is intended for quick sharing and archival summaries.

Generated plots include:

- GC distribution
- sequence or read length distribution
- per-position quality profile for FASTQ inputs

## Example Workflow

One practical workflow looks like this:

1. Run `autobio detect` to confirm the file type and compression handling look correct.
2. Run `autobio run` to generate reports and plots.
3. Inspect the terminal summary for immediate findings.
4. Use the JSON output for programmatic downstream use.
5. Share the PDF or CSV output with collaborators if needed.

## CLI

`autobio run <input_file> [--output-dir DIR] [--config FILE] [--verbose]`

- runs the complete pipeline
- returns exit code `0` for `PASS` or `WARNING`
- returns exit code `1` for runtime or input errors
- returns exit code `2` when the pipeline completes but the overall decision is `FAIL`

`autobio detect <input_file>`

- prints detected file type, compression status, confidence, and inferred sequencing layout

`autobio report <json_file>`

- re-renders a previously generated JSON report in the terminal

`autobio version`

- prints the AutoBioPipe version and Python version

## Decision Rules

The decision engine currently includes rules `QC001` through `QC015`.

Covered finding areas include:

- low average quality
- GC anomaly or extreme GC content
- minimum length and abnormal length spread
- high ambiguous-base fraction
- suspiciously small datasets
- low complexity and homopolymer-heavy sequences
- weak coding potential or ORF absence
- suspicious duplication
- GC skew imbalance
- unusually short sequences

Findings are reported with weighted severities:

- `INFO = 1`
- `WARNING = 2`
- `CRITICAL = 3`

This makes the final status more informative than a simple single-threshold pass/fail.

## Configuration

AutoBioPipe uses TOML configuration files.

Examples:

- [examples/config.toml](examples/config.toml)
- [examples/config_strict.toml](examples/config_strict.toml)

Main config sections:

- `[qc]` for core QC thresholds
- `[biology]` for biological heuristics and advanced rule thresholds
- `[visualization]` for plot settings
- `[pipeline]` for runtime and output behavior

Run with a custom config:

```bash
autobio run examples/sample.fastq --config examples/config.toml
```

Representative configuration knobs include:

- `qc.min_avg_quality`
- `qc.gc_content_min`
- `qc.gc_content_max`
- `qc.min_read_length`
- `qc.small_dataset_threshold`
- `biology.low_complexity_fraction_threshold`
- `biology.homopolymer_fraction_threshold`
- `biology.duplicate_fraction_warning_threshold`
- `visualization.dpi`
- `pipeline.output_dir`
- `pipeline.max_records`

Minimal example:

```toml
[qc]
min_avg_quality = 25.0
gc_content_min = 30.0
gc_content_max = 70.0

[pipeline]
output_dir = "results"
max_records = 5000
```

Use the stricter example config if you want a harsher QC posture for demos or experiments:

- [examples/config_strict.toml](examples/config_strict.toml)

## Output Structure

The generated JSON report includes top-level sections for:

- `detection`
- `qc`
- `decision`
- `visualizations`

That makes it suitable for:

- scripting
- notebook analysis
- downstream dashboards
- regression checks in automated workflows

The terminal report is intended for fast review, while JSON is the best format for automation.

## Exit Codes

AutoBioPipe uses exit codes deliberately:

- `0` means the run succeeded and the overall status was `PASS` or `WARNING`
- `1` means there was an input or runtime error
- `2` means the pipeline completed, but the rule engine classified the result as `FAIL`

That makes it easy to integrate into shell scripts and CI pipelines.

## Design Notes

Some implementation choices are intentional:

- parsers are generator-based for memory efficiency
- configuration is TOML-based and human-editable
- decision rules are registry-driven so new checks can be added cleanly
- biological interpretation is heuristic and interpretable rather than black-box ML
- visualization uses a non-interactive matplotlib backend for CLI and CI safety

## Project Structure

- [autobiopipe/](autobiopipe) contains the CLI and pipeline modules
- [tests/](tests) contains parser, QC, decision, pipeline, visualization, biology, and CLI coverage
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) explains the stage design
- [docs/CONTRIBUTING.md](docs/CONTRIBUTING.md) documents contribution guidance

## Verification

The test suite is exercised with `pytest`, and the repository includes GitHub Actions CI in [.github/workflows/test.yml](.github/workflows/test.yml).

Recent local verification:

- `112 passed`
- coverage threshold met at `90.56%`

Run locally:

```bash
pytest
```

## Author

AutoBioPipe is maintained by [NischayMehta403](https://github.com/NischayMehta403).

## License

This project is released under the MIT License.

See [LICENSE](LICENSE) for the full text.

## Notes

- Empty files, malformed FASTQ records, multiline FASTA sequences, gzip compression, and paired-end filename hints are handled explicitly.
- The bundled example files are intentionally small so the project stays lightweight for demos and tests.
- `examples/sample.fastq` is intentionally harsh under the advanced rule set and may return `FAIL`; treat it as a regression fixture, not as a biologically representative "good" sample.
