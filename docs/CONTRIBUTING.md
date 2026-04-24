# Contributing

## Setup

```bash
git clone https://github.com/example/autobiopipe.git
cd autobiopipe
pip install -e .[dev]
```

## Development Guidelines

- Keep public APIs typed and documented.
- Preserve streaming behavior in parsers unless a stage explicitly requires materialized records.
- Add tests for both normal and failure paths.
- Route thresholds through configuration rather than embedding values directly in decision rules.

## Test Workflow

```bash
pytest
```

Coverage is configured with a `fail_under` target of 90%.

## Pull Requests

- Keep changes focused.
- Update docs when CLI, config, or output formats change.
- Include representative tests for new rules, edge cases, or formats.
