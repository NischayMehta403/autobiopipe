"""Optional AI-assisted explanation hooks for AutoBioPipe reports."""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Any, Mapping, Optional


@dataclass
class AIExplanationResult:
    """Structured output for optional AI-generated explanations."""

    provider: str
    prompt: str
    explanation: str


def build_explanation_prompt(report_payload: Mapping[str, Any]) -> str:
    """Build a prompt for a generative model using report JSON.

    Args:
        report_payload: Parsed AutoBioPipe JSON report.

    Returns:
        Prompt text suitable for a downstream LLM request.
    """
    return (
        "Explain the following AutoBioPipe QC report for a biologist. "
        "Summarize the main risks, note whether the sample is usable, and keep the answer "
        "grounded in the report values.\n\n"
        f"{json.dumps(report_payload, indent=2)}"
    )


def explain_with_gemini(
    report_payload: Mapping[str, Any],
    model: str = "gemini-1.5-flash",
    api_key: Optional[str] = None,
) -> AIExplanationResult:
    """Placeholder Gemini integration.

    Args:
        report_payload: Parsed AutoBioPipe JSON report.
        model: Gemini model name to use in a future implementation.
        api_key: Optional Gemini API key override.

    Returns:
        AIExplanationResult containing a placeholder explanation.
    """
    resolved_api_key = api_key or os.getenv("GEMINI_API_KEY")
    prompt = build_explanation_prompt(report_payload)
    if not resolved_api_key:
        return AIExplanationResult(
            provider="gemini-placeholder",
            prompt=prompt,
            explanation=(
                "Gemini integration is not configured. Set GEMINI_API_KEY and replace the "
                "placeholder implementation in autobiopipe.ai_explain.explain_with_gemini."
            ),
        )
    return AIExplanationResult(
        provider=f"gemini-placeholder:{model}",
        prompt=prompt,
        explanation=(
            "Gemini API key detected, but live API calls are intentionally not implemented in "
            "this placeholder module yet."
        ),
    )
