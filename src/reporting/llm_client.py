"""
LLM Client for Assay Draft Generation

Handles Groq API communication with retries, timeouts, and provenance logging.
"""

import hashlib
import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Any, Optional

from groq import Groq

from ..config import (
    GROQ_API_KEY, LLM_MODEL, LLM_TEMP, LLM_MAX_TOKENS,
    logger, validate_llm_config, RUN_ID
)


@dataclass
class LLMConfig:
    """Configuration for LLM calls."""
    model: str
    temperature: float
    max_tokens: int
    api_key: str


def get_prompt_hash(prompt: str) -> str:
    """Generate SHA256 hash of prompt for provenance tracking."""
    return hashlib.sha256(prompt.encode()).hexdigest()[:16]


def create_llm_config() -> LLMConfig:
    """Create validated LLM configuration."""
    validate_llm_config()
    return LLMConfig(
        model=LLM_MODEL,
        temperature=LLM_TEMP,
        max_tokens=LLM_MAX_TOKENS,
        api_key=GROQ_API_KEY
    )


class AssayDraftError(Exception):
    """Custom exception for assay draft generation errors."""
    pass


def generate_assay_markdown(variant: Dict[str, Any], cfg: LLMConfig) -> tuple[str, str, int]:
    """
    Generate assay protocol markdown using Groq LLM.

    Args:
        variant: Dictionary with variant information
        cfg: LLM configuration

    Returns:
        Tuple of (generated_markdown, prompt_hash, tokens_used)

    Raises:
        AssayDraftError: On API failures or validation errors
    """
    from jinja2 import Environment, FileSystemLoader
    from ..config import CAMPAIGN_ROOT

    # Load and render prompt template
    template_dir = CAMPAIGN_ROOT / "src" / "reporting" / "templates"
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template("assay_prompt.md.jinja")
    prompt = template.render(**variant)

    # Generate prompt hash for provenance
    prompt_hash = get_prompt_hash(prompt)

    # Initialize Groq client
    client = Groq(api_key=cfg.api_key)

    start_time = time.time()
    max_retries = 3
    backoff_factor = 2

    for attempt in range(max_retries):
        try:
            logger.info(f"Calling Groq API for variant {variant['variant_id']} (attempt {attempt + 1})")

            response = client.chat.completions.create(
                model=cfg.model,
                messages=[
                    {
                        "role": "system",
                        "content": "You are an expert in designing functional assays for genetic variants. "
                                 "Generate concise, feasible assay protocols following the exact format specified."
                    },
                    {"role": "user", "content": prompt}
                ],
                temperature=cfg.temperature,
                max_tokens=cfg.max_tokens,
                timeout=30.0  # 30 second timeout
            )

            # Extract response
            if not response.choices:
                raise AssayDraftError("No response choices returned from LLM")

            content = response.choices[0].message.content
            if not content:
                raise AssayDraftError("Empty response content from LLM")

            # Validate word count (<300 words)
            word_count = len(content.split())
            if word_count > 300:
                logger.warning(f"LLM response exceeded word limit: {word_count} words")
                # Truncate to ~250 words to stay safe
                words = content.split()[:250]
                content = " ".join(words) + "... [truncated for length]"

            # Validate response format
            if not validate_assay_response(content):
                raise AssayDraftError("LLM response failed validation checks - invalid assay format")

            # Log provenance
            duration = time.time() - start_time
            tokens_used = response.usage.total_tokens if response.usage else 0

            logger.info(
                f"LLM call successful: model={cfg.model}, "
                f"temp={cfg.temperature}, tokens={tokens_used}, "
                f"duration={duration:.2f}s, prompt_hash={prompt_hash}"
            )

            return content, prompt_hash, tokens_used

        except Exception as e:
            if attempt == max_retries - 1:  # Last attempt
                error_msg = f"LLM API call failed after {max_retries} attempts: {str(e)}"
                logger.error(error_msg)
                raise AssayDraftError(error_msg) from e

            # Exponential backoff
            wait_time = backoff_factor ** attempt
            logger.warning(f"LLM API call failed (attempt {attempt + 1}), retrying in {wait_time}s: {str(e)}")
            time.sleep(wait_time)


def validate_assay_response(response: str) -> bool:
    """
    Validate that LLM response contains required assay components.

    Args:
        response: Raw LLM response text

    Returns:
        True if response appears valid
    """
    required_patterns = [
        "assay type:",
        "cell line:",
        "construct design:",
        "readout:",
        "controls:",
        "expected wt vs mutant:",
        "effort:",
        "budget note:",
        "protocol:"
    ]

    response_lower = response.lower()
    missing_patterns = []

    for pattern in required_patterns:
        if pattern not in response_lower:
            missing_patterns.append(pattern)

    if missing_patterns:
        logger.warning(f"LLM response missing required patterns: {missing_patterns}")
        return False

    # Check for allowed assay types (more flexible matching)
    allowed_assays = ["minigene", "trafficking", "wb", "thermal shift", "activity"]
    has_allowed_assay = any(assay in response_lower for assay in allowed_assays)

    if not has_allowed_assay:
        logger.warning("LLM response does not contain allowed assay type")
        return False

    # Check word count is reasonable (not too short)
    word_count = len(response.split())
    if word_count < 50:
        logger.warning(f"LLM response too short: {word_count} words")
        return False

    return True
