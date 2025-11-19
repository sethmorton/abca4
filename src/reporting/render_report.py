#!/usr/bin/env python3
"""
Render the compact ABCA4 HTML report from a Jinja template.

This focuses on the narrative, overview-style report in
``data_processed/reports/abca4_report.html``. Figures and the manifest
are generated separately by ``src.reporting.plot_figures``.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

from ..config import REPORTS_DIR, logger
from .generate_pdf import generate_compact_report_pdf


CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]
TEMPLATES_DIR = CAMPAIGN_ROOT / "templates"
TEMPLATE_NAME = "abca4_report.html.j2"


def _load_manifest() -> dict:
    """Load manifest produced by plot_figures."""
    manifest_path = REPORTS_DIR / "abca4_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"Manifest not found at {manifest_path}. "
            "Run src.reporting.plot_figures.generate_figures_and_manifest() first."
        )

    with manifest_path.open(encoding="utf-8") as handle:
        manifest = json.load(handle)

    # Minimal sanity checks
    required_keys = ["metadata", "total_variants", "k", "clusters_covered", "total_clusters"]
    missing = [key for key in required_keys if key not in manifest]
    if missing:
        raise ValueError(f"Manifest missing required keys: {missing}")

    return manifest


def render_html_report(manifest: dict) -> Path:
    """Render the HTML report using the Jinja template."""
    env = Environment(
        loader=FileSystemLoader(str(TEMPLATES_DIR)),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = env.get_template(TEMPLATE_NAME)

    html = template.render(manifest=manifest)

    output_path = REPORTS_DIR / "abca4_report.html"
    output_path.write_text(html, encoding="utf-8")
    logger.info("Wrote HTML report to %s", output_path)
    return output_path


def main() -> None:
    """CLI entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
    )
    manifest = _load_manifest()
    render_html_report(manifest)

    # Generate PDF from the HTML report
    try:
        from .generate_pdf import generate_full_report_pdf
        pdf_path = generate_full_report_pdf()
        logger.info("Generated PDF report: %s", pdf_path)
    except Exception as exc:
        logger.error("Failed to generate PDF: %s", exc)
        # Don't fail the entire process if PDF generation fails


if __name__ == "__main__":
    main()

