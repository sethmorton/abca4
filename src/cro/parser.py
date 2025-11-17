"""Stage 1: Markdown Report Parser

Parses ABCA4 variant selection reports into structured VariantPanel with proper typing.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Dict, List, Literal

import pandas as pd
import yaml

try:
    # When run as module
    from .cro_types import (
        AssayModuleId, ConsequenceType, DomainType, VariantPanel, VariantPanelRow,
        VariantPanelRowData, VariantReportMetadata
    )
except ImportError:
    # When run as standalone script
    import sys
    from pathlib import Path
    # Add the src directory to sys.path
    src_dir = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(src_dir))
    from cro.cro_types import (
        AssayModuleId, ConsequenceType, DomainType, VariantPanel, VariantPanelRow
    )

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
)
LOGGER = logging.getLogger(__name__)

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]
CRO_DATA_DIR = CAMPAIGN_ROOT / "data_processed" / "cro"
CRO_DATA_DIR.mkdir(parents=True, exist_ok=True)


def _parse_variant_id(variant_str: str) -> Dict[str, str | int]:
    """Parse variant ID like '1:94112975:T/A' into components."""
    match = re.match(r"(\d+):(\d+):([A-Z]+)/([A-Z]+)", variant_str)
    if not match:
        raise ValueError(f"Invalid variant format: {variant_str}")
    chrom, pos, ref, alt = match.groups()
    return {
        "chrom": chrom,
        "pos": int(pos),
        "ref": ref,
        "alt": alt,
        "variant_id": variant_str,
    }


def _normalize_consequence(consequence_str: str) -> ConsequenceType:
    """Normalize consequence string to controlled vocabulary."""
    # Handle compound consequences (take first one for primary type)
    primary = consequence_str.split(",")[0].split("_")[0]

    mapping = {
        "missense": "missense",
        "nonsense": "nonsense",
        "synonymous": "synonymous",
        "splice": "splice_region",
        "intron": "intronic_splice",
        "frameshift": "frameshift",
        "inframe": "inframe_indel",
        "utr": "utr",
        "regulatory": "regulatory",
    }

    for key, value in mapping.items():
        if key in primary.lower():
            return value
    return "other"


def _normalize_domain(cluster_str: str) -> DomainType:
    """Extract domain from cluster string."""
    if "NBD1" in cluster_str:
        return "NBD1"
    elif "NBD2" in cluster_str:
        return "NBD2"
    elif "TMD" in cluster_str:
        return "TMD"
    elif "CTD" in cluster_str:
        return "CTD"
    else:
        return "other"


def _normalize_assay_hint(assay_str: str) -> AssayModuleId | Literal["UNKNOWN"]:
    """Normalize assay hint to controlled vocabulary."""
    assay_lower = assay_str.lower()

    if "differential scanning fluorimetry" in assay_lower or "dsf" in assay_lower:
        if "size exclusion" in assay_lower or "sec" in assay_lower:
            return "DSF_SEC"
        return "DSF_SEC"  # DSF implies SEC in this context

    if "functional" in assay_lower:
        return "FUNCTIONAL"

    if "trafficking" in assay_lower:
        return "TRAFFICKING"

    if "splice minigene" in assay_lower:
        return "SPLICE_MINIGENE"

    if "rna-seq" in assay_lower:
        return "RNA_SEQ_REPORTER"

    if "luciferase" in assay_lower or "reporter" in assay_lower:
        return "TRANSCRIPTIONAL_REPORTER"

    return "UNKNOWN"


def _parse_markdown_table(content: str) -> List[VariantPanelRowData]:
    """Parse markdown table into list of dictionaries."""
    lines = content.split("\n")
    table_start = None

    # Find table start (header row)
    for i, line in enumerate(lines):
        if "| Rank | Variant |" in line:
            table_start = i
            break

    if table_start is None:
        raise ValueError("Could not find variant table in markdown")

    # Parse table rows
    variants = []
    for line in lines[table_start + 2:]:  # Skip header and separator
        line = line.strip()
        if not line or not line.startswith("|"):
            continue
        if line == "|" + "-" * (len(line) - 2) + "|":  # Skip separator lines
            continue

        # Parse table row
        parts = [col.strip() for col in line.split("|")[1:-1]]
        if len(parts) < 8:
            continue

        rank, variant, consequence, impact, model, cluster, coverage, assay = parts[:8]

        try:
            variant_data = _parse_variant_id(variant)
            variants.append({
                "rank": int(rank),
                "variant_id": variant_data["variant_id"],
                "chrom": variant_data["chrom"],
                "pos": variant_data["pos"],
                "ref": variant_data["ref"],
                "alt": variant_data["alt"],
                "consequence_raw": consequence,
                "consequence": _normalize_consequence(consequence),
                "impact_score": float(impact),
                "model_score": float(model),
                "cluster_id": cluster,
                "domain": _normalize_domain(cluster),
                "coverage": float(coverage),
                "assay_hint_raw": assay,
                "assay_hint": _normalize_assay_hint(assay),
            })
        except (ValueError, IndexError) as e:
            LOGGER.warning(f"Failed to parse variant row: {line} - {e}")
            continue

    return variants


def _extract_metadata(content: str) -> VariantReportMetadata:
    """Extract metadata from markdown header."""
    metadata = {
        "gene": "ABCA4",  # Default
        "selection_strategy": "unknown",
        "panel_size": 0,
        "coverage_target": 0.0,
    }

    # Extract strategy
    strategy_match = re.search(r"\*\*Strategy:\*\*\s*(.+)", content)
    if strategy_match:
        metadata["selection_strategy"] = strategy_match.group(1).strip()

    # Extract panel size
    size_match = re.search(r"\*\*Panel Size \(K\):\*\*\s*(\d+)", content)
    if size_match:
        metadata["panel_size"] = int(size_match.group(1))

    # Extract coverage
    coverage_match = re.search(r"\*\*Coverage λ:\*\*\s*([\d.]+)", content)
    if coverage_match:
        metadata["coverage_target"] = float(coverage_match.group(1))

    return metadata


def parse_report_to_variant_panel(
    report_path: Path,
    panel_id: str = "abca4_selected"
) -> VariantPanel:
    """Parse markdown report into typed VariantPanel.

    Args:
        report_path: Path to the markdown report file
        panel_id: Identifier for this variant panel

    Returns:
        List of VariantPanelRow objects
    """
    LOGGER.info("Parsing report from %s", report_path)

    with open(report_path, "r") as f:
        content = f.read()

    # Extract metadata
    metadata = _extract_metadata(content)

    # Parse variants table
    variants_data = _parse_markdown_table(content)

    if not variants_data:
        raise ValueError("No variants found in report")

    LOGGER.info("Parsed %d variants", len(variants_data))

    # Create VariantPanelRow objects
    panel = []
    for row in variants_data:
        panel_row = VariantPanelRow(
            variant_id=row["variant_id"],
            gene=metadata["gene"],
            rank=row["rank"],
            chrom=row["chrom"],
            pos=row["pos"],
            ref=row["ref"],
            alt=row["alt"],
            consequence=row["consequence"],
            domain=row["domain"],
            model_score=row["model_score"],
            impact_score=row["impact_score"],
            cluster_id=row["cluster_id"],
            coverage=row["coverage"],
            assay_hint=row["assay_hint"],
            selection_strategy=metadata["selection_strategy"],
            panel_size=metadata["panel_size"],
            panel_id=panel_id,
            source_md_path=report_path,
        )
        panel.append(panel_row)

    return panel


def save_variant_panel(panel: VariantPanel, output_dir: Path = CRO_DATA_DIR) -> Path:
    """Save VariantPanel to parquet and JSON schema files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Convert to DataFrame
    df_data = []
    for row in panel:
        df_data.append({
            "variant_id": row.variant_id,
            "gene": row.gene,
            "rank": row.rank,
            "chrom": row.chrom,
            "pos": row.pos,
            "ref": row.ref,
            "alt": row.alt,
            "consequence": row.consequence,
            "domain": row.domain,
            "model_score": row.model_score,
            "impact_score": row.impact_score,
            "cluster_id": row.cluster_id,
            "coverage": row.coverage,
            "assay_hint": row.assay_hint,
            "selection_strategy": row.selection_strategy,
            "panel_size": row.panel_size,
            "panel_id": row.panel_id,
            "source_md_path": str(row.source_md_path),
        })

    df = pd.DataFrame(df_data)

    # Save parquet
    parquet_path = output_dir / "variant_panel.parquet"
    df.to_parquet(parquet_path)
    LOGGER.info("Saved variant panel to %s", parquet_path)

    # Generate and save JSON schema
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {
            "variant_id": {"type": "string"},
            "gene": {"type": "string"},
            "rank": {"type": "integer"},
            "chrom": {"type": "string"},
            "pos": {"type": "integer"},
            "ref": {"type": "string"},
            "alt": {"type": "string"},
            "consequence": {
                "type": "string",
                "enum": ["missense", "nonsense", "synonymous", "splice_region",
                        "intronic_splice", "frameshift", "inframe_indel", "utr",
                        "regulatory", "other"]
            },
            "domain": {
                "type": "string",
                "enum": ["NBD1", "NBD2", "TMD", "CTD", "other"]
            },
            "model_score": {"type": "number"},
            "impact_score": {"type": "number"},
            "cluster_id": {"type": "string"},
            "coverage": {"type": "number"},
            "assay_hint": {
                "type": "string",
                "enum": ["DSF_SEC", "FUNCTIONAL", "TRAFFICKING", "SPLICE_MINIGENE",
                        "RNA_SEQ_REPORTER", "TRANSCRIPTIONAL_REPORTER", "UNKNOWN"]
            },
            "selection_strategy": {"type": "string"},
            "panel_size": {"type": "integer"},
            "panel_id": {"type": "string"},
            "source_md_path": {"type": "string"},
        },
        "required": ["variant_id", "gene", "rank", "chrom", "pos", "ref", "alt",
                    "consequence", "domain", "model_score", "impact_score",
                    "cluster_id", "coverage", "assay_hint", "selection_strategy",
                    "panel_size", "panel_id", "source_md_path"]
    }

    schema_path = output_dir / "variant_panel.schema.json"
    with open(schema_path, "w") as f:
        json.dump(schema, f, indent=2)
    LOGGER.info("Saved schema to %s", schema_path)

    return parquet_path


def load_variant_panel() -> VariantPanel:
    """Load variant panel from saved parquet file."""
    parquet_path = CRO_DATA_DIR / "variant_panel.parquet"
    if not parquet_path.exists():
        raise FileNotFoundError(f"Variant panel not found: {parquet_path}. Run parser first.")

    df = pd.read_parquet(parquet_path)

    # Convert back to VariantPanelRow objects
    panel = []
    for _, row in df.iterrows():
        panel_row = VariantPanelRow(
            variant_id=row["variant_id"],
            gene=row["gene"],
            rank=row["rank"],
            chrom=row["chrom"],
            pos=row["pos"],
            ref=row["ref"],
            alt=row["alt"],
            consequence=row["consequence"],
            domain=row["domain"],
            model_score=row["model_score"],
            impact_score=row["impact_score"],
            cluster_id=row["cluster_id"],
            coverage=row["coverage"],
            assay_hint=row["assay_hint"],
            selection_strategy=row["selection_strategy"],
            panel_size=row["panel_size"],
            panel_id=row["panel_id"],
            source_md_path=Path(row["source_md_path"]),
        )
        panel.append(panel_row)

    return panel


def main() -> None:
    """Parse the current report snapshot into variant panel."""
    report_path = CAMPAIGN_ROOT / "data_processed" / "reports" / "report_snapshot.md"
    if not report_path.exists():
        raise FileNotFoundError(f"Report not found: {report_path}")

    panel = parse_report_to_variant_panel(report_path)
    save_variant_panel(panel)

    LOGGER.info("Successfully parsed %d variants into variant panel", len(panel))


if __name__ == "__main__":
    main()
