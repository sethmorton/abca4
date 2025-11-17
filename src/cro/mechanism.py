"""Stage 2: Mechanism Annotator

Tags variants with molecular mechanisms using rule-based system with optional LLM enhancement.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Union

import pandas as pd
import yaml

try:
    # When run as module
    from .cro_types import (
        MechanismAnnotation, MechanismAnnotationData, MechanismPanel, MechanismRuleCondition,
        MechanismRulesConfig, MechanismTag, VariantPanel, VariantPanelRow
    )
except ImportError:
    # When run as standalone script
    import sys
    from pathlib import Path
    # Add the src directory to sys.path
    src_dir = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(src_dir))
    from cro.cro_types import (
        MechanismAnnotation, MechanismPanel, MechanismTag, VariantPanel, VariantPanelRow
    )

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
)
LOGGER = logging.getLogger(__name__)

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]
CRO_DATA_DIR = CAMPAIGN_ROOT / "data_processed" / "cro"


def _load_mechanism_rules(gene: str = "ABCA4") -> MechanismRulesConfig:
    """Load mechanism tagging rules for a specific gene."""
    rules_path = CAMPAIGN_ROOT / "src" / "cro" / "catalog" / f"{gene.lower()}_mechanisms.yaml"
    if not rules_path.exists():
        raise FileNotFoundError(f"Mechanism rules not found: {rules_path}")

    with open(rules_path, "r") as f:
        rules = yaml.safe_load(f)

    return rules  # type: ignore


def _matches_rule_condition(variant: VariantPanelRow, condition: MechanismRuleCondition) -> bool:
    """Check if a variant matches a rule condition."""
    for field, expected in condition.items():
        if not hasattr(variant, field):
            continue

        actual_value = getattr(variant, field)

        # Handle different condition types
        if isinstance(expected, list):
            # Range condition like [0.5, 1.0]
            if len(expected) == 2 and all(isinstance(x, (int, float)) for x in expected):
                min_val, max_val = expected
                if not (min_val <= actual_value <= max_val):
                    return False
            # List membership
            else:
                if actual_value not in expected:
                    return False
        elif isinstance(expected, str):
            # String matching
            if actual_value != expected:
                return False
        else:
            # Exact match
            if actual_value != expected:
                return False

    return True


def _apply_mechanism_rules(
    variant: VariantPanelRow,
    rules: Dict[str, Any]
) -> tuple[Set[MechanismTag], str]:
    """Apply mechanism tagging rules to a variant."""
    mechanisms: Set[MechanismTag] = set()
    rationales: List[str] = []

    # Apply each rule
    for rule in rules.get("rules", []):
        condition = rule.get("condition", {})
        mechanism = rule.get("mechanism")
        rationale = rule.get("rationale", "")

        if _matches_rule_condition(variant, condition):
            if isinstance(mechanism, str):
                mechanisms.add(mechanism)  # type: ignore
                rationales.append(rationale)
            elif isinstance(mechanism, list):
                mechanisms.update(mechanism)  # type: ignore
                rationales.append(rationale)

    # Apply default if no mechanisms matched
    if not mechanisms:
        default_mechanism = rules.get("default_mechanism", "folding_stability")
        default_rationale = rules.get("default_rationale", "Conservative assignment")
        mechanisms.add(default_mechanism)  # type: ignore
        rationales.append(default_rationale)

    combined_rationale = "; ".join(rationales)
    return mechanisms, combined_rationale


def _enhance_with_llm(
    variant: VariantPanelRow,
    mechanisms: Set[MechanismTag],
    rationale: str,
    llm_client: Optional[Any] = None
) -> str:
    """Optionally enhance rationale with LLM insights."""
    if llm_client is None:
        return rationale

    # This would integrate with an LLM service for enhanced rationales
    # For now, return the original rationale
    return rationale


def annotate_mechanisms(
    variant_panel: VariantPanel,
    gene: str = "ABCA4",
    llm_client: Optional[Any] = None
) -> MechanismPanel:
    """Annotate variants with molecular mechanisms.

    Args:
        variant_panel: List of VariantPanelRow objects
        gene: Gene name for rule selection
        llm_client: Optional LLM client for enhanced rationales

    Returns:
        List of MechanismAnnotation objects
    """
    LOGGER.info("Loading mechanism rules for %s", gene)
    rules = _load_mechanism_rules(gene)

    LOGGER.info("Annotating %d variants with mechanisms", len(variant_panel))

    mechanism_panel = []
    for variant in variant_panel:
        mechanisms, rationale = _apply_mechanism_rules(variant, rules)

        # Optionally enhance with LLM
        enhanced_rationale = _enhance_with_llm(variant, mechanisms, rationale, llm_client)

        annotation = MechanismAnnotation(
            variant_id=variant.variant_id,
            mechanism_tags=mechanisms,
            rationale=enhanced_rationale,
        )
        mechanism_panel.append(annotation)

    LOGGER.info("Annotated variants with %d unique mechanism combinations",
                len(set(frozenset(ann.mechanism_tags) for ann in mechanism_panel)))

    return mechanism_panel


def load_variant_panel() -> VariantPanel:
    """Load variant panel from Stage 1 output."""
    parquet_path = CRO_DATA_DIR / "variant_panel.parquet"
    if not parquet_path.exists():
        raise FileNotFoundError(f"Variant panel not found: {parquet_path}. Run Stage 1 first.")

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


def save_mechanism_panel(panel: MechanismPanel, output_dir: Path = CRO_DATA_DIR) -> Path:
    """Save MechanismPanel to parquet and JSON files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Convert to DataFrame
    df_data = []
    for ann in panel:
        df_data.append({
            "variant_id": ann.variant_id,
            "mechanism_tags": list(ann.mechanism_tags),
            "rationale": ann.rationale,
        })

    df = pd.DataFrame(df_data)

    # Save parquet
    parquet_path = output_dir / "mechanism_panel.parquet"
    df.to_parquet(parquet_path)
    LOGGER.info("Saved mechanism panel to %s", parquet_path)

    # Save JSON
    json_path = output_dir / "mechanism_panel.json"
    with open(json_path, "w") as f:
        json.dump([{
            "variant_id": ann.variant_id,
            "mechanism_tags": list(ann.mechanism_tags),
            "rationale": ann.rationale,
        } for ann in panel], f, indent=2)
    LOGGER.info("Saved mechanism panel to %s", json_path)

    return parquet_path


def main() -> None:
    """Run mechanism annotation on current variant panel."""
    variant_panel = load_variant_panel()
    mechanism_panel = annotate_mechanisms(variant_panel)
    save_mechanism_panel(mechanism_panel)

    LOGGER.info("Successfully annotated %d variants with mechanisms", len(mechanism_panel))


if __name__ == "__main__":
    main()
