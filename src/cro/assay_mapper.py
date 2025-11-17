"""Stage 3: Assay Assignment Mapper

Maps molecular mechanisms to appropriate assay modules with rationales.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Set

import pandas as pd
import yaml

try:
    # When run as module
    from .cro_types import (
        AssayAssignment, AssayAssignmentData, AssayAssignments, AssayCatalog,
        AssayCatalogConfig, AssayModule, AssayModuleId, MechanismAnnotation,
        MechanismPanel, MechanismTag
    )
except ImportError:
    # When run as standalone script
    import sys
    from pathlib import Path
    # Add the src directory to sys.path
    src_dir = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(src_dir))
    from cro.cro_types import (
        AssayAssignment, AssayAssignments, AssayCatalog, AssayModule, AssayModuleId,
        MechanismAnnotation, MechanismPanel, MechanismTag
    )

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
)
LOGGER = logging.getLogger(__name__)

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]
CRO_DATA_DIR = CAMPAIGN_ROOT / "data_processed" / "cro"


def load_assay_catalog() -> AssayCatalog:
    """Load assay module catalog."""
    catalog_path = CAMPAIGN_ROOT / "src" / "cro" / "catalog" / "assay_modules.yaml"
    if not catalog_path.exists():
        raise FileNotFoundError(f"Assay catalog not found: {catalog_path}")

    with open(catalog_path, "r") as f:
        catalog_data: AssayCatalogConfig = yaml.safe_load(f)

    catalog = {}
    for module_id, module_data in catalog_data.items():
        module = AssayModule(
            id=module_id,  # type: ignore
            name=module_data["name"],
            description=module_data["description"],
            inputs=module_data["inputs"],
            outputs=module_data["outputs"],
            estimated_duration_weeks=module_data["estimated_duration_weeks"],
            cost_category=module_data["cost_category"],
            mechanism_mapping=set(module_data["mechanism_mapping"]),  # type: ignore
        )
        catalog[module_id] = module  # type: ignore

    return catalog


def load_mechanism_panel() -> MechanismPanel:
    """Load mechanism panel from Stage 2 output."""
    json_path = CRO_DATA_DIR / "mechanism_panel.json"
    if not json_path.exists():
        raise FileNotFoundError(f"Mechanism panel not found: {json_path}. Run Stage 2 first.")

    with open(json_path, "r") as f:
        data = json.load(f)

    panel = []
    for item in data:
        annotation = MechanismAnnotation(
            variant_id=item["variant_id"],
            mechanism_tags=set(item["mechanism_tags"]),  # type: ignore
            rationale=item["rationale"],
        )
        panel.append(annotation)

    return panel


def _select_assay_modules(
    mechanism_tags: Set[MechanismTag],
    catalog: AssayCatalog
) -> Set[AssayModuleId]:
    """Select appropriate assay modules for given mechanism tags."""
    selected_modules = set()

    for module_id, module in catalog.items():
        # Check if any of the variant's mechanisms match the module's capabilities
        if mechanism_tags & module.mechanism_mapping:
            selected_modules.add(module_id)

    return selected_modules


def _generate_assignment_rationale(
    mechanism_tags: Set[MechanismTag],
    selected_modules: Set[AssayModuleId],
    catalog: AssayCatalog
) -> str:
    """Generate rationale for assay module assignments."""
    if not selected_modules:
        return "No assay modules selected for mechanism tags"

    rationales = []
    for module_id in selected_modules:
        module = catalog[module_id]
        matching_mechanisms = mechanism_tags & module.mechanism_mapping
        rationale = f"{module.name} selected for mechanisms: {', '.join(matching_mechanisms)}"
        rationales.append(rationale)

    return "; ".join(rationales)


def assign_assays(
    mechanism_panel: MechanismPanel,
    catalog: AssayCatalog
) -> AssayAssignments:
    """Assign assay modules to variants based on mechanism annotations.

    Args:
        mechanism_panel: List of MechanismAnnotation objects
        catalog: Assay module catalog

    Returns:
        List of AssayAssignment objects
    """
    LOGGER.info("Assigning assays for %d variants", len(mechanism_panel))

    assignments = []
    for annotation in mechanism_panel:
        assay_modules = _select_assay_modules(annotation.mechanism_tags, catalog)
        rationale = _generate_assignment_rationale(
            annotation.mechanism_tags, assay_modules, catalog
        )

        assignment = AssayAssignment(
            variant_id=annotation.variant_id,
            assay_modules=assay_modules,
            rationale=rationale,
        )
        assignments.append(assignment)

    # Log assignment statistics
    module_counts = {}
    for assignment in assignments:
        for module_id in assignment.assay_modules:
            module_counts[module_id] = module_counts.get(module_id, 0) + 1

    LOGGER.info("Assay module assignments: %s", module_counts)

    # Check for variants with no assignments
    unassigned = [a for a in assignments if not a.assay_modules]
    if unassigned:
        LOGGER.warning("%d variants have no assay assignments", len(unassigned))
        # Create lookup for mechanism annotations
        mech_lookup = {ann.variant_id: ann for ann in mechanism_panel}
        for assignment in unassigned:
            mech_annotation = mech_lookup.get(assignment.variant_id)
            mech_tags = list(mech_annotation.mechanism_tags) if mech_annotation else []
            LOGGER.warning("Unassigned variant: %s (mechanisms: %s)",
                         assignment.variant_id, mech_tags)

    return assignments


def save_assay_assignments(
    assignments: AssayAssignments,
    output_dir: Path = CRO_DATA_DIR
) -> Path:
    """Save AssayAssignments to JSON file."""
    output_dir.mkdir(parents=True, exist_ok=True)

    json_path = output_dir / "assay_assignments.json"
    with open(json_path, "w") as f:
        json.dump([{
            "variant_id": a.variant_id,
            "assay_modules": list(a.assay_modules),
            "rationale": a.rationale,
        } for a in assignments], f, indent=2)

    LOGGER.info("Saved assay assignments to %s", json_path)
    return json_path


def main() -> None:
    """Run assay assignment mapping."""
    catalog = load_assay_catalog()
    mechanism_panel = load_mechanism_panel()
    assignments = assign_assays(mechanism_panel, catalog)
    save_assay_assignments(assignments)

    LOGGER.info("Successfully assigned assays to %d variants", len(assignments))


if __name__ == "__main__":
    main()
