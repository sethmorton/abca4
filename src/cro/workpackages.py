"""Stage 4: Work Package Aggregator

Groups variants into work packages by gene × assay_module combinations.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Set

import pandas as pd

try:
    # When run as module
    from .assay_mapper import load_assay_catalog
    from .cro_types import (
        AssayAssignment, AssayAssignments, AssayCatalog, AssayModuleId,
        VariantPanel, WorkPackage, WorkPackages
    )
    from .parser import load_variant_panel
except ImportError:
    # When run as standalone script
    import sys
    from pathlib import Path
    # Add the src directory to sys.path
    src_dir = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(src_dir))
    from cro.assay_mapper import load_assay_catalog
    from cro.cro_types import (
        AssayAssignment, AssayAssignments, AssayCatalog, AssayModuleId,
        VariantPanel, WorkPackage, WorkPackages
    )
    from cro.parser import load_variant_panel

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
)
LOGGER = logging.getLogger(__name__)

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]
CRO_DATA_DIR = CAMPAIGN_ROOT / "data_processed" / "cro"


def load_assay_assignments() -> AssayAssignments:
    """Load assay assignments from Stage 3 output."""
    json_path = CRO_DATA_DIR / "assay_assignments.json"
    if not json_path.exists():
        raise FileNotFoundError(f"Assay assignments not found: {json_path}. Run Stage 3 first.")

    with open(json_path, "r") as f:
        data = json.load(f)

    assignments = []
    for item in data:
        assignment = AssayAssignment(
            variant_id=item["variant_id"],
            assay_modules=set(item["assay_modules"]),  # type: ignore
            rationale=item["rationale"],
        )
        assignments.append(assignment)

    return assignments


def _generate_work_package_id(
    gene: str,
    assay_module: AssayModuleId,
    variant_count: int
) -> str:
    """Generate a unique work package ID."""
    return f"{gene}_{assay_module}_{variant_count:02d}variants"


def _get_work_package_objective(
    gene: str,
    assay_module: AssayModuleId,
    catalog: AssayCatalog,
    variant_count: int
) -> str:
    """Generate objective description for work package."""
    module = catalog[assay_module]
    return f"Evaluate {variant_count} {gene} variants using {module.name} to assess {', '.join(module.outputs)}"


def _get_materials_provided(
    gene: str,
    assay_module: AssayModuleId,
    catalog: AssayCatalog,
    variant_count: int
) -> List[str]:
    """Get list of materials that will be provided by the customer."""
    materials = [
        f"DNA constructs for {variant_count} variants",
        f"WT {gene} control construct",
    ]

    # Note: Instruments, reagents, and specialized materials are provided by CRO
    # Customer provides only the variant DNA constructs and controls

    return materials


def _get_materials_needed(
    assay_module: AssayModuleId,
    catalog: AssayCatalog
) -> List[str]:
    """Get list of materials that CRO needs to provide."""
    materials = [
        "Assay buffers and reagents",
        "Standard laboratory equipment",
    ]

    # Add module-specific instrumentation and reagents provided by CRO
    if assay_module == "DSF_SEC":
        materials.extend([
            "Thermal cycler for DSF",
            "Size exclusion chromatography system",
            "Cell lines or biochemical reagents as specified",
        ])
    elif assay_module == "FUNCTIONAL":
        materials.extend([
            "Functional assay detection system",
            "Cell lines or biochemical reagents as specified",
        ])
    elif assay_module == "TRAFFICKING":
        materials.extend([
            "Microscopy equipment for localization studies",
            "Cell lines or biochemical reagents as specified",
        ])
    elif assay_module == "SPLICE_MINIGENE":
        materials.extend([
            "Minigene vectors",
            "Transfection reagents",
        ])
    elif assay_module == "RNA_SEQ_REPORTER":
        materials.extend([
            "RNA extraction kit",
            "RNA-seq library preparation kit",
            "Illumina sequencing service",
            "Cell lines or biochemical reagents as specified",
        ])
    elif assay_module == "TRANSCRIPTIONAL_REPORTER":
        materials.extend([
            "Cell lines or biochemical reagents as specified",
        ])

    return materials


def _get_cro_notes(
    assay_module: AssayModuleId,
    catalog: AssayCatalog,
    variant_count: int
) -> str:
    """Generate CRO-specific notes and requirements."""
    module = catalog[assay_module]

    notes = [
        f"Estimated duration: {module.estimated_duration_weeks} weeks",
        f"Cost category: {module.cost_category}",
        f"Expected output: {', '.join(module.outputs)}",
    ]

    if variant_count > 20:
        notes.append("High-throughput screening recommended for large variant set")
    elif variant_count < 5:
        notes.append("Small variant set - consider pooling with other work packages")

    return ". ".join(notes)


def create_work_packages(
    variant_panel: VariantPanel,
    assignments: AssayAssignments,
    catalog: AssayCatalog
) -> WorkPackages:
    """Create work packages by grouping variants by gene × assay_module.

    Args:
        variant_panel: List of VariantPanelRow objects
        assignments: List of AssayAssignment objects
        catalog: Assay module catalog

    Returns:
        List of WorkPackage objects
    """
    LOGGER.info("Creating work packages from %d variants and %d assignments",
                len(variant_panel), len(assignments))

    # Create lookup dictionaries
    variant_lookup = {v.variant_id: v for v in variant_panel}
    assignment_lookup = {a.variant_id: a for a in assignments}

    # Group variants by gene × assay_module
    work_package_groups: Dict[tuple[str, AssayModuleId], List[str]] = {}

    for assignment in assignments:
        variant = variant_lookup.get(assignment.variant_id)
        if not variant:
            LOGGER.warning("Variant not found in panel: %s", assignment.variant_id)
            continue

        for assay_module in assignment.assay_modules:
            key = (variant.gene, assay_module)
            if key not in work_package_groups:
                work_package_groups[key] = []
            work_package_groups[key].append(assignment.variant_id)

    # Create work packages
    work_packages = []
    for (gene, assay_module), variant_ids in work_package_groups.items():
        variant_count = len(variant_ids)

        wp_id = _generate_work_package_id(gene, assay_module, variant_count)
        objective = _get_work_package_objective(gene, assay_module, catalog, variant_count)
        materials_provided = _get_materials_provided(gene, assay_module, catalog, variant_count)
        materials_needed = _get_materials_needed(assay_module, catalog)
        cro_notes = _get_cro_notes(assay_module, catalog, variant_count)

        work_package = WorkPackage(
            wp_id=wp_id,
            gene=gene,
            assay_module=assay_module,
            objective=objective,
            variant_ids=sorted(variant_ids),  # Sort for consistency
            materials_provided=materials_provided,
            materials_needed=materials_needed,
            cro_notes=cro_notes,
        )
        work_packages.append(work_package)

    # Sort work packages by ID for consistency
    work_packages.sort(key=lambda wp: wp.wp_id)

    LOGGER.info("Created %d work packages", len(work_packages))

    # Log work package statistics
    gene_counts = {}
    module_counts = {}
    for wp in work_packages:
        gene_counts[wp.gene] = gene_counts.get(wp.gene, 0) + 1
        module_counts[wp.assay_module] = module_counts.get(wp.assay_module, 0) + 1

    LOGGER.info("Work packages by gene: %s", gene_counts)
    LOGGER.info("Work packages by module: %s", module_counts)

    return work_packages


def save_work_packages(
    work_packages: WorkPackages,
    output_dir: Path = CRO_DATA_DIR
) -> Path:
    """Save WorkPackages to JSONL file."""
    output_dir.mkdir(parents=True, exist_ok=True)

    jsonl_path = output_dir / "work_packages.jsonl"
    with open(jsonl_path, "w") as f:
        for wp in work_packages:
            json.dump({
                "wp_id": wp.wp_id,
                "gene": wp.gene,
                "assay_module": wp.assay_module,
                "objective": wp.objective,
                "variant_ids": wp.variant_ids,
                "materials_provided": wp.materials_provided,
                "materials_needed": wp.materials_needed,
                "cro_notes": wp.cro_notes,
            }, f)
            f.write("\n")

    LOGGER.info("Saved %d work packages to %s", len(work_packages), jsonl_path)
    return jsonl_path


def load_work_packages() -> WorkPackages:
    """Load work packages from saved JSONL file."""
    jsonl_path = CRO_DATA_DIR / "work_packages.jsonl"
    if not jsonl_path.exists():
        raise FileNotFoundError(f"Work packages not found: {jsonl_path}. Run workpackages first.")

    work_packages = []
    with open(jsonl_path, "r") as f:
        for line in f:
            if line.strip():
                data = json.loads(line)
                wp = WorkPackage(
                    wp_id=data["wp_id"],
                    gene=data["gene"],
                    assay_module=data["assay_module"],
                    objective=data["objective"],
                    variant_ids=data["variant_ids"],
                    materials_provided=data["materials_provided"],
                    materials_needed=data["materials_needed"],
                    cro_notes=data["cro_notes"],
                )
                work_packages.append(wp)

    return work_packages


def main() -> None:
    """Run work package aggregation."""
    variant_panel = load_variant_panel()
    assignments = load_assay_assignments()
    catalog = load_assay_catalog()

    work_packages = create_work_packages(variant_panel, assignments, catalog)
    save_work_packages(work_packages)

    LOGGER.info("Successfully created %d work packages", len(work_packages))


if __name__ == "__main__":
    main()
