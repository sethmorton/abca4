"""Type definitions for CRO study plan generation.

All data shapes are gene-agnostic and use strict typing without Any types.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, List, Literal, Set, TypedDict


# Core controlled vocabularies
ConsequenceType = Literal[
    "missense", "nonsense", "synonymous", "splice_region", "intronic_splice",
    "frameshift", "inframe_indel", "utr", "regulatory", "other"
]

DomainType = Literal["NBD1", "NBD2", "TMD", "CTD", "other"]

AssayModuleId = Literal[
    "DSF_SEC", "FUNCTIONAL", "TRAFFICKING", "SPLICE_MINIGENE",
    "RNA_SEQ_REPORTER", "TRANSCRIPTIONAL_REPORTER"
]

MechanismTag = Literal[
    "folding_stability", "trafficking_localization", "transport_activity",
    "splicing_defect", "regulatory_transcriptional"
]

DesignType = Literal["full_factorial", "one_factor", "fractional"]


@dataclass(frozen=True)
class VariantPanelRow:
    """Row in the variant panel from Stage 1 parsing."""
    variant_id: str
    gene: str
    rank: int
    chrom: str
    pos: int
    ref: str
    alt: str
    consequence: ConsequenceType
    domain: DomainType
    model_score: float
    impact_score: float
    cluster_id: str
    coverage: float
    assay_hint: Literal["DSF_SEC", "FUNCTIONAL", "SPLICE_MINIGENE", "RNA_SEQ_REPORTER", "TRAFFICKING", "TRANSCRIPTIONAL_REPORTER", "UNKNOWN"]
    selection_strategy: str
    panel_size: int
    panel_id: str
    source_md_path: Path


@dataclass(frozen=True)
class MechanismAnnotation:
    """Mechanism annotation from Stage 2."""
    variant_id: str
    mechanism_tags: Set[MechanismTag]
    rationale: str


@dataclass(frozen=True)
class AssayAssignment:
    """Assay assignment from Stage 3."""
    variant_id: str
    assay_modules: Set[AssayModuleId]
    rationale: str


@dataclass(frozen=True)
class WorkPackage:
    """Work package from Stage 4."""
    wp_id: str
    gene: str
    assay_module: AssayModuleId
    objective: str
    variant_ids: list[str]
    materials_provided: list[str]
    materials_needed: list[str]
    cro_notes: str


@dataclass(frozen=True)
class Factor:
    """Experimental factor definition."""
    name: str
    levels: list[str | float]


@dataclass(frozen=True)
class DesignSummary:
    """Experimental design from Stage 5."""
    wp_id: str
    factors: list[Factor]
    design_type: DesignType
    tech_reps: int
    bio_reps: int
    controls: list[str]


@dataclass(frozen=True)
class DeliverableSpec:
    """Deliverable specification from Stage 6."""
    wp_id: str
    primary_metrics: list[str]
    raw_returns: list[str]
    summary_columns: list[str]
    qc_expectations: list[str]


@dataclass(frozen=True)
class AssayModule:
    """Assay module definition from catalog."""
    id: AssayModuleId
    name: str
    description: str
    inputs: list[str]
    outputs: list[str]
    estimated_duration_weeks: int
    cost_category: Literal["low", "medium", "high"]
    mechanism_mapping: Set[MechanismTag]


# TypedDict definitions for parsed/config data structures
class MechanismRuleCondition(TypedDict, total=False):
    """Condition for mechanism tagging rules."""
    consequence: ConsequenceType | List[ConsequenceType]
    domain: DomainType | List[DomainType]
    impact_score: List[float]  # [min, max] range

class MechanismRule(TypedDict):
    """Mechanism tagging rule configuration."""
    condition: MechanismRuleCondition
    mechanism: MechanismTag | List[MechanismTag]
    rationale: str

class MechanismRulesConfig(TypedDict):
    """Complete mechanism rules configuration."""
    rules: List[MechanismRule]
    default_mechanism: MechanismTag
    default_rationale: str

class AssayModuleConfig(TypedDict):
    """Assay module configuration from YAML."""
    name: str
    description: str
    inputs: List[str]
    outputs: List[str]
    estimated_duration_weeks: int
    cost_category: Literal["low", "medium", "high"]
    mechanism_mapping: List[MechanismTag]

class AssayCatalogConfig(TypedDict):
    """Complete assay catalog configuration."""
    DSF_SEC: AssayModuleConfig
    FUNCTIONAL: AssayModuleConfig
    TRAFFICKING: AssayModuleConfig
    SPLICE_MINIGENE: AssayModuleConfig
    RNA_SEQ_REPORTER: AssayModuleConfig
    TRANSCRIPTIONAL_REPORTER: AssayModuleConfig

class WorkPackageData(TypedDict):
    """Work package data from JSONL."""
    wp_id: str
    gene: str
    assay_module: AssayModuleId
    objective: str
    variant_ids: List[str]
    materials_provided: List[str]
    materials_needed: List[str]
    cro_notes: str

class DesignSummaryData(TypedDict):
    """Design summary data from JSON."""
    wp_id: str
    factors: List[FactorData]  # Factor definitions with name/levels
    design_type: DesignType
    tech_reps: int
    bio_reps: int
    controls: List[str]

class DeliverableSpecData(TypedDict):
    """Deliverable specification data from JSON."""
    wp_id: str
    primary_metrics: List[str]
    raw_returns: List[str]
    summary_columns: List[str]
    qc_expectations: List[str]

class FactorData(TypedDict):
    """Factor definition with name and levels."""
    name: str
    levels: List[str | float | int]

class AssayAssignmentData(TypedDict):
    """Assay assignment data from JSON."""
    variant_id: str
    assay_modules: List[AssayModuleId]
    rationale: str

class MechanismAnnotationData(TypedDict):
    """Mechanism annotation data from JSON."""
    variant_id: str
    mechanism_tags: List[MechanismTag]
    rationale: str

class VariantPanelRowData(TypedDict):
    """Variant panel row data from CSV/parquet."""
    variant_id: str
    gene: str
    rank: int
    chrom: str
    pos: int
    ref: str
    alt: str
    consequence: ConsequenceType
    domain: DomainType
    model_score: float
    impact_score: float
    cluster_id: str
    coverage: float
    assay_hint: AssayModuleId | Literal["UNKNOWN"]
    selection_strategy: str
    panel_size: int
    panel_id: str
    source_md_path: str

class VariantReportMetadata(TypedDict):
    """Metadata extracted from variant selection report."""
    gene: str
    selection_strategy: str
    panel_size: int
    coverage_target: float

# Type aliases for collections
VariantPanel = list[VariantPanelRow]
MechanismPanel = list[MechanismAnnotation]
AssayAssignments = list[AssayAssignment]
WorkPackages = list[WorkPackage]
DesignSummaries = list[DesignSummary]
DeliverableSpecs = list[DeliverableSpec]
AssayCatalog = dict[AssayModuleId, AssayModule]

# Type aliases for parsed data
WorkPackagesData = List[WorkPackageData]
DesignSummariesData = List[DesignSummaryData]
DeliverableSpecsData = List[DeliverableSpecData]
AssayAssignmentsData = List[AssayAssignmentData]
MechanismPanelData = List[MechanismAnnotationData]
