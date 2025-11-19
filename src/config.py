"""
Shared configuration for ABCA4 notebook pipeline.
Centralizes paths, constants, and validation logic.
"""

from pathlib import Path
import pandas as pd
import logging
from typing import Optional, Dict, Any
import os
import hashlib
import subprocess

# Initialize logging once when module is imported
LOGGER_NAME = "abca4"
_root_logger = logging.getLogger()
if not _root_logger.handlers:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        force=True
    )
logger = logging.getLogger(LOGGER_NAME)

# Core paths - all notebooks should import from here
def get_campaign_root() -> Path:
    """Get the campaign/repository root directory."""
    return Path(__file__).resolve().parents[1]


CAMPAIGN_ROOT = get_campaign_root()
DATA_RAW_DIR = CAMPAIGN_ROOT / "data_raw"
DATA_PROCESSED_DIR = CAMPAIGN_ROOT / "data_processed"

# Data directories
ANNOTATIONS_DIR = DATA_PROCESSED_DIR / "annotations"
FEATURES_DIR = DATA_PROCESSED_DIR / "features"
REPORTS_DIR = DATA_PROCESSED_DIR / "reports"

# Standard file names
ANNOTATED_VARIANTS_FILE = "abca4_vus_annotated.parquet"
SCORED_VARIANTS_FILE = "variants_scored.parquet"

# Demo mode settings
DEMO_MODE_ENABLED = False  # Set to True to enable synthetic data fallbacks


# LLM Configuration for Assay Draft Generation
GROQ_API_KEY = os.getenv("GROQ_API_KEY")
LLM_MODEL = os.getenv("LLM_MODEL", "llama-3.3-70b-versatile")
LLM_TEMP = float(os.getenv("LLM_TEMP", "0.2"))
LLM_MAX_TOKENS = int(os.getenv("LLM_MAX_TOKENS", "600"))
LLM_MAX_VARIANTS = int(os.getenv("LLM_MAX_VARIANTS", "12"))

# Provenance tracking
def get_git_sha() -> str:
    """Get current git SHA for provenance tracking."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            cwd=CAMPAIGN_ROOT
        )
        return result.stdout.strip() if result.returncode == 0 else "unknown"
    except Exception:
        return "unknown"

RUN_ID = os.getenv("RUN_ID", get_git_sha())

# Validate required LLM environment variables
def validate_llm_config():
    """Validate LLM configuration and fail fast if required vars missing."""
    if not GROQ_API_KEY:
        raise ValueError(
            "GROQ_API_KEY environment variable is required for assay draft generation.\n"
            "Please set GROQ_API_KEY in your environment and try again."
        )

    if not (0.1 <= LLM_TEMP <= 0.5):
        raise ValueError(f"LLM_TEMP must be between 0.1 and 0.5, got {LLM_TEMP}")

    if LLM_MAX_TOKENS > 1000:
        raise ValueError(f"LLM_MAX_TOKENS must be <= 1000 for cost control, got {LLM_MAX_TOKENS}")

    if LLM_MAX_VARIANTS > 20:
        raise ValueError(f"LLM_MAX_VARIANTS must be <= 20 for cost control, got {LLM_MAX_VARIANTS}")

    logger.info(f"LLM config validated: model={LLM_MODEL}, temp={LLM_TEMP}, max_tokens={LLM_MAX_TOKENS}, max_variants={LLM_MAX_VARIANTS}")

# Data contract for assay draft input
REQUIRED_VARIANT_COLUMNS = [
    "variant_id", "gene", "vep_consequence", "cluster_id",
    "impact_score", "gnomad_max_af"
]

ASSAY_DRAFTS_DIR = REPORTS_DIR / "assay_drafts"


def get_annotated_variants_path() -> Path:
    """Get path to annotated variants parquet file."""
    return ANNOTATIONS_DIR / ANNOTATED_VARIANTS_FILE


def get_scored_variants_path() -> Path:
    """Get path to scored variants parquet file."""
    return FEATURES_DIR / SCORED_VARIANTS_FILE


def validate_file_exists(path: Path, context: str = "") -> None:
    """Validate that a file exists, raise informative error if not."""
    if not path.exists():
        context_msg = f" ({context})" if context else ""
        raise FileNotFoundError(
            f"Required file missing: {path}{context_msg}\n"
            f"Please ensure upstream processing has completed."
        )


def validate_dataframe(df: pd.DataFrame, name: str, required_columns: Optional[list] = None) -> None:
    """Validate dataframe is not empty and has required columns."""
    if df.empty:
        raise ValueError(f"DataFrame '{name}' is empty")

    if required_columns:
        missing_cols = set(required_columns) - set(df.columns)
        if missing_cols:
            raise ValueError(f"DataFrame '{name}' missing required columns: {missing_cols}")


def load_parquet_safely(path: Path, name: str) -> pd.DataFrame:
    """Load parquet file with validation."""
    validate_file_exists(path, f"for {name}")
    df = pd.read_parquet(path)
    validate_dataframe(df, name)
    return df


# Conservation column mapping for consistency
CONSERVATION_COLUMN_MAP = {
    'phyloP100way': 'phylop_score',
    'phastCons100way': 'phastcons_score',
    'conservation_score': 'phylop_score'  # Use phyloP as primary conservation score
}


def standardize_conservation_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize conservation column names to phylop_score."""
    df = df.copy()
    for old_name, new_name in CONSERVATION_COLUMN_MAP.items():
        if old_name in df.columns and new_name not in df.columns:
            df = df.rename(columns={old_name: new_name})
    return df


def setup_logging():
    """Get a configured logger instance."""
    return logging.getLogger(__name__)


# ============================================================================
# Gene-Specific Configuration Management
# ============================================================================

def load_gene_config(gene_name: str) -> Dict[str, Any]:
    """
    Load gene configuration from YAML file.
    
    Args:
        gene_name: Gene symbol (e.g., "ABCA4", "CFTR")
        
    Returns:
        Dictionary containing all gene-specific settings
        
    Raises:
        FileNotFoundError: If config file for gene does not exist
        ValueError: If config file is invalid or missing required fields
    """
    config_dir = CAMPAIGN_ROOT / "config"
    config_path = config_dir / f"{gene_name.lower()}.yaml"
    
    if not config_path.exists():
        raise FileNotFoundError(
            f"Gene configuration not found: {config_path}\n"
            f"Available config directory: {config_dir}\n"
            f"Please create {config_path} with gene-specific settings."
        )
    
    try:
        import yaml
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    except Exception as e:
        raise ValueError(
            f"Failed to parse YAML config for {gene_name}: {e}\n"
            f"Config file: {config_path}"
        )
    
    if not isinstance(config, dict):
        raise ValueError(
            f"Invalid config format for {gene_name}: expected dict, got {type(config)}\n"
            f"Config file: {config_path}"
        )
    
    # Validate required top-level fields
    required_fields = ['gene_name', 'mode', 'features']
    missing_fields = [f for f in required_fields if f not in config]
    if missing_fields:
        raise ValueError(
            f"Gene config missing required fields: {missing_fields}\n"
            f"Config file: {config_path}"
        )
    
    logger.info(f"Loaded configuration for {config['gene_name']}")
    return config


def get_gene_config(gene_name: str, cache: Optional[Dict] = None) -> Dict[str, Any]:
    """
    Load or retrieve cached gene configuration.
    
    Args:
        gene_name: Gene symbol
        cache: Optional cache dict to store/retrieve configs
        
    Returns:
        Gene configuration dictionary
    """
    if cache is not None and gene_name in cache:
        return cache[gene_name]
    
    config = load_gene_config(gene_name)
    
    if cache is not None:
        cache[gene_name] = config
    
    return config


def validate_config_for_mode(config: Dict[str, Any]) -> bool:
    """
    Validate that config has required fields for its mode.
    
    Args:
        config: Gene configuration dictionary
        
    Returns:
        True if valid, False otherwise
    """
    mode = config.get('mode', 'clinvar')
    
    if mode == 'clinvar':
        required_clinvar = ['domains', 'scoring_weights', 'clustering', 'selection']
        missing = [f for f in required_clinvar if f not in config]
        if missing:
            logger.warning(
                f"ClinVar mode config missing fields: {missing}\n"
                f"Gene: {config.get('gene_name', 'unknown')}"
            )
            return False
    
    elif mode == 'mave':
        required_mave = ['scoring_weights', 'clustering', 'selection']
        missing = [f for f in required_mave if f not in config]
        if missing:
            logger.warning(
                f"MAVE mode config missing fields: {missing}\n"
                f"Gene: {config.get('gene_name', 'unknown')}"
            )
            return False
    
    return True
