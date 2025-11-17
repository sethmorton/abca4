"""
Shared configuration for ABCA4 notebook pipeline.
Centralizes paths, constants, and validation logic.
"""

from pathlib import Path
import pandas as pd
import logging
from typing import Optional, Dict, Any

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
