"""Constants for CRO workflow slice."""

from pathlib import Path

# Paths
# Note: CAMPAIGN_ROOT is usually determined dynamically, but we can define relative paths here
# or expect CAMPAIGN_ROOT to be passed/imported. 
# For simplicity in this refactor, we'll keep the dynamic root resolution in the modules
# but define the relative paths here.

CATALOG_DIR_NAME = "catalog"
MECHANISM_RULES_FILENAME = "{gene}_mechanisms.yaml"
ASSAY_CATALOG_FILENAME = "assay_modules.yaml"

# Mechanism Defaults
DEFAULT_MECHANISM = "folding_stability"
DEFAULT_RATIONALE = "Conservative assignment"

# Data Directory Names
CRO_DATA_SUBDIR = "cro"
