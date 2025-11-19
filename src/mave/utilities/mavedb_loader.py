"""
Load real MaveDB CSV files from the Zenodo dump and convert to our canonical schema.
Handles MaveDB's CSV format directly without needing individual downloads.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import re

logger = logging.getLogger(__name__)


def load_mavedb_csv(csv_path: Path) -> pd.DataFrame:
    """
    Load a single MaveDB CSV file and parse it.
    
    Args:
        csv_path: Path to a .scores.csv file from MaveDB dump
    
    Returns:
        DataFrame with columns: hgvs_nt, hgvs_pro, score, and metadata
    """
    try:
        df = pd.read_csv(csv_path, dtype={'hgvs_pro': str, 'hgvs_nt': str})
        
        # Filter to rows with valid HGVS protein notation
        if 'hgvs_pro' in df.columns:
            df = df[df['hgvs_pro'].notna()].copy()
        
        return df
    except Exception as e:
        logger.error(f"Failed to load {csv_path}: {e}")
        return None


def parse_hgvs_protein(hgvs_str: str) -> dict:
    """
    Parse HGVS protein notation like p.Met1Val or M1V.
    
    Returns: {pos, wt_aa, mut_aa} or None if unparseable
    """
    if not hgvs_str or not isinstance(hgvs_str, str):
        return None
    
    # Remove 'p.' prefix if present
    hgvs_str = hgvs_str.replace('p.', '').strip()
    
    three_letter_to_one = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'STOP': '*'
    }
    
    # Try three-letter format: MetXXXVal -> M XXX V
    for wt_3, wt_1 in three_letter_to_one.items():
        for mut_3, mut_1 in three_letter_to_one.items():
            pattern = f"^{wt_3}(\\d+){mut_3}$"
            match = re.match(pattern, hgvs_str, re.IGNORECASE)
            if match:
                pos_str = match.group(1)
                try:
                    pos = int(pos_str)
                    return {'pos': pos, 'wt_aa': wt_1, 'mut_aa': mut_1}
                except:
                    pass
    
    # Try one-letter format: M1V
    match = re.match(r"^([A-Z\*])(\d+)([A-Z\*])$", hgvs_str)
    if match:
        wt_aa, pos_str, mut_aa = match.groups()
        try:
            pos = int(pos_str)
            # Normalize
            wt_aa = wt_aa.replace('*', '*').upper()
            mut_aa = mut_aa.replace('*', '*').upper()
            return {'pos': pos, 'wt_aa': wt_aa, 'mut_aa': mut_aa}
        except:
            pass
    
    return None


def ingest_mavedb_csv(
    csv_path: Path,
    dataset_id: str,
    gene: str,
    mavedb_accession: str,
    assay_type: str = "unknown",
) -> pd.DataFrame:
    """
    Convert a MaveDB CSV to canonical schema.
    
    Args:
        csv_path: Path to .scores.csv file
        dataset_id: Name for this dataset (e.g., "TP53_Saturation2016")
        gene: Gene symbol (e.g., "TP53")
        mavedb_accession: URN or accession ID
        assay_type: Type of assay
    
    Returns:
        DataFrame in canonical schema
    """
    logger.info(f"Loading {csv_path.name}...")
    df = load_mavedb_csv(csv_path)
    
    if df is None or len(df) == 0:
        logger.error(f"No data loaded from {csv_path}")
        return None
    
    logger.info(f"Loaded {len(df)} rows")
    
    # Parse HGVS
    parsed = df['hgvs_pro'].apply(parse_hgvs_protein)
    valid_mask = parsed.notna()
    
    logger.info(f"Parsed {valid_mask.sum()} / {len(df)} variants successfully")
    
    # Extract parsed data
    rows = []
    for idx, (row_idx, row) in enumerate(df[valid_mask].iterrows()):
        p = parsed[row_idx]
        
        if pd.isna(row['score']) or row['score'] == '':
            continue
        
        try:
            score = float(row['score'])
        except:
            continue
        
        rows.append({
            'dataset_id': dataset_id,
            'gene': gene,
            'wt_aa': p['wt_aa'],
            'pos': p['pos'],
            'mut_aa': p['mut_aa'],
            'variant_id': f"{gene}:{p['pos']}:{p['wt_aa']}>{p['mut_aa']}",
            'functional_score_raw': score,
            'functional_score': score,  # Will be normalized later
            'assay_type': assay_type,
            'mavedb_accession': mavedb_accession,
            'metadata_json': None,
        })
    
    if not rows:
        logger.warning(f"No valid variants extracted from {csv_path}")
        return None
    
    result_df = pd.DataFrame(rows)
    logger.info(f"Extracted {len(result_df)} canonical variants")
    
    return result_df


def find_mavedb_files(
    dump_dir: Path = Path("data_raw/mave/mavedb-dump.20250612164404"),
    limit: int = None,
) -> list:
    """
    Find all score CSV files in MaveDB dump directory.
    
    Args:
        dump_dir: Path to mavedb-dump directory
        limit: Max number of files to return (for testing)
    
    Returns:
        List of (csv_path, dataset_id) tuples
    """
    csv_dir = dump_dir / "csv"
    
    if not csv_dir.exists():
        logger.error(f"MaveDB dump not found at {csv_dir}")
        return []
    
    # Find all score files
    score_files = sorted(list(csv_dir.glob("*.scores.csv")))
    
    logger.info(f"Found {len(score_files)} score files in {csv_dir}")
    
    if limit:
        score_files = score_files[:limit]
        logger.info(f"Limited to {limit} files for testing")
    
    # Return as dataset descriptors
    datasets = []
    for csv_file in score_files:
        # Extract URN from filename
        urn_match = re.search(r'urn-mavedb-(\d+)-[a-z0-9]+-\d+', csv_file.name)
        if urn_match:
            dataset_id = f"MAVEDB_{urn_match.group(1)}"
            datasets.append({
                'csv_path': csv_file,
                'dataset_id': dataset_id,
                'accession': f"urn:mavedb:{urn_match.group(1)}",
            })
    
    return datasets


def auto_generate_config(
    dump_dir: Path = Path("data_raw/mave/mavedb-dump.20250612164404"),
    limit: int = 10,
) -> str:
    """
    Scan MaveDB dump and generate config YAML for top datasets by size.
    """
    datasets = find_mavedb_files(dump_dir, limit=None)
    
    # Score each file by size
    ranked = []
    for ds in datasets:
        try:
            df = pd.read_csv(ds['csv_path'], nrows=1)
            size = len(pd.read_csv(ds['csv_path']))
            ranked.append((ds, size))
        except:
            continue
    
    # Sort by size
    ranked.sort(key=lambda x: x[1], reverse=True)
    
    # Take top N
    top_datasets = [ds for ds, size in ranked[:limit]]
    
    # Generate YAML
    yaml_lines = ["datasets:"]
    
    for ds in top_datasets:
        yaml_lines.append(f"\n  - dataset_id: {ds['dataset_id']}")
        yaml_lines.append(f"    gene: UNKNOWN  # TODO: detect from CSV")
        yaml_lines.append(f"    mavedb_accession: {ds['accession']}")
        # Use absolute path or relative from project root
        rel_path = ds['csv_path']
        try:
            rel_path = ds['csv_path'].relative_to(Path.cwd())
        except ValueError:
            rel_path = ds['csv_path']
        yaml_lines.append(f"    raw_data_path: {rel_path}")
        yaml_lines.append(f"    assay_type: unknown")
        yaml_lines.append(f"    raw_score_higher_is_better: false  # TODO: verify")
        yaml_lines.append(f"    normalization:")
        yaml_lines.append(f"      method: zscore")
        yaml_lines.append(f"      flip_sign: false")
        yaml_lines.append(f"      clip_low: null")
        yaml_lines.append(f"      clip_high: null")
        yaml_lines.append(f"    hit_definition:")
        yaml_lines.append(f"      direction: low_is_loss_of_function")
        yaml_lines.append(f"      percentile_threshold: 0.2")
        yaml_lines.append(f"    notes: \"MaveDB accession {ds['accession']}\"")
    
    return "\n".join(yaml_lines)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    # Test loading a few datasets
    datasets = find_mavedb_files(limit=3)
    print(f"Found {len(datasets)} datasets\n")
    
    for ds in datasets[:2]:
        print(f"\n{ds['dataset_id']}:")
        df = ingest_mavedb_csv(
            ds['csv_path'],
            ds['dataset_id'],
            "UNKNOWN",
            ds['accession'],
        )
        if df is not None:
            print(f"  Shape: {df.shape}")
            print(f"  Score range: [{df['functional_score_raw'].min():.2f}, {df['functional_score_raw'].max():.2f}]")

