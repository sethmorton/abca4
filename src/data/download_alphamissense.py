#!/usr/bin/env python3
"""
Download AlphaMissense scores for ABCA4.

Downloads AlphaMissense pathogenicity predictions.
"""

import subprocess
import sys
from pathlib import Path
from typing import Optional
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]

class AlphaMissenseDownloader:
    """Download AlphaMissense data."""

    # AlphaMissense lookup bucket (mirrors the upstream release)
    ALPHAMISSENSE_URL = "https://storage.googleapis.com/spliceai-lookup-reference-data/AlphaMissense_hg38.tsv.gz"
    COLUMN_NAMES = [
        'chrom', 'pos', 'ref', 'alt', 'genome', 'uniprot_id',
        'transcript_id', 'protein_variant', 'am_pathogenicity', 'am_class'
    ]
    ABCA4_UNIPROT = "P78363"
    ABCA4_TRANSCRIPTS = ("ENST00000370225",)

    def __init__(self, output_dir: Optional[Path] = None):
        default_dir = CAMPAIGN_ROOT / "data_raw" / "alphamissense"
        self.output_dir = output_dir or default_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def download_file(self, url: str, output_path: Path) -> bool:
        """Download a file using curl."""
        logger.info(f"Downloading {url} to {output_path}")

        try:
            cmd = [
                "curl", "-L", "-o", str(output_path),
                "--progress-bar", url
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"Download failed: {result.stderr}")
                return False

            logger.info(f"Downloaded {output_path.name} ({output_path.stat().st_size} bytes)")
            return True

        except Exception as e:
            logger.error(f"Download failed for {url}: {e}")
            return False

    def extract_abca4_data(self, tsv_path: Path, output_path: Path) -> bool:
        """Extract ABCA4-related data from AlphaMissense TSV."""
        logger.info("Extracting ABCA4 data from AlphaMissense...")

        try:
            import pandas as pd

            # Load AlphaMissense data (it's a large file, so we stream it)
            logger.info("Loading AlphaMissense TSV (this may take a while)...")

            # Read in chunks to handle large file
            chunk_size = 100000
            abca4_chunks = []

            for chunk in pd.read_csv(
                tsv_path,
                sep='\t',
                chunksize=chunk_size,
                comment='#',
                names=self.COLUMN_NAMES,
                header=None,
                dtype={
                    'chrom': str,
                    'pos': int,
                    'ref': str,
                    'alt': str,
                    'genome': str,
                    'uniprot_id': str,
                    'transcript_id': str,
                    'protein_variant': str,
                    'am_pathogenicity': float,
                    'am_class': str,
                }
            ):
                mask = chunk['uniprot_id'].astype(str).str.upper().eq(self.ABCA4_UNIPROT)
                transcript_mask = chunk['transcript_id'].astype(str).str.startswith(self.ABCA4_TRANSCRIPTS)
                mask = mask | transcript_mask

                abca4_chunk = chunk[mask]
                if not abca4_chunk.empty:
                    abca4_chunks.append(abca4_chunk)

            if not abca4_chunks:
                logger.warning("No ABCA4 data found in AlphaMissense")
                return False

            # Combine chunks
            abca4_df = pd.concat(abca4_chunks, ignore_index=True)
            logger.info(f"Found {len(abca4_df)} ABCA4 variants in AlphaMissense")

            # Save filtered data
            abca4_df.to_csv(output_path, sep='\t', index=False)
            logger.info(f"Saved ABCA4 AlphaMissense data to {output_path}")

            return True

        except Exception as e:
            logger.error(f"ABCA4 data extraction failed: {e}")
            return False

    def run(self) -> bool:
        """Run the complete AlphaMissense download and extraction process."""
        logger.info("Starting AlphaMissense download and extraction process...")

        # Download full AlphaMissense dataset
        alphamissense_path = self.output_dir / "AlphaMissense_hg38.tsv.gz"
        abca4_output_path = self.output_dir / "alphamissense_abca4.tsv"

        if not self.download_file(self.ALPHAMISSENSE_URL, alphamissense_path):
            logger.error("AlphaMissense download failed")
            return False

        # Extract ABCA4 data
        if not self.extract_abca4_data(alphamissense_path, abca4_output_path):
            logger.error("ABCA4 data extraction failed")
            return False

        logger.info("AlphaMissense download and extraction process completed successfully!")
        logger.info(f"AlphaMissense scores for ABCA4 saved to: {abca4_output_path}")
        return True

def main():
    """Main entry point."""
    downloader = AlphaMissenseDownloader()
    success = downloader.run()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
