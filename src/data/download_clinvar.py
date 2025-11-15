#!/usr/bin/env python3
"""
Download ClinVar data for ABCA4 campaign.

Downloads the latest GRCh38 ClinVar VCF and TSV files.
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

class ClinVarDownloader:
    """Download and verify ClinVar data files."""

    # Parameterized URLs - can be overridden via config
    CLINVAR_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar"

    def __init__(self, output_dir: Optional[Path] = None,
                 release_date: str = "20251109",
                 expected_vcf_size: Optional[int] = None,
                 expected_tsv_size: Optional[int] = None):
        default_dir = CAMPAIGN_ROOT / "data_raw" / "clinvar"
        self.output_dir = output_dir or default_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.release_date = release_date
        self.expected_vcf_size = expected_vcf_size or 181299219  # Default for 20251109
        self.expected_tsv_size = expected_tsv_size

        # Build URLs from parameters
        self.CLINVAR_VCF_URL = f"{self.CLINVAR_BASE_URL}/vcf_GRCh38/clinvar_{self.release_date}.vcf.gz"
        self.CLINVAR_VCF_TBI_URL = f"{self.CLINVAR_BASE_URL}/vcf_GRCh38/clinvar_{self.release_date}.vcf.gz.tbi"
        self.CLINVAR_TSV_URL = f"{self.CLINVAR_BASE_URL}/tab_delimited/variant_summary.txt.gz"

    def download_file(self, url: str, output_path: Path,
                     expected_size: Optional[int] = None,
                     checksum: Optional[str] = None,
                     resume: bool = True) -> bool:
        """Download a file with progress tracking, size verification, checksum validation, and resume support."""
        logger.info(f"Downloading {url} to {output_path}")

        try:
            # Use curl with resume support if file exists
            cmd = ["curl", "-L"]

            if resume and output_path.exists():
                cmd.extend(["-C", "-"])  # Resume from where it left off
                logger.info(f"Resuming download from {output_path.stat().st_size} bytes")

            cmd.extend([
                "-o", str(output_path),
                "--progress-bar",
                "--retry", "3",
                "--retry-delay", "5",
                url
            ])

            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"Download failed: {result.stderr}")
                return False

            # Verify file size if expected
            actual_size = output_path.stat().st_size
            if expected_size and actual_size != expected_size:
                logger.warning(f"File size mismatch for {output_path.name}: "
                              f"expected {expected_size}, got {actual_size}")
                # Don't fail on size mismatch - some files may have updated sizes

            # Verify checksum if provided
            if checksum:
                import hashlib
                sha256 = hashlib.sha256()
                with open(output_path, 'rb') as f:
                    for chunk in iter(lambda: f.read(4096), b""):
                        sha256.update(chunk)
                actual_checksum = sha256.hexdigest()

                if actual_checksum != checksum:
                    logger.error(f"Checksum mismatch for {output_path.name}: "
                                f"expected {checksum}, got {actual_checksum}")
                    return False
                else:
                    logger.info(f"Checksum verified for {output_path.name}")

            logger.info(f"Successfully downloaded {output_path.name} "
                       f"({actual_size} bytes)")
            return True

        except Exception as e:
            logger.error(f"Download failed for {url}: {e}")
            return False

    def download_vcf(self, vcf_checksum: Optional[str] = None,
                    tbi_checksum: Optional[str] = None) -> bool:
        """Download ClinVar VCF and index."""
        vcf_filename = f"clinvar_{self.release_date}.vcf.gz"
        tbi_filename = f"clinvar_{self.release_date}.vcf.gz.tbi"

        vcf_path = self.output_dir / vcf_filename
        tbi_path = self.output_dir / tbi_filename

        # Download VCF
        if not self.download_file(self.CLINVAR_VCF_URL, vcf_path,
                                self.expected_vcf_size, vcf_checksum):
            return False

        # Download TBI index
        if not self.download_file(self.CLINVAR_VCF_TBI_URL, tbi_path,
                                checksum=tbi_checksum):
            logger.warning("TBI index download failed, but VCF was successful")

        return True

    def download_tsv(self, tsv_checksum: Optional[str] = None) -> bool:
        """Download ClinVar TSV summary."""
        tsv_path = self.output_dir / "variant_summary.txt.gz"

        return self.download_file(self.CLINVAR_TSV_URL, tsv_path,
                                self.expected_tsv_size, tsv_checksum)

    def verify_downloads(self) -> bool:
        """Verify all downloads completed successfully."""
        required_files = [
            f"clinvar_{self.release_date}.vcf.gz",
            "variant_summary.txt.gz"
        ]

        missing_files = []
        for filename in required_files:
            filepath = self.output_dir / filename
            if not filepath.exists():
                missing_files.append(filename)

        if missing_files:
            logger.error(f"Missing required files: {missing_files}")
            return False

        # Check file sizes
        vcf_path = self.output_dir / f"clinvar_{self.release_date}.vcf.gz"
        actual_vcf_size = vcf_path.stat().st_size
        if self.expected_vcf_size and actual_vcf_size != self.expected_vcf_size:
            logger.warning(f"VCF size mismatch: {actual_vcf_size} vs {self.expected_vcf_size}")

        # Quick validation of VCF format
        try:
            import gzip
            with gzip.open(vcf_path, 'rt') as f:
                header_lines = []
                for line in f:
                    if line.startswith('#'):
                        header_lines.append(line.strip())
                    else:
                        # Found first data line, check basic format
                        fields = line.strip().split('\t')
                        if len(fields) >= 8:
                            logger.info(f"VCF validation passed - {len(header_lines)} header lines")
                            break
                        else:
                            logger.error("VCF format validation failed")
                            return False
        except Exception as e:
            logger.error(f"VCF validation failed: {e}")
            return False

        logger.info("All ClinVar downloads verified!")
        return True

    def run(self, vcf_checksum: Optional[str] = None,
            tbi_checksum: Optional[str] = None,
            tsv_checksum: Optional[str] = None) -> bool:
        """Run the complete download process."""
        logger.info(f"Starting ClinVar download process for release {self.release_date}...")

        success = True

        # Download VCF files
        if not self.download_vcf(vcf_checksum, tbi_checksum):
            logger.error("VCF download failed")
            success = False

        # Download TSV summary
        if not self.download_tsv(tsv_checksum):
            logger.error("TSV download failed")
            success = False

        # Verify downloads
        if not self.verify_downloads():
            logger.error("Download verification failed")
            success = False

        if success:
            logger.info("ClinVar download process completed successfully!")
            logger.info(f"Files saved to: {self.output_dir}")
            logger.info(f"Release date: {self.release_date}")
        else:
            logger.error("ClinVar download process failed!")

        return success

def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Download ClinVar data")
    parser.add_argument("--release-date", default="20251109",
                       help="ClinVar release date (default: 20251109)")
    parser.add_argument("--output-dir", type=Path,
                       help="Output directory (default: data_raw/clinvar)")
    parser.add_argument("--vcf-checksum", help="Expected SHA256 checksum for VCF file")
    parser.add_argument("--tbi-checksum", help="Expected SHA256 checksum for TBI index")
    parser.add_argument("--tsv-checksum", help="Expected SHA256 checksum for TSV file")
    parser.add_argument("--expected-vcf-size", type=int,
                       help="Expected VCF file size in bytes")
    parser.add_argument("--expected-tsv-size", type=int,
                       help="Expected TSV file size in bytes")

    args = parser.parse_args()

    downloader = ClinVarDownloader(
        output_dir=args.output_dir,
        release_date=args.release_date,
        expected_vcf_size=args.expected_vcf_size,
        expected_tsv_size=args.expected_tsv_size
    )

    success = downloader.run(
        vcf_checksum=args.vcf_checksum,
        tbi_checksum=args.tbi_checksum,
        tsv_checksum=args.tsv_checksum
    )
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
