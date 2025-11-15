#!/usr/bin/env python3
"""
Download and extract gnomAD data for ABCA4 region.

Downloads gnomAD v4.1.0 data and extracts variants in the ABCA4 region (chr1:94.4M-95.2M).
"""

import logging
import subprocess
import sys
from pathlib import Path
from typing import Optional
from urllib.parse import quote

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

CAMPAIGN_ROOT = Path(__file__).resolve().parents[2]

class GnomADDownloader:
    """Download and extract gnomAD data for specified region."""

    # Default ABCA4 region: chr1:93,500,000-95,000,000 (±500kb from gene)
    # ABCA4 gene coordinates (GRCh38): chr1:94,040,717-94,406,815
    DEFAULT_CHROM = "1"
    DEFAULT_START = 93500000  # 93,500,000 (500kb upstream)
    DEFAULT_END = 95000000    # 95,000,000 (500kb downstream)

    # gnomAD base URLs
    GNOMAD_BASE_URL = "https://storage.googleapis.com/download/storage/v1/b/gcp-public-data--gnomad/o"

    def __init__(self, output_dir: Optional[Path] = None,
                 version: str = "4.1",
                 chrom: str = DEFAULT_CHROM,
                 start: Optional[int] = None,
                 end: Optional[int] = None,
                 stream_only: bool = True,
                 chr_prefix: bool = True):
        default_dir = CAMPAIGN_ROOT / "data_raw" / "gnomad"
        self.output_dir = output_dir or default_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.version = version
        self.chr_prefix = chr_prefix
        self.chrom = chrom
        self.chrom_label = chrom if chrom.startswith("chr") else (
            f"chr{chrom}" if chr_prefix else chrom
        )
        self.start = start or self.DEFAULT_START
        self.end = end or self.DEFAULT_END
        self.region = f"{self.chrom_label}:{self.start}-{self.end}"
        self.stream_only = stream_only

        # Build URLs from parameters
        genome_object = (
            f"release/{self.version}/vcf/genomes/"
            f"gnomad.genomes.v{self.version}.sites.{self.chrom_label}.vcf.bgz"
        )
        exome_object = (
            f"release/{self.version}/vcf/exomes/"
            f"gnomad.exomes.v{self.version}.sites.{self.chrom_label}.vcf.bgz"
        )

        self.GNOMAD_GENOME_URL = self._build_gcs_url(genome_object)
        self.GNOMAD_GENOME_TBI_URL = self._build_gcs_url(f"{genome_object}.tbi")

        self.GNOMAD_EXOME_URL = self._build_gcs_url(exome_object)
        self.GNOMAD_EXOME_TBI_URL = self._build_gcs_url(f"{exome_object}.tbi")

    def _build_gcs_url(self, object_path: str) -> str:
        """Return a media download URL for the requested GCS object."""
        encoded = quote(object_path, safe="")
        return f"{self.GNOMAD_BASE_URL}/{encoded}?alt=media"

    def download_file(self, url: str, output_path: Path,
                     expected_size: Optional[int] = None,
                     checksum: Optional[str] = None,
                     resume: bool = True) -> bool:
        """Download a file using curl with resume and verification support."""
        logger.info(f"Downloading {url} to {output_path}")

        try:
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

            logger.info(f"Downloaded {output_path.name} ({actual_size} bytes)")
            return True

        except Exception as e:
            logger.error(f"Download failed for {url}: {e}")
            return False

    def stream_region(self, source_url: str, output_path: Path,
                      dataset_type: str) -> bool:
        """Stream a genomic window directly from a remote bgzipped VCF via tabix."""
        logger.info(f"Streaming {self.region} from {source_url} ({dataset_type})")

        cmd = [
            "tabix", "-h", source_url, self.region
        ]

        try:
            with open(output_path, 'w') as handle:
                result = subprocess.run(cmd, stdout=handle, stderr=subprocess.PIPE, text=True)

            if result.returncode != 0:
                logger.error(f"tabix streaming failed: {result.stderr}")
                output_path.unlink(missing_ok=True)
                return False

            data_lines = 0
            with open(output_path, 'r') as handle:
                for line in handle:
                    if not line.startswith('#'):
                        data_lines += 1

            logger.info(f"Streamed {data_lines} variants for {dataset_type}")
            return True
        except FileNotFoundError:
            logger.error("tabix is required but not found on PATH")
            return False
        except Exception as exc:
            logger.error(f"Streaming failed: {exc}")
            return False

    def extract_region(self, vcf_path: Path, output_path: Path, dataset_type: str) -> bool:
        """Extract specified region from gnomAD VCF using tabix/bcftools."""
        logger.info(f"Extracting {self.region} from {vcf_path.name}")

        try:
            # Use tabix to extract region (assuming TBI index exists)
            cmd = [
                "tabix", "-h", str(vcf_path), self.region
            ]

            with open(output_path, 'w') as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)

            if result.returncode != 0:
                logger.error(f"Region extraction failed: {result.stderr}")
                # Try alternative: bcftools view if tabix fails
                logger.info("Trying bcftools view as fallback...")
                cmd = [
                    "bcftools", "view", "-r", self.region,
                    "-o", str(output_path), str(vcf_path)
                ]
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    logger.error(f"bcftools extraction also failed: {result.stderr}")
                    return False

            # Check if we got any variants
            with open(output_path, 'r') as f:
                lines = f.readlines()
                data_lines = [line for line in lines if not line.startswith('#')]
                logger.info(f"Extracted {len(data_lines)} variants from {dataset_type}")

                if len(data_lines) == 0:
                    logger.warning(f"No variants found in region {self.region} for {dataset_type}")
                    return True  # Not an error, just no variants in region

            logger.info(f"Successfully extracted region to {output_path.name}")
            return True

        except Exception as e:
            logger.error(f"Region extraction failed: {e}")
            return False

    def download_and_extract_genome(self, vcf_checksum: Optional[str] = None,
                                   tbi_checksum: Optional[str] = None) -> bool:
        """Fetch genome data either by streaming or downloading the full chromosome."""
        logger.info("Processing gnomAD genome data...")

        if self.stream_only:
            output_path = self.output_dir / f"gnomad_v{self.version}_abca4_genome.vcf"
            return self.stream_region(self.GNOMAD_GENOME_URL, output_path, "genome")

        genome_vcf_filename = (
            f"gnomad.genomes.v{self.version}.sites.{self.chrom_label}.vcf.bgz"
        )
        genome_tbi_filename = f"{genome_vcf_filename}.tbi"
        genome_vcf = self.output_dir / genome_vcf_filename
        genome_tbi = self.output_dir / genome_tbi_filename

        if not self.download_file(self.GNOMAD_GENOME_URL, genome_vcf, checksum=vcf_checksum):
            return False

        if not self.download_file(self.GNOMAD_GENOME_TBI_URL, genome_tbi, checksum=tbi_checksum):
            logger.warning("TBI index download failed, attempting extraction anyway")

        output_path = self.output_dir / f"gnomad_v{self.version}_abca4_genome.vcf"
        return self.extract_region(genome_vcf, output_path, "genome")

    def download_and_extract_exome(self, vcf_checksum: Optional[str] = None,
                                  tbi_checksum: Optional[str] = None) -> bool:
        """Fetch exome data either by streaming or downloading the full chromosome."""
        logger.info("Processing gnomAD exome data...")

        if self.stream_only:
            output_path = self.output_dir / f"gnomad_v{self.version}_abca4_exome.vcf"
            return self.stream_region(self.GNOMAD_EXOME_URL, output_path, "exome")

        exome_vcf_filename = (
            f"gnomad.exomes.v{self.version}.sites.{self.chrom_label}.vcf.bgz"
        )
        exome_tbi_filename = f"{exome_vcf_filename}.tbi"
        exome_vcf = self.output_dir / exome_vcf_filename
        exome_tbi = self.output_dir / exome_tbi_filename

        if not self.download_file(self.GNOMAD_EXOME_URL, exome_vcf, checksum=vcf_checksum):
            return False

        if not self.download_file(self.GNOMAD_EXOME_TBI_URL, exome_tbi, checksum=tbi_checksum):
            logger.warning("TBI index download failed, attempting extraction anyway")

        output_path = self.output_dir / f"gnomad_v{self.version}_abca4_exome.vcf"
        return self.extract_region(exome_vcf, output_path, "exome")

    def verify_downloads(self) -> bool:
        """Verify all extractions completed successfully."""
        required_files = []

        if self.stream_only:
            required_files = [
                f"gnomad_v{self.version}_abca4_genome.vcf",
                f"gnomad_v{self.version}_abca4_exome.vcf"
            ]
        else:
            required_files = [
                f"gnomad.v{self.version}.sites.chr{self.chrom}.vcf.bgz",
                f"gnomad_exome.v{self.version}.sites.chr{self.chrom}.vcf.bgz"
            ]

        missing_files = []
        for filename in required_files:
            filepath = self.output_dir / filename
            if not filepath.exists():
                missing_files.append(filename)

        if missing_files:
            logger.error(f"Missing required files: {missing_files}")
            return False

        # Check file sizes and variant counts
        for filename in required_files:
            filepath = self.output_dir / filename
            try:
                if filename.endswith('.vcf'):
                    data_lines = 0
                    with open(filepath, 'r') as handle:
                        for line in handle:
                            if not line.startswith('#'):
                                data_lines += 1
                    logger.info(f"{filename}: {data_lines} variants")
                else:
                    logger.info(f"{filename}: {filepath.stat().st_size} bytes")
            except Exception as e:
                logger.warning(f"Could not read {filename}: {e}")

        logger.info("All gnomAD downloads verified!")
        return True

    def run(self, genome_vcf_checksum: Optional[str] = None,
            genome_tbi_checksum: Optional[str] = None,
            exome_vcf_checksum: Optional[str] = None,
            exome_tbi_checksum: Optional[str] = None) -> bool:
        """Run the complete gnomAD download and extraction process."""
        logger.info("Starting gnomAD download and extraction process...")
        logger.info(f"Version: {self.version}, Chromosome: {self.chrom}, Region: {self.region}")

        success = True

        # Process genome data
        if not self.download_and_extract_genome(genome_vcf_checksum, genome_tbi_checksum):
            logger.error("Genome data processing failed")
            success = False

        # Process exome data
        if not self.download_and_extract_exome(exome_vcf_checksum, exome_tbi_checksum):
            logger.error("Exome data processing failed")
            success = False

        # Verify results
        if not self.verify_downloads():
            logger.error("Download verification failed")
            success = False

        if success:
            logger.info("gnomAD download and extraction process completed successfully!")
            if self.stream_only:
                logger.info(f"Region-level VCFs stored under: {self.output_dir}")
            else:
                logger.info(f"Full chromosome archives stored under: {self.output_dir}")
        else:
            logger.error("gnomAD download and extraction process failed!")

        return success

def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Download and extract gnomAD data")
    parser.add_argument("--version", default="4.1",
                       help="gnomAD version (default: 4.1)")
    parser.add_argument("--chrom", default="1",
                       help="Chromosome number (default: 1)")
    parser.add_argument("--start", type=int,
                       help="Region start position (default: 94400000 for ABCA4)")
    parser.add_argument("--end", type=int,
                       help="Region end position (default: 95200000 for ABCA4)")
    parser.add_argument("--output-dir", type=Path,
                       help="Output directory (default: campaigns/abca4/data_raw/gnomad)")
    parser.add_argument("--download-full", action="store_true",
                       help="Download entire chromosomes instead of streaming the region")
    parser.add_argument("--no-chr-prefix", action="store_true",
                       help="Use contigs without implicit 'chr' prefix")
    parser.add_argument("--genome-vcf-checksum", help="Expected SHA256 checksum for genome VCF")
    parser.add_argument("--genome-tbi-checksum", help="Expected SHA256 checksum for genome TBI")
    parser.add_argument("--exome-vcf-checksum", help="Expected SHA256 checksum for exome VCF")
    parser.add_argument("--exome-tbi-checksum", help="Expected SHA256 checksum for exome TBI")

    args = parser.parse_args()

    downloader = GnomADDownloader(
        output_dir=args.output_dir,
        version=args.version,
        chrom=args.chrom,
        start=args.start,
        end=args.end,
        stream_only=not args.download_full,
        chr_prefix=not args.no_chr_prefix,
    )

    success = downloader.run(
        genome_vcf_checksum=args.genome_vcf_checksum,
        genome_tbi_checksum=args.genome_tbi_checksum,
        exome_vcf_checksum=args.exome_vcf_checksum,
        exome_tbi_checksum=args.exome_tbi_checksum
    )
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
