"""Constants for data ingestion and processing."""

from pathlib import Path

# ClinVar Configuration
CLINVAR_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar"
CLINVAR_RELEASE_DATE = "20251116"
CLINVAR_EXPECTED_VCF_SIZE = 173000000  # ~173MB
CLINVAR_VCF_FILENAME = f"clinvar_{CLINVAR_RELEASE_DATE}.vcf.gz"
CLINVAR_TSV_FILENAME = "variant_summary.txt.gz"

# gnomAD Configuration
GNOMAD_VERSION = "4.1"
GNOMAD_BASE_URL = "https://storage.googleapis.com/download/storage/v1/b/gcp-public-data--gnomad/o"
# ABCA4 Region: chr1:93,500,000-95,000,000 (±500kb from gene)
GNOMAD_DEFAULT_CHROM = "1"
GNOMAD_DEFAULT_START = 93500000
GNOMAD_DEFAULT_END = 95000000

# Gene Configuration
ABCA4_GENE_SYMBOL = "ABCA4"
