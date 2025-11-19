"""Constants for conservation feature calculation."""

# UCSC Genome Browser Config
UCSC_BASE_URL = "https://api.genome.ucsc.edu"
UCSC_TRACKS = ["phyloP100way", "phastCons100way"]
UCSC_CHROM = "chr1"
UCSC_CHUNK_SIZE = 1000

# Conservation Weights
PHYLOP_WEIGHT = 0.6
PHASTCONS_WEIGHT = 0.4

# Thresholds
PHYLOP_HIGH_THRESHOLD = 2.0
PHASTCONS_HIGH_THRESHOLD = 0.8
