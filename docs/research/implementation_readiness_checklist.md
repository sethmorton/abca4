# ABCA4 Campaign - Implementation Readiness Checklist

**Date:** November 14, 2025
**Status:** ✅ **Documentation Finalized - Ready for Implementation**

## 🎯 Executive Summary

All documentation, tooling decisions, and architectural choices have been finalized for the ABCA4 campaign implementation. We are adopting a **modern Python development stack** centered on `uv` + Marimo for reproducible, interactive genomics research.

**Key Decisions Finalized:**
- **Package Manager**: `uv` (10-100x faster than pip, reproducible environments)
- **Interactive Development**: Marimo (reactive notebooks, SQL support, git-friendly)
- **Orchestration**: Hybrid invoke + Marimo workflow
- **Data Sources**: ClinVar, gnomAD, SpliceAI/AlphaMissense verified and accessible
- **Package Versions**: All pinned and compatible (November 2025)

---

## 📋 Pre-Implementation Checklist

### ✅ Environment & Tooling
- [x] `uv` adopted as primary package manager across all documentation
- [x] Marimo integrated for interactive development
- [x] Environment setup documented (`docs/research/abca4_campaign/env_notes.md`)
- [x] All README, tutorials, and guides updated to use `uv` commands

### ✅ Data Sources Verified
- [x] ClinVar GRCh38 VCF: `clinvar_20251109.vcf.gz` (181MB)
- [x] gnomAD v4.1.0: Genome/exome VCFs with ABCA4 region extraction
- [x] SpliceAI scores: Public BigQuery/TSV dumps available
- [x] AlphaMissense: Public database accessible
- [x] Protein domains: UniProt/ClinGen sources identified
- [x] DNA FM weights: Enformer, HyenaDNA, Evo2 download paths verified

### ✅ Package Ecosystem Finalized
| Component | Package | Version | Status |
|-----------|---------|---------|--------|
| Package Manager | uv | 0.9.9 | ✅ Adopted |
| Genomics IO | cyvcf2, pysam, biopython | 0.31.4, 0.23.3, 1.86 | ✅ Verified |
| Interactive Dev | marimo | 0.17.8 | ✅ Adopted |
| ML Framework | torch, transformers | 2.9.1, 4.35.0 | ✅ Compatible |
| Orchestration | invoke, mlflow | 2.2.1, 3.6.0 | ✅ Integrated |

### ✅ Architecture Decisions
- [x] **Hybrid Workflow**: invoke (automation) + Marimo (exploration)
- [x] **Storage Strategy**: Parquet for features, git-friendly notebooks
- [x] **Caching**: PyArrow/Parquet with MLflow provenance tracking
- [x] **Interactive Notebooks**: 3 concrete examples ready (`notebooks/01_*.py`)
- [x] **Environment Management**: `uv venv` + `uv pip sync` for reproducibility

### ✅ Implementation Assets Ready
- [x] **Marimo Notebooks**: 3 complete templates with reactive UIs
- [x] **Data Download Plans**: wget/gsutil commands verified
- [x] **Feature Engineering**: Modular scripts outlined
- [x] **Orchestration**: invoke task stubs defined in plan
- [x] **Documentation**: Complete setup guides and tutorials

---

## 🚀 Implementation Plan (Ready to Execute)

### Phase 1: Environment & Scaffolding (Week 1)
```bash
# Install uv and create environment
curl -LsSf https://astral.sh/uv/install.sh | sh
uv venv
uv pip sync requirements-dev.txt pyproject.toml
uv pip install -e .[interactive]

# Test Marimo integration
uv run marimo edit notebooks/01_data_exploration.py
```

### Phase 2: Data Ingestion (Week 1-2)
- Implement `src/data/download_clinvar.py`
- Implement `src/data/download_gnomad.py`
- Create `src/data/filter_abca4_variants.py`
- Test with sample data before full downloads

### Phase 3: Feature Computation (Week 2-3)
- Build modular feature scripts (`src/features/*.py`)
- Integrate with Marimo notebook for interactive development
- Test feature correlations and distributions

### Phase 4: Strand Integration (Week 3-4)
- Implement ABCA4-specific reward environment
- Create optimization configurations
- Build interactive dashboards

### Phase 5: Reporting & Validation (Week 4-5)
- Generate variant intelligence snapshots
- Create stakeholder dashboards
- Validate end-to-end pipeline

---

## 🛠️ Key Commands (Ready to Use)

### Environment Setup
```bash
# Complete setup (one-time)
curl -LsSf https://astral.sh/uv/install.sh | sh
uv venv
uv pip sync requirements-dev.txt pyproject.toml
uv pip install -e .[interactive]
```

### Development Workflow
```bash
# Activate environment
source .venv/bin/activate

# Edit notebooks
uv run marimo edit notebooks/01_data_exploration.py

# Run as web apps
uv run marimo run notebooks/01_data_exploration.py

# Execute as scripts
uv run python notebooks/01_data_exploration.py
```

### Package Management
```bash
# Add new dependency
uv pip install new-package
uv pip sync requirements-dev.txt  # Update lockfile

# Install extras
uv pip install -e .[variant-triage]
uv pip install -e .[interactive]
```

---

## 🔍 Quality Assurance

### Documentation Consistency
- [x] All `pip install` references converted to `uv pip install`
- [x] Environment setup instructions standardized
- [x] Marimo integration documented across notebooks and guides

### Technical Verification
- [x] Package versions verified current (November 2025)
- [x] Data source URLs tested and accessible
- [x] Import paths validated in existing codebase
- [x] GPU/CUDA compatibility documented

### Implementation Readiness
- [x] Concrete code examples provided (Marimo notebooks)
- [x] File/directory structure defined
- [x] Integration points identified (invoke + Marimo)
- [x] Error handling and logging strategies outlined

---

## 🎯 Success Criteria Met

✅ **Modern Tooling**: `uv` + Marimo stack adopted for speed and interactivity
✅ **Reproducible Environments**: Lockfile-based dependency management
✅ **Interactive Development**: Reactive notebooks for genomics exploration
✅ **Comprehensive Documentation**: All guides updated and consistent
✅ **Implementation Assets**: Concrete examples and clear next steps
✅ **Data Access Verified**: All required sources accessible and documented

## 🚦 Go/No-Go Decision

**✅ GO for Implementation**

All documentation finalized, tooling decisions made, and concrete examples provided. The ABCA4 campaign is ready to proceed with implementation using the modern `uv` + Marimo development stack.

**Next Action**: Begin Phase 1 implementation or dive into any specific component ready to build!

---

*This checklist ensures we start implementation with complete clarity on tooling, architecture, and requirements. No more "figure it out as we go" - everything is specified and ready.*
