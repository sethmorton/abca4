# ABCA4 Campaign - Environment Setup Notes

**Date:** November 14, 2025
**Status:** Finalized for implementation

## uv Adoption Rationale

We're adopting `uv` as the primary package manager for the ABCA4 campaign due to:

- **10-100x faster** dependency resolution than pip
- **Reproducible environments** via lockfile-based dependency management
- **Integrated virtual environment** management
- **Modern Python tooling** aligned with the project's forward-looking approach
- **Seamless integration** with pyproject.toml and modern Python packaging

## Core Setup Commands

### Development Environment
```bash
# Install uv (one-time)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create virtual environment
uv venv

# Sync dependencies from requirements-dev.txt
uv pip sync requirements-dev.txt

# Install project in editable mode
uv pip install -e .

# Install interactive dependencies (Marimo)
uv pip install -e .[interactive]
```

### Production Environment
```bash
# Create isolated environment
uv venv --python 3.11

# Sync only production dependencies
uv pip sync requirements.txt pyproject.toml

# Install with specific extras
uv pip install -e .[variant-triage]
```

### GPU/CUDA Environment
```bash
# For PyTorch with CUDA support
uv venv
uv pip sync requirements-dev.txt
uv pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# Verify CUDA availability
uv run python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"
```

## uv Command Reference

| Command | Purpose | Example |
|---------|---------|---------|
| `uv venv` | Create virtual environment | `uv venv --python 3.11` |
| `uv pip sync` | Install exact versions from lockfile | `uv pip sync requirements-dev.txt` |
| `uv pip install` | Install packages | `uv pip install marimo` |
| `uv run` | Run commands in environment | `uv run python script.py` |
| `uv lock` | Update lockfile | `uv lock --upgrade` |
| `uv export` | Export requirements.txt | `uv export --format requirements-txt` |

## Environment Variables

### For GPU Workloads
```bash
export CUDA_VISIBLE_DEVICES=0,1  # Specify GPU devices
export TORCH_USE_CUDA_DSA=1      # Enable CUDA device-side assertions
```

### For Large Dataset Downloads
```bash
export GCS_REQUEST_TIMEOUT=300   # Longer timeout for gnomAD downloads
export HTTP_TIMEOUT=300          # General HTTP timeout
```

## Lockfile Management

### Creating Lockfiles
```bash
# Update lockfile with latest compatible versions
uv lock

# Update specific package
uv lock --upgrade-package torch

# Check for outdated packages
uv lock --dry-run
```

### Exporting for CI/CD
```bash
# Export to requirements.txt for legacy systems
uv export --format requirements-txt > requirements-ci.txt

# Export with hashes for security
uv export --format requirements-txt --hashes > requirements-secure.txt
```

## Troubleshooting

### Common Issues

1. **"uv command not found"**
   ```bash
   # Reinstall uv
   curl -LsSf https://astral.sh/uv/install.sh | sh
   source ~/.cargo/env  # or restart shell
   ```

2. **Environment not activated**
   ```bash
   # Use uv run instead of activating
   uv run python script.py

   # Or activate manually
   source .venv/bin/activate
   ```

3. **CUDA compatibility issues**
   ```bash
   # Check PyTorch CUDA version
   uv run python -c "import torch; print(torch.version.cuda)"

   # Reinstall with correct CUDA version
   uv pip install torch --index-url https://download.pytorch.org/whl/cu121
   ```

### Performance Tips

- **Use `uv pip sync`** instead of `pip install` for reproducible builds
- **Run `uv lock --upgrade`** monthly to stay current with security updates
- **Use `uv run`** for one-off commands to avoid environment activation
- **Cache uv's download directory** in CI for faster builds

## Package Version Pins

All packages are pinned in `requirements-dev.txt` and `pyproject.toml` with versions verified as of November 2025:

- **Core**: uv (0.9.9), Python 3.11+
- **Genomics**: cyvcf2 (0.31.4), pysam (0.23.3), biopython (1.86)
- **ML**: torch (2.9.1), transformers (4.35.0), accelerate (0.24.0)
- **Interactive**: marimo (0.17.8)
- **Orchestration**: invoke (2.2.1), mlflow (3.6.0)

## Migration from pip

If migrating an existing environment:

```bash
# Export current environment
pip freeze > requirements-old.txt

# Create uv environment
uv venv

# Install from old requirements (may need manual resolution)
uv pip install -r requirements-old.txt

# Sync with project requirements
uv pip sync requirements-dev.txt pyproject.toml
```

## Integration with CI/CD

### GitHub Actions Example
```yaml
- name: Setup uv
  run: curl -LsSf https://astral.sh/uv/install.sh | sh

- name: Setup Python environment
  run: uv venv --python 3.11

- name: Install dependencies
  run: uv pip sync requirements-dev.txt pyproject.toml

- name: Run tests
  run: uv run pytest
```

This setup provides **reproducible**, **fast**, and **modern** Python environment management for the ABCA4 campaign.
