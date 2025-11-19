#!/usr/bin/env python3
import sys
from pathlib import Path

print("Current directory:", Path.cwd())
print("Python path:")
for i, p in enumerate(sys.path[:5]):
    print(f"  {i}: {p}")

# Test the path setup like in the notebook
try:
    _file_path = Path(__file__).resolve()
    print(f"__file__ resolved: {_file_path}")
except NameError:
    print("No __file__ available")
    _file_path = Path.cwd() / "notebooks" / "01_data_exploration.py"
    print(f"Fallback path: {_file_path}")

_project_root = _file_path.parent.parent
print(f"Project root: {_project_root}")

if str(_project_root) not in sys.path:
    sys.path.insert(0, str(_project_root))
    print("Added project root to path")

print("Testing src import...")
try:
    import src
    print("✅ src module found")
except ImportError as e:
    print(f"❌ src module not found: {e}")

try:
    from src.config import logger
    print("✅ src.config imported successfully")
except ImportError as e:
    print(f"❌ src.config import failed: {e}")
