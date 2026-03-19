import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_PATH = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC_PATH))

from lstspline import MyArray

arr = MyArray([1, 2, 3, 4])

print("Value at 1:", arr.get(1))
print("sum:", arr.sum())