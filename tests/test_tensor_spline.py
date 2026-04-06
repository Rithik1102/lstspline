import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_PATH = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC_PATH))

from lstspline import TensorSpline

knots = [
    [0.0, 0.5, 1.0],
    [0.0, 0.5, 1.0],
]

data = [
    [0.2, 0.3, 0.06],
    [0.7, 0.4, 0.28],
    [0.5, 0.5, 0.25],
    [0.9, 0.8, 0.72],
]

exactdata = [
    [0.0, 0.0, 0.0],
    [1.0, 1.0, 1.0],
]

s1 = TensorSpline(kind=0, dim=2, knots=knots, data=data, exactdata=exactdata)
s2 = TensorSpline(kind=0, dim=2, knots=knots, data=data, exactdata=exactdata)

print("s1 value:", s1.value([0.5, 0.5]))
print("s2 value:", s2.value([0.25, 0.75]))
print("s1 d/dx:", s1.value_der([0.5, 0.5], 0))
print("s1 d/dy:", s1.value_der([0.5, 0.5], 1))