import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_PATH = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC_PATH))

from lstspline import CubicSpline

x = [0,1,2,3,4,5]
y = [0,0.1,0.2,0.3,0.4,0.5]

#x = [0,1,2,3]
#y = [0,1,0,1]

s = CubicSpline(x, y)

print(s.value(2.5))