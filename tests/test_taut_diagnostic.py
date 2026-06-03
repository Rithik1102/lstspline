import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lstspline import CubicSpline


def test_dataset(name, x, y):
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)

    x_grid = np.linspace(x.min(), x.max(), 300)

    plt.figure(figsize=(9, 5))
    plt.scatter(x, y, label="Input data")

    for tau in [0.01, 0.05, 0.1, 0.5, 1.0]:
        try:
            s = CubicSpline(x, y, type=2, tau=tau)
            y_grid = np.array([s.value(float(xx)) for xx in x_grid])

            print(name, "tau", tau)
            print("min:", np.nanmin(y_grid), "max:", np.nanmax(y_grid))
            print("value at middle:", s.value(float(x_grid[len(x_grid)//2])))

            plt.plot(x_grid, y_grid, label=f"taut tau={tau}")

        except Exception as e:
            print(name, "tau", tau, "FAILED:", e)

    plt.title(name)
    plt.legend()
    plt.grid(True)
    plt.show()


# Very simple monotone linear data
test_dataset(
    "Linear monotone data",
    np.array([0, 1, 2, 3, 4, 5], dtype=float),
    np.array([0, 1, 2, 3, 4, 5], dtype=float),
)

# Smooth nonlinear data
test_dataset(
    "Smooth quadratic-like data",
    np.array([0, 1, 2, 3, 4, 5], dtype=float),
    np.array([0, 1, 1.5, 2.2, 3.4, 5.0], dtype=float),
)

# Professor-like data
test_dataset(
    "Plateau then sharp increase",
    np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=float),
    np.array([10, 10, 10, 10, 11, 25, 27, 28, 40, 50], dtype=float),
)