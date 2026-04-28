import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lstspline import CubicSpline


def safe_values(spline, x_grid, y_min=-2, y_max=2):
    values = []
    for xx in x_grid:
        try:
            yy = spline.value(float(xx))
            if not np.isfinite(yy) or yy < y_min or yy > y_max:
                values.append(np.nan)
            else:
                values.append(yy)
        except Exception:
            values.append(np.nan)
    return np.array(values)


def plot_comparison(x, y, dataset_name):
    x_grid = np.linspace(float(x.min()), float(x.max()), 300)

    spline_types = {
        1: "Cubic spline",
        3: "Monotone spline",
    }

    plt.figure(figsize=(9, 6))
    plt.scatter(x, y, label="Input data points")

    for spline_type, label in spline_types.items():
        spline = CubicSpline(x, y, type=spline_type)
        y_grid = safe_values(spline, x_grid, y_min=-1, y_max=2)
        plt.plot(x_grid, y_grid, label=label)

    plt.title(f"2D Spline Approximation: {dataset_name}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.grid(True)

    output_path = PROJECT_ROOT / "plots" / f"safe_2d_{dataset_name.lower().replace(' ', '_')}.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()

    print(f"Saved: {output_path}")


def main():
    x_mono = np.array([0, 1, 2, 3, 4, 5], dtype=float)
    y_mono = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5], dtype=float)

    x_curve = np.array([0, 1, 2, 3, 4, 5], dtype=float)
    y_curve = np.array([0.0, 0.2, 0.7, 0.5, 0.9, 0.6], dtype=float)

    plot_comparison(x_mono, y_mono, "Monotone linear data")
    plot_comparison(x_curve, y_curve, "Non-monotone curved data")


if __name__ == "__main__":
    main()