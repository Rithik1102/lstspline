import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from sklearn.isotonic import IsotonicRegression
from scipy.interpolate import CubicSpline as SciPyCubicSpline
from scipy.interpolate import PchipInterpolator

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lstspline import CubicSpline

OUT_DIR = PROJECT_ROOT / "plots" / "midterm_2d_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def safe_values(spline, x_grid):
    values = []
    for t in x_grid:
        try:
            y = spline.value(float(t))
            values.append(y if np.isfinite(y) else np.nan)
        except Exception:
            values.append(np.nan)
    return np.array(values)


def clipped(values, lower=-200, upper=200):
    return np.clip(values, lower, upper)


def monotonicity_violations(y_values):
    dy = np.diff(y_values)
    return int(np.sum(dy < -1e-8))


def plot_interpolation_dataset(name, x, y, zoom=None):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x_grid = np.linspace(x.min(), x.max(), 700)

    our_cubic = CubicSpline(x, y, type=1)
    our_taut = CubicSpline(x, y, type=2, tau=10.0)
    our_mono = CubicSpline(x, y, type=3)

    y_our_cubic = safe_values(our_cubic, x_grid)
    y_our_taut = clipped(safe_values(our_taut, x_grid))
    y_our_mono = safe_values(our_mono, x_grid)

    scipy_cubic = SciPyCubicSpline(x, y)
    scipy_pchip = PchipInterpolator(x, y)

    plt.figure(figsize=(10, 6))
    plt.scatter(x, y, label="Input data", s=35)

    plt.plot(x_grid, y_our_cubic, label="Our cubic spline")
    plt.plot(x_grid, y_our_taut, label="Our taut spline (tau=10.0)")
    plt.plot(x_grid, y_our_mono, label="Our monotone spline")
    plt.plot(x_grid, scipy_cubic(x_grid), "--", label="SciPy cubic spline")
    plt.plot(x_grid, scipy_pchip(x_grid), "--", label="SciPy PCHIP")

    if zoom:
        plt.xlim(zoom[0], zoom[1])

    plt.title(name)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()

    fname = name.lower().replace(" ", "_").replace("/", "_")
    if zoom:
        fname += "_zoom"

    path = OUT_DIR / f"{fname}.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()

    print("Saved:", path)
    print(name)
    print("  Cubic violations:", monotonicity_violations(y_our_cubic))
    print("  Taut violations:", monotonicity_violations(y_our_taut))
    print("  Monotone violations:", monotonicity_violations(y_our_mono))


def plot_noisy_with_isotonic(name, x, y, true_y=None):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    ir = IsotonicRegression(increasing=True)
    y_iso = ir.fit_transform(x, y)

    x_grid = np.linspace(x.min(), x.max(), 700)

    cubic_raw = CubicSpline(x, y, type=1)
    taut_raw = CubicSpline(x, y, type=2, tau=10.0)
    mono_iso = CubicSpline(x, y_iso, type=3)

    y_cubic_raw = safe_values(cubic_raw, x_grid)
    y_taut_raw = clipped(safe_values(taut_raw, x_grid))
    y_mono_iso = safe_values(mono_iso, x_grid)

    pchip_iso = PchipInterpolator(x, y_iso)

    plt.figure(figsize=(10, 6))
    plt.scatter(x, y, label="Original noisy data", s=25)
    plt.plot(x, y_iso, marker="o", linewidth=1.5, label="Isotonic regression / PAVA data")
    plt.plot(x_grid, y_cubic_raw, label="Cubic spline on noisy data")
    plt.plot(x_grid, y_taut_raw, label="Taut spline on noisy data (tau=10.0)")
    plt.plot(x_grid, y_mono_iso, label="Monotone spline after isotonic regression")
    plt.plot(x_grid, pchip_iso(x_grid), "--", label="SciPy PCHIP after isotonic regression")

    if true_y is not None:
        true_grid = np.interp(x_grid, x, true_y)
        plt.plot(x_grid, true_grid, ":", label="True underlying function")

    plt.title(name)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()

    fname = name.lower().replace(" ", "_").replace("/", "_")
    path = OUT_DIR / f"isotonic_{fname}.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()

    print("Saved:", path)
    print(name)
    print("  Cubic raw violations:", monotonicity_violations(y_cubic_raw))
    print("  Taut raw violations:", monotonicity_violations(y_taut_raw))
    print("  Monotone after isotonic violations:", monotonicity_violations(y_mono_iso))


def main():
    x1 = np.arange(0, 10)
    y1 = [10, 10, 10, 10, 11, 25, 27, 28, 40, 50]
    plot_interpolation_dataset("Plateau then sharp increase", x1, y1)
    plot_interpolation_dataset("Plateau then sharp increase", x1, y1, zoom=(3, 7))

    x2 = [0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15]
    y2 = [10, 10, 10, 10, 10.5, 15, 50, 60, 85, 90, 95]
    plot_interpolation_dataset("Irregular monotone increase", x2, y2)
    plot_interpolation_dataset("Irregular monotone increase", x2, y2, zoom=(7, 10))

    x3 = [0, 1, 2, 3, 4, 5, 6, 7]
    y3 = [0, 1, 2, 2.2, 2.2, 3, 5, 8]
    plot_interpolation_dataset("Slow monotone growth", x3, y3)
    plot_interpolation_dataset("Slow monotone growth", x3, y3, zoom=(2, 5))

    np.random.seed(0)
    x = np.linspace(0, 20, 80)
    y_true = 1 / (1 + np.exp(-(x - 10)))
    y_noisy = y_true + np.random.normal(0, 0.02, size=len(x))
    plot_noisy_with_isotonic("Noisy sigmoid data", x, y_noisy, y_true)

    np.random.seed(1)
    x = np.linspace(0, 30, 120)
    y_true = np.piecewise(
        x,
        [x < 8, (x >= 8) & (x < 18), x >= 18],
        [
            lambda x: 0.05 * x,
            lambda x: 0.4 + 0.15 * (x - 8),
            lambda x: 2.0 + 0.02 * (x - 18) ** 2,
        ],
    )
    y_noisy = y_true + np.random.normal(0, 0.03, size=len(x))
    plot_noisy_with_isotonic("Noisy piecewise growth", x, y_noisy, y_true)

    np.random.seed(2)
    x = np.linspace(0, 50, 150)
    y_true = np.floor(x / 5)
    y_noisy = y_true + np.random.normal(0, 0.4, size=len(x))
    plot_noisy_with_isotonic("Noisy step data", x, y_noisy, y_true)

    print("Done. 2D plots saved in:", OUT_DIR)


if __name__ == "__main__":
    main()