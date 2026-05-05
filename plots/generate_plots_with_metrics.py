import sys
from pathlib import Path
import csv

import numpy as np
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lstspline import CubicSpline

try:
    from scipy.interpolate import CubicSpline as SciPyCubicSpline
    from scipy.interpolate import PchipInterpolator
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


PLOTS_DIR = PROJECT_ROOT / "plots" / "metrics_outputs"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

METRICS_FILE = PLOTS_DIR / "mse_comparison_results.csv"


def safe_values(model, x_grid):
    values = []
    for xx in x_grid:
        try:
            yy = model.value(float(xx))
            if np.isfinite(yy):
                values.append(yy)
            else:
                values.append(np.nan)
        except Exception:
            values.append(np.nan)
    return np.array(values)


def mse(y_true, y_pred):
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    if np.sum(mask) == 0:
        return np.nan
    return np.mean((y_true[mask] - y_pred[mask]) ** 2)


def mae(y_true, y_pred):
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    if np.sum(mask) == 0:
        return np.nan
    return np.mean(np.abs(y_true[mask] - y_pred[mask]))


def monotonicity_violations(y_pred):
    diffs = np.diff(y_pred)
    diffs = diffs[np.isfinite(diffs)]
    return int(np.sum(diffs < -1e-6))


def plot_and_measure_dataset(name, x, y_noisy, y_true_func=None):
    x = np.array(x, dtype=float)
    y_noisy = np.array(y_noisy, dtype=float)

    x_grid = np.linspace(x.min(), x.max(), 500)

    if y_true_func is not None:
        y_true = y_true_func(x_grid)
    else:
        y_true = None

    results = []

    plt.figure(figsize=(10, 6))
    plt.scatter(x, y_noisy, label="Input data points", s=30)

    # cubic spline
    our_cubic = CubicSpline(x, y_noisy, type=1)
    y_our_cubic = safe_values(our_cubic, x_grid)
    plt.plot(x_grid, y_our_cubic, label="Our cubic spline")

    # monotone spline
    our_mono = CubicSpline(x, y_noisy, type=3)
    y_our_mono = safe_values(our_mono, x_grid)
    plt.plot(x_grid, y_our_mono, label="Our monotone spline")

    if y_true is not None:
        plt.plot(x_grid, y_true, linestyle="--", label="True function")

    results.append({
        "dataset": name,
        "method": "Our cubic spline",
        "mse": mse(y_true, y_our_cubic) if y_true is not None else "",
        "mae": mae(y_true, y_our_cubic) if y_true is not None else "",
        "monotonicity_violations": monotonicity_violations(y_our_cubic),
    })

    results.append({
        "dataset": name,
        "method": "Our monotone spline",
        "mse": mse(y_true, y_our_mono) if y_true is not None else "",
        "mae": mae(y_true, y_our_mono) if y_true is not None else "",
        "monotonicity_violations": monotonicity_violations(y_our_mono),
    })

    # SciPy comparison
    if SCIPY_AVAILABLE:
        scipy_cubic = SciPyCubicSpline(x, y_noisy)
        y_scipy_cubic = scipy_cubic(x_grid)
        plt.plot(x_grid, y_scipy_cubic, label="SciPy cubic spline")

        scipy_pchip = PchipInterpolator(x, y_noisy)
        y_pchip = scipy_pchip(x_grid)
        plt.plot(x_grid, y_pchip, label="SciPy PCHIP monotone")

        results.append({
            "dataset": name,
            "method": "SciPy cubic spline",
            "mse": mse(y_true, y_scipy_cubic) if y_true is not None else "",
            "mae": mae(y_true, y_scipy_cubic) if y_true is not None else "",
            "monotonicity_violations": monotonicity_violations(y_scipy_cubic),
        })

        results.append({
            "dataset": name,
            "method": "SciPy PCHIP monotone",
            "mse": mse(y_true, y_pchip) if y_true is not None else "",
            "mae": mae(y_true, y_pchip) if y_true is not None else "",
            "monotonicity_violations": monotonicity_violations(y_pchip),
        })

    plt.title(f"Spline Comparison and Benchmarking: {name}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.grid(True)

    safe_name = name.lower().replace(" ", "_").replace("/", "_")
    plot_path = PLOTS_DIR / f"{safe_name}_comparison.png"
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved plot: {plot_path}")

    return results


def main():
    all_results = []

    # Dataset 1: noisy sigmoid with known true function
    np.random.seed(0)
    x = np.linspace(0, 20, 80)

    def sigmoid_true(x_grid):
        return 1 / (1 + np.exp(-(x_grid - 10)))

    y_true = sigmoid_true(x)
    y_noisy = y_true + np.random.normal(0, 0.02, size=len(x))

    all_results.extend(
        plot_and_measure_dataset(
            "Noisy sigmoid data",
            x,
            y_noisy,
            y_true_func=sigmoid_true,
        )
    )

    # Dataset 2: noisy piecewise growth with known true function
    np.random.seed(1)
    x = np.linspace(0, 30, 120)

    def piecewise_true(x_grid):
        return np.piecewise(
            x_grid,
            [x_grid < 8, (x_grid >= 8) & (x_grid < 18), x_grid >= 18],
            [
                lambda x: 0.05 * x,
                lambda x: 0.4 + 0.15 * (x - 8),
                lambda x: 2.0 + 0.02 * (x - 18) ** 2,
            ],
        )

    y_true = piecewise_true(x)
    y_noisy = y_true + np.random.normal(0, 0.03, size=len(x))

    all_results.extend(
        plot_and_measure_dataset(
            "Noisy piecewise growth",
            x,
            y_noisy,
            y_true_func=piecewise_true,
        )
    )

    # Dataset 3: noisy step data with known true function
    np.random.seed(2)
    x = np.linspace(0, 50, 150)

    def step_true(x_grid):
        return np.floor(x_grid / 5)

    y_true = step_true(x)
    y_noisy = y_true + np.random.normal(0, 0.4, size=len(x))

    all_results.extend(
        plot_and_measure_dataset(
            "Noisy step data",
            x,
            y_noisy,
            y_true_func=step_true,
        )
    )

    # Dataset 4: plateau dataset, no known true function
    x = np.arange(0, 10)
    y = np.array([10, 10, 10, 10, 11, 25, 27, 28, 40, 50], dtype=float)

    all_results.extend(
        plot_and_measure_dataset(
            "Plateau then sharp increase",
            x,
            y,
            y_true_func=None,
        )
    )

    # Save CSV metrics
    with open(METRICS_FILE, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "dataset",
                "method",
                "mse",
                "mae",
                "monotonicity_violations",
            ],
        )
        writer.writeheader()
        writer.writerows(all_results)

    print(f"Saved metrics: {METRICS_FILE}")

    if not SCIPY_AVAILABLE:
        print("\nSciPy is not installed.")
        print("To include SciPy benchmarking, run:")
        print("pip install scipy")


if __name__ == "__main__":
    main()