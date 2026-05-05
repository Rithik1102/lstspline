import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lstspline import CubicSpline, TensorSpline


PLOTS_DIR = PROJECT_ROOT / "plots" / "outputs"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)


def safe_spline_values(spline, x_grid, y_min=None, y_max=None):
    values = []

    for x in x_grid:
        try:
            y = spline.value(float(x))

            if not np.isfinite(y):
                values.append(np.nan)
                continue

            if y_min is not None and y < y_min:
                values.append(np.nan)
                continue

            if y_max is not None and y > y_max:
                values.append(np.nan)
                continue

            values.append(y)

        except Exception:
            values.append(np.nan)

    return np.array(values)


def plot_1d_dataset(name, x, y):
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)

    x_grid = np.linspace(x.min(), x.max(), 400)

    y_range = y.max() - y.min()
    y_min = y.min() - 2 * y_range - 1
    y_max = y.max() + 2 * y_range + 1

    spline_types = {
        1: "Cubic spline",
        3: "Monotone spline",
    }

    plt.figure(figsize=(9, 5))
    plt.scatter(x, y, label="Input data points")

    for spline_type, label in spline_types.items():
        spline = CubicSpline(x, y, type=spline_type)
        y_grid = safe_spline_values(spline, x_grid, y_min=y_min, y_max=y_max)
        plt.plot(x_grid, y_grid, label=label)

    plt.title(f"1D Spline Approximation: {name}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.grid(True)

    filename = name.lower().replace(" ", "_").replace("/", "_")
    path = PLOTS_DIR / f"1d_{filename}.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {path}")


def plot_all_1d():
    datasets = []

    datasets.append((
        "Plateau then sharp increase",
        np.arange(0, 10),
        [10, 10, 10, 10, 11, 25, 27, 28, 40, 50],
    ))

    datasets.append((
        "Irregular monotone increase",
        [0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15],
        [10, 10, 10, 10, 10.5, 15, 50, 60, 85, 90, 95],
    ))

    datasets.append((
        "Slow monotone growth",
        [0, 1, 2, 3, 4, 5, 6, 7],
        [0, 1, 2, 2.2, 2.2, 3, 5, 8],
    ))

    y_values = [
        .644, .622, .638, .649, .652, .639, .646, .657,
        .652, .655, .644, .663, .663, .668, .676, .676,
        .686, .679, .678, .683, .694, .699, .710, .730,
        .763, .812, .907, 1.044, 1.336, 1.881, 2.169, 2.075
    ]
    x_values = np.arange(1, len(y_values) + 1)
    datasets.append((
        "Provided y-only dataset",
        x_values,
        y_values,
    ))

    np.random.seed(0)
    x = np.linspace(0, 20, 80)
    y = 1 / (1 + np.exp(-(x - 10)))
    y += np.random.normal(0, 0.02, size=len(x))
    datasets.append((
        "Noisy sigmoid data",
        x,
        y,
    ))

    np.random.seed(1)
    x = np.linspace(0, 30, 120)
    y = np.piecewise(
        x,
        [x < 8, (x >= 8) & (x < 18), x >= 18],
        [
            lambda x: 0.05 * x,
            lambda x: 0.4 + 0.15 * (x - 8),
            lambda x: 2.0 + 0.02 * (x - 18) ** 2,
        ],
    )
    y += np.random.normal(0, 0.03, size=len(x))
    datasets.append((
        "Noisy piecewise growth",
        x,
        y,
    ))

    np.random.seed(2)
    x = np.linspace(0, 50, 150)
    y = np.floor(x / 5)
    y = y + np.random.normal(0, 0.4, size=len(x))
    datasets.append((
        "Noisy step data",
        x,
        y,
    ))

    for name, x, y in datasets:
        plot_1d_dataset(name, x, y)


def tensor_function(x, y):
    if x * y > 0:
        return abs(x) * y
    return 0.0


def plot_tensor_surface(noisy=False):
    np.random.seed(3 if noisy else 4)

    x_samples = np.linspace(-1, 1, 6)
    y_samples = np.linspace(-1, 1, 6)

    data = []
    for x in x_samples:
        for y in y_samples:
            z = tensor_function(x, y)
            if noisy:
                z += np.random.normal(0, 0.05)
            data.append([float(x), float(y), float(z)])

    knots = [
        [-1.0, -0.5, 0.0, 0.5, 1.0],
        [-1.0, -0.5, 0.0, 0.5, 1.0],
    ]

    exactdata = [
        [-1.0, -1.0, tensor_function(-1.0, -1.0)],
        [0.0, 0.0, tensor_function(0.0, 0.0)],
        [1.0, 1.0, tensor_function(1.0, 1.0)],
    ]

    spline = TensorSpline(
        kind=0,
        dim=2,
        knots=knots,
        data=data,
        exactdata=exactdata,
    )

    grid = np.linspace(-1, 1, 50)
    X, Y = np.meshgrid(grid, grid)
    Z = np.zeros_like(X)

    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            try:
                Z[i, j] = spline.value([float(X[i, j]), float(Y[i, j])])
            except Exception:
                Z[i, j] = np.nan

    data_np = np.array(data)

    suffix = "noisy" if noisy else "clean"

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, Z, alpha=0.75)
    ax.scatter(data_np[:, 0], data_np[:, 1], data_np[:, 2], s=35, label="Input data points")
    ax.set_title(f"3D Tensor Spline Surface Approximation ({suffix})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("f(x, y)")
    ax.legend()

    path = PLOTS_DIR / f"3d_tensor_surface_{suffix}.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2)
    ax.scatter(data_np[:, 0], data_np[:, 1], data_np[:, 2], s=35, label="Input data points")
    ax.set_title(f"3D Tensor Spline Wireframe Approximation ({suffix})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("f(x, y)")
    ax.legend()

    path = PLOTS_DIR / f"3d_tensor_wireframe_{suffix}.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")


def plot_all_3d():
    plot_tensor_surface(noisy=False)
    plot_tensor_surface(noisy=True)


def main():
    print("Generating 1D spline plots...")
    plot_all_1d()

    print("Generating 3D tensor spline plots...")
    plot_all_3d()

    print("All plots generated in:", PLOTS_DIR)


if __name__ == "__main__":
    main()