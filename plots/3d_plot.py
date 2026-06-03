import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lstspline import TensorSpline

OUT_DIR = PROJECT_ROOT / "plots" / "midterm_3d_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def target_function(x, y):
    return abs(x) * y if x * y > 0 else 0.0


def build_tensor_data(noisy=False):
    np.random.seed(10 if not noisy else 11)

    xs = np.linspace(-1, 1, 7)
    ys = np.linspace(-1, 1, 7)

    data = []
    for x in xs:
        for y in ys:
            z = target_function(x, y)
            if noisy:
                z += np.random.normal(0, 0.05)
            data.append([float(x), float(y), float(z)])

    exactdata = [
        [-1.0, -1.0, target_function(-1.0, -1.0)],
        [0.0, 0.0, target_function(0.0, 0.0)],
        [1.0, 1.0, target_function(1.0, 1.0)],
    ]

    knots = [
        [-1.0, -0.5, 0.0, 0.5, 1.0],
        [-1.0, -0.5, 0.0, 0.5, 1.0],
    ]

    return knots, data, exactdata


def evaluate_tensor_spline(spline):
    grid = np.linspace(-1, 1, 70)
    X, Y = np.meshgrid(grid, grid)
    Z = np.zeros_like(X)

    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            try:
                Z[i, j] = spline.value([float(X[i, j]), float(Y[i, j])])
            except Exception:
                Z[i, j] = np.nan

    return X, Y, Z


def plot_surface(X, Y, Z, data, title, filename):
    data = np.asarray(data)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot_surface(X, Y, Z, alpha=0.75)
    ax.scatter(data[:, 0], data[:, 1], data[:, 2], s=35, label="Input data points")

    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("f(x,y)")
    ax.legend()

    path = OUT_DIR / filename
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()

    print("Saved:", path)


def plot_wireframe(X, Y, Z, data, title, filename):
    data = np.asarray(data)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot_wireframe(X, Y, Z, rstride=3, cstride=3)
    ax.scatter(data[:, 0], data[:, 1], data[:, 2], s=35, label="Input data points")

    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("f(x,y)")
    ax.legend()

    path = OUT_DIR / filename
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()

    print("Saved:", path)


def run_tensor_experiment(noisy=False):
    suffix = "noisy" if noisy else "clean"

    knots, data, exactdata = build_tensor_data(noisy=noisy)

    spline = TensorSpline(
        kind=0,
        dim=2,
        knots=knots,
        data=data,
        exactdata=exactdata,
    )

    X, Y, Z = evaluate_tensor_spline(spline)

    plot_surface(
        X,
        Y,
        Z,
        data,
        f"Tensor-product spline surface ({suffix})",
        f"tensor_surface_{suffix}.png",
    )

    plot_wireframe(
        X,
        Y,
        Z,
        data,
        f"Tensor-product spline wireframe ({suffix})",
        f"tensor_wireframe_{suffix}.png",
    )


def main():
    run_tensor_experiment(noisy=False)
    run_tensor_experiment(noisy=True)
    print("Done. 3D plots saved in:", OUT_DIR)


if __name__ == "__main__":
    main()