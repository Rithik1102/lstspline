import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lstspline import TensorSpline


def target_function(x, y):
    # Simple monotone surface: f(x, y) = x * y
    return x * y


def main():
    plots_dir = PROJECT_ROOT / "plots"
    plots_dir.mkdir(exist_ok=True)

    # Knots for 2D tensor spline
    knots = [
        [0.0, 0.5, 1.0],
        [0.0, 0.5, 1.0],
    ]

    # Training data: [x, y, z]
    data = [
        [0.2, 0.3, target_function(0.2, 0.3)],
        [0.7, 0.4, target_function(0.7, 0.4)],
        [0.5, 0.5, target_function(0.5, 0.5)],
        [0.9, 0.8, target_function(0.9, 0.8)],
        [0.3, 0.9, target_function(0.3, 0.9)],
        [0.8, 0.2, target_function(0.8, 0.2)],
    ]

    # Exact boundary/interpolation points
    exactdata = [
        [0.0, 0.0, target_function(0.0, 0.0)],
        [1.0, 1.0, target_function(1.0, 1.0)],
    ]

    spline = TensorSpline(
        kind=0,
        dim=2,
        knots=knots,
        data=data,
        exactdata=exactdata,
    )

    # Grid for surface
    x_grid = np.linspace(0.0, 1.0, 40)
    y_grid = np.linspace(0.0, 1.0, 40)

    X, Y = np.meshgrid(x_grid, y_grid)
    Z = np.zeros_like(X)

    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            Z[i, j] = spline.value([float(X[i, j]), float(Y[i, j])])

    data_np = np.array(data)

    # 1. Surface plot
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot_surface(X, Y, Z, alpha=0.75)
    ax.scatter(data_np[:, 0], data_np[:, 1], data_np[:, 2], s=50, label="Input data points")

    ax.set_title("3D Tensor Spline Surface Approximation")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("f(x, y)")
    ax.legend()

    surface_path = plots_dir / "tensor_spline_3d_surface.png"
    plt.savefig(surface_path, dpi=300, bbox_inches="tight")
    plt.show()

    print(f"Saved: {surface_path}")

    # 2. Wireframe plot
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2)
    ax.scatter(data_np[:, 0], data_np[:, 1], data_np[:, 2], s=50, label="Input data points")

    ax.set_title("3D Tensor Spline Wireframe Approximation")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("f(x, y)")
    ax.legend()

    wireframe_path = plots_dir / "tensor_spline_3d_wireframe.png"
    plt.savefig(wireframe_path, dpi=300, bbox_inches="tight")
    plt.show()

    print(f"Saved: {wireframe_path}")


if __name__ == "__main__":
    main()