import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.interpolate import CubicSpline as SciPyCubicSpline
from scipy.interpolate import PchipInterpolator, UnivariateSpline
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lstspline import CubicSpline


# ============================================================
# Paths
# ============================================================

DATA_PATH = PROJECT_ROOT / "data" / "melb_data.csv"

OUT_DIR = PROJECT_ROOT / "benchmarks" / "housing_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PLOTS_DIR = OUT_DIR / "plots"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

RESULTS_CSV = OUT_DIR / "melbourne_housing_benchmark_results.csv"


# ============================================================
# Helper functions
# ============================================================

def load_housing_data():
    df = pd.read_csv(DATA_PATH)

    required_cols = ["Distance", "Price"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    df = df[required_cols].copy()
    df = df.dropna()

    df["Distance"] = pd.to_numeric(df["Distance"], errors="coerce")
    df["Price"] = pd.to_numeric(df["Price"], errors="coerce")
    df = df.dropna()

    # Remove unrealistic or zero values
    df = df[(df["Distance"] >= 0) & (df["Price"] > 0)]

    return df


def create_distance_price_curve(df, bins=35):
    """
    Housing data is noisy and has repeated distance values.
    For spline approximation, we aggregate median price by distance bins.
    This converts scattered housing data into a smoother regression curve.
    """

    df = df.copy()

    df["distance_bin"] = pd.cut(df["Distance"], bins=bins)

    curve = (
        df.groupby("distance_bin", observed=True)
        .agg(
            distance=("Distance", "mean"),
            median_price=("Price", "median"),
            count=("Price", "size")
        )
        .reset_index(drop=True)
    )

    # Keep bins with enough samples
    curve = curve[curve["count"] >= 10]
    curve = curve.dropna()
    curve = curve.sort_values("distance")

    x = curve["distance"].to_numpy(dtype=float)
    y = curve["median_price"].to_numpy(dtype=float)

    # Remove duplicate x values if any
    unique_x, unique_indices = np.unique(x, return_index=True)
    x = unique_x
    y = y[unique_indices]

    return curve, x, y


def normalise_x(x):
    return x - np.min(x)


def safe_values(model, x_grid):
    values = []
    for x in x_grid:
        try:
            y = model.value(float(x))
            values.append(y if np.isfinite(y) else np.nan)
        except Exception:
            values.append(np.nan)
    return np.array(values)


def regression_metrics(y_true, y_pred):
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    y_true = y_true[mask]
    y_pred = y_pred[mask]

    mse = mean_squared_error(y_true, y_pred)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(y_true, y_pred)

    return mse, rmse, mae


def monotonicity_violations_decreasing(y_pred):
    """
    For Distance vs Price, expected trend is generally decreasing:
    price should not increase as distance increases.
    Therefore violation = positive difference.
    """

    diffs = np.diff(y_pred)
    diffs = diffs[np.isfinite(diffs)]
    return int(np.sum(diffs > 1e-8))


def runtime_seconds(func):
    start = time.perf_counter()
    result = func()
    end = time.perf_counter()
    return result, end - start


def clip_extreme(y_values, y_min, y_max):
    return np.clip(y_values, y_min, y_max)


# ============================================================
# Main benchmark
# ============================================================

def main():
    df = load_housing_data()
    curve, x_raw, y = create_distance_price_curve(df, bins=35)

    x = normalise_x(x_raw)

    x_grid = np.linspace(x.min(), x.max(), 700)
    distance_grid = x_grid + x_raw.min()

    # Reference curve for metrics
    y_reference = np.interp(x_grid, x, y)

    # Plot clipping range for unstable methods
    y_min = max(0, np.nanmin(y) * 0.3)
    y_max = np.nanmax(y) * 2.5

    results = []

    # ============================================================
    # Method 1: Our cubic spline
    # ============================================================

    def run_our_cubic():
        model = CubicSpline(x, y, type=1)
        return safe_values(model, x_grid)

    y_our_cubic, rt = runtime_seconds(run_our_cubic)
    mse, rmse, mae = regression_metrics(y_reference, y_our_cubic)

    results.append({
        "dataset": "Melbourne housing Distance vs Price",
        "method": "Our cubic spline",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "monotonicity_violations_decreasing": monotonicity_violations_decreasing(y_our_cubic),
        "runtime_seconds": rt,
    })

    # ============================================================
    # Method 2: Our monotone spline after isotonic regression
    # ============================================================

    
    ir = IsotonicRegression(increasing=False)
    y_iso = ir.fit_transform(x, y)

    def run_iso_monotone():
        model = CubicSpline(x, y_iso, type=3)
        return safe_values(model, x_grid)

    y_iso_monotone, rt = runtime_seconds(run_iso_monotone)
    mse, rmse, mae = regression_metrics(y_reference, y_iso_monotone)

    results.append({
        "dataset": "Melbourne housing Distance vs Price",
        "method": "Isotonic + our monotone spline",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "monotonicity_violations_decreasing": monotonicity_violations_decreasing(y_iso_monotone),
        "runtime_seconds": rt,
    })

    # ============================================================
    # Method 3: Our taut spline experimental
    # ============================================================

    def run_our_taut():
        model = CubicSpline(x, y, type=2, tau=10.0)
        return safe_values(model, x_grid)

    y_our_taut, rt = runtime_seconds(run_our_taut)
    mse, rmse, mae = regression_metrics(y_reference, y_our_taut)

    results.append({
        "dataset": "Melbourne housing Distance vs Price",
        "method": "Our taut spline tau=10.0",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "monotonicity_violations_decreasing": monotonicity_violations_decreasing(y_our_taut),
        "runtime_seconds": rt,
    })

    # ============================================================
    # Method 4: SciPy CubicSpline
    # ============================================================

    def run_scipy_cubic():
        model = SciPyCubicSpline(x, y)
        return model(x_grid)

    y_scipy_cubic, rt = runtime_seconds(run_scipy_cubic)
    mse, rmse, mae = regression_metrics(y_reference, y_scipy_cubic)

    results.append({
        "dataset": "Melbourne housing Distance vs Price",
        "method": "SciPy CubicSpline",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "monotonicity_violations_decreasing": monotonicity_violations_decreasing(y_scipy_cubic),
        "runtime_seconds": rt,
    })

    # ============================================================
    # Method 5: SciPy PCHIP
    # ============================================================

    def run_scipy_pchip():
        model = PchipInterpolator(x, y)
        return model(x_grid)

    y_scipy_pchip, rt = runtime_seconds(run_scipy_pchip)
    mse, rmse, mae = regression_metrics(y_reference, y_scipy_pchip)

    results.append({
        "dataset": "Melbourne housing Distance vs Price",
        "method": "SciPy PCHIP",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "monotonicity_violations_decreasing": monotonicity_violations_decreasing(y_scipy_pchip),
        "runtime_seconds": rt,
    })

    # ============================================================
    # Method 6: Smoothing spline
    # ============================================================

    def run_smoothing_spline():
       
        model = UnivariateSpline(x, y, s=len(x) * np.var(y) * 0.2)
        return model(x_grid)

    y_smoothing, rt = runtime_seconds(run_smoothing_spline)
    mse, rmse, mae = regression_metrics(y_reference, y_smoothing)

    results.append({
        "dataset": "Melbourne housing Distance vs Price",
        "method": "SciPy smoothing spline",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "monotonicity_violations_decreasing": monotonicity_violations_decreasing(y_smoothing),
        "runtime_seconds": rt,
    })

    # ============================================================
    # Method 7: Polynomial regression baseline
    # ============================================================

    def run_polynomial_regression():
        model = make_pipeline(
            PolynomialFeatures(degree=3),
            LinearRegression()
        )

        model.fit(x.reshape(-1, 1), y)
        return model.predict(x_grid.reshape(-1, 1))

    y_poly, rt = runtime_seconds(run_polynomial_regression)
    mse, rmse, mae = regression_metrics(y_reference, y_poly)

    results.append({
        "dataset": "Melbourne housing Distance vs Price",
        "method": "Polynomial regression degree 3",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "monotonicity_violations_decreasing": monotonicity_violations_decreasing(y_poly),
        "runtime_seconds": rt,
    })

    # ============================================================
    # Save results
    # ============================================================

    results_df = pd.DataFrame(results)
    results_df.to_csv(RESULTS_CSV, index=False)

    print("\nBenchmark results:")
    print(results_df)
    print("\nSaved results to:", RESULTS_CSV)

    # ============================================================
    # Plot 1: raw housing scatter
    # ============================================================

    plt.figure(figsize=(10, 6))
    plt.scatter(df["Distance"], df["Price"], alpha=0.25, s=12, label="Individual properties")
    plt.scatter(x_raw, y, s=45, label="Median price by distance bin")
    plt.title("Melbourne Housing: Raw Price vs Distance from CBD")
    plt.xlabel("Distance from CBD")
    plt.ylabel("Price")
    plt.grid(True)
    plt.legend()

    path = PLOTS_DIR / "melbourne_housing_raw_scatter_and_binned_median.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 2: full method comparison
    # ============================================================

    plt.figure(figsize=(11, 6))
    plt.scatter(x_raw, y, label="Median price by distance bin", s=40)

    plt.plot(distance_grid, y_our_cubic, label="Our cubic spline")
    plt.plot(distance_grid, y_iso_monotone, label="Isotonic + our monotone spline")
    plt.plot(distance_grid, clip_extreme(y_our_taut, y_min, y_max), label="Our taut spline tau=10.0")
    plt.plot(distance_grid, y_scipy_cubic, "--", label="SciPy CubicSpline")
    plt.plot(distance_grid, y_scipy_pchip, "--", label="SciPy PCHIP")
    plt.plot(distance_grid, y_smoothing, label="SciPy smoothing spline")
    plt.plot(distance_grid, y_poly, label="Polynomial regression degree 3")

    plt.title("Melbourne Housing: Distance vs Price Approximation")
    plt.xlabel("Distance from CBD")
    plt.ylabel("Median Price")
    plt.grid(True)
    plt.legend()

    path = PLOTS_DIR / "melbourne_housing_distance_price_method_comparison.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 3: clean comparison without taut
    # ============================================================

    plt.figure(figsize=(11, 6))
    plt.scatter(x_raw, y, label="Median price by distance bin", s=40)

    plt.plot(distance_grid, y_our_cubic, label="Our cubic spline")
    plt.plot(distance_grid, y_iso_monotone, label="Isotonic + our monotone spline")
    plt.plot(distance_grid, y_scipy_cubic, "--", label="SciPy CubicSpline")
    plt.plot(distance_grid, y_scipy_pchip, "--", label="SciPy PCHIP")
    plt.plot(distance_grid, y_smoothing, label="SciPy smoothing spline")
    plt.plot(distance_grid, y_poly, label="Polynomial regression degree 3")

    plt.title("Melbourne Housing: Stable Method Comparison")
    plt.xlabel("Distance from CBD")
    plt.ylabel("Median Price")
    plt.grid(True)
    plt.legend()

    path = PLOTS_DIR / "melbourne_housing_stable_method_comparison.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 4: monotonicity violations
    # ============================================================

    plt.figure(figsize=(10, 6))
    plt.bar(
        results_df["method"],
        results_df["monotonicity_violations_decreasing"]
    )
    plt.xticks(rotation=35, ha="right")
    plt.title("Melbourne Housing: Decreasing Trend Violations by Method")
    plt.ylabel("Violations of expected decreasing trend")
    plt.grid(axis="y")

    path = PLOTS_DIR / "melbourne_housing_monotonicity_violations.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 5: MAE comparison
    # ============================================================

    plt.figure(figsize=(10, 6))
    plt.bar(
        results_df["method"],
        results_df["mae"]
    )
    plt.xticks(rotation=35, ha="right")
    plt.title("Melbourne Housing: MAE Comparison")
    plt.ylabel("Mean Absolute Error")
    plt.grid(axis="y")

    path = PLOTS_DIR / "melbourne_housing_mae_comparison.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 6: runtime comparison
    # ============================================================

    plt.figure(figsize=(10, 6))
    plt.bar(
        results_df["method"],
        results_df["runtime_seconds"]
    )
    plt.xticks(rotation=35, ha="right")
    plt.title("Melbourne Housing: Runtime Comparison")
    plt.ylabel("Runtime seconds")
    plt.grid(axis="y")

    path = PLOTS_DIR / "melbourne_housing_runtime_comparison.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    print("\nAll plots saved in:", PLOTS_DIR)


if __name__ == "__main__":
    main()