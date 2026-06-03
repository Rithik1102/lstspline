import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.interpolate import CubicSpline as SciPyCubicSpline
from scipy.interpolate import PchipInterpolator, UnivariateSpline
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lstspline import CubicSpline


DATA_PATH = PROJECT_ROOT / "data" / "API_FP.CPI.TOTL.ZG_DS2_en_csv_v2_115367.csv"

OUT_DIR = PROJECT_ROOT / "benchmarks" / "inflation_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PLOTS_DIR = OUT_DIR / "plots"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

RESULTS_CSV = OUT_DIR / "australia_inflation_benchmark_results.csv"


def load_world_bank_inflation(country_name="Australia"):
    df = pd.read_csv(DATA_PATH, skiprows=4)

    row = df[df["Country Name"] == country_name]

    if row.empty:
        raise ValueError(f"Country not found: {country_name}")

    year_cols = [col for col in df.columns if col.isdigit()]
    values = row[year_cols].iloc[0]

    data = pd.DataFrame({
        "year": year_cols,
        "inflation": values
    })

    data["year"] = data["year"].astype(int)
    data["inflation"] = pd.to_numeric(data["inflation"], errors="coerce")

    data = data.dropna()
    data = data.sort_values("year")

    return data


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


def turning_points(y_pred):
    """
    Counts approximate turning points.
    Useful for non-monotone economic data such as inflation.
    More turning points usually means more oscillatory behaviour.
    """

    dy = np.diff(y_pred)
    dy = dy[np.isfinite(dy)]

    signs = np.sign(dy)
    signs = signs[signs != 0]

    if len(signs) < 2:
        return 0

    return int(np.sum(signs[1:] != signs[:-1]))


def runtime_seconds(func):
    start = time.perf_counter()
    result = func()
    end = time.perf_counter()
    return result, end - start


def main():
    data = load_world_bank_inflation("Australia")

    years = data["year"].to_numpy(dtype=float)
    inflation = data["inflation"].to_numpy(dtype=float)

    x = years - years.min()
    y = inflation

    x_grid = np.linspace(x.min(), x.max(), 800)
    year_grid = x_grid + years.min()

    y_reference = np.interp(x_grid, x, y)

    results = []

    # ============================================================
    # Our cubic spline
    # ============================================================

    def run_our_cubic():
        model = CubicSpline(x, y, type=1)
        return safe_values(model, x_grid)

    y_our_cubic, rt = runtime_seconds(run_our_cubic)
    mse, rmse, mae = regression_metrics(y_reference, y_our_cubic)

    results.append({
        "dataset": "Australia inflation",
        "method": "Our cubic spline",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "turning_points": turning_points(y_our_cubic),
        "runtime_seconds": rt,
    })

    # ============================================================
    # Our taut spline experimental
    # ============================================================

    def run_our_taut():
        model = CubicSpline(x, y, type=2, tau=10.0)
        return safe_values(model, x_grid)

    y_our_taut, rt = runtime_seconds(run_our_taut)
    mse, rmse, mae = regression_metrics(y_reference, y_our_taut)

    results.append({
        "dataset": "Australia inflation",
        "method": "Our taut spline tau=10.0",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "turning_points": turning_points(y_our_taut),
        "runtime_seconds": rt,
    })

    # ============================================================
    # SciPy CubicSpline
    # ============================================================

    def run_scipy_cubic():
        model = SciPyCubicSpline(x, y)
        return model(x_grid)

    y_scipy_cubic, rt = runtime_seconds(run_scipy_cubic)
    mse, rmse, mae = regression_metrics(y_reference, y_scipy_cubic)

    results.append({
        "dataset": "Australia inflation",
        "method": "SciPy CubicSpline",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "turning_points": turning_points(y_scipy_cubic),
        "runtime_seconds": rt,
    })

    # ============================================================
    # SciPy PCHIP
    # ============================================================

    def run_scipy_pchip():
        model = PchipInterpolator(x, y)
        return model(x_grid)

    y_scipy_pchip, rt = runtime_seconds(run_scipy_pchip)
    mse, rmse, mae = regression_metrics(y_reference, y_scipy_pchip)

    results.append({
        "dataset": "Australia inflation",
        "method": "SciPy PCHIP",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "turning_points": turning_points(y_scipy_pchip),
        "runtime_seconds": rt,
    })

    # ============================================================
    # Smoothing spline
    # ============================================================

    def run_smoothing_spline():
        model = UnivariateSpline(x, y, s=len(x) * np.var(y) * 0.25)
        return model(x_grid)

    y_smoothing, rt = runtime_seconds(run_smoothing_spline)
    mse, rmse, mae = regression_metrics(y_reference, y_smoothing)

    results.append({
        "dataset": "Australia inflation",
        "method": "SciPy smoothing spline",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "turning_points": turning_points(y_smoothing),
        "runtime_seconds": rt,
    })

    # ============================================================
    # Polynomial regression baseline
    # ============================================================

    def run_polynomial_regression():
        model = make_pipeline(
            PolynomialFeatures(degree=4),
            LinearRegression()
        )
        model.fit(x.reshape(-1, 1), y)
        return model.predict(x_grid.reshape(-1, 1))

    y_poly, rt = runtime_seconds(run_polynomial_regression)
    mse, rmse, mae = regression_metrics(y_reference, y_poly)

    results.append({
        "dataset": "Australia inflation",
        "method": "Polynomial regression degree 4",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "turning_points": turning_points(y_poly),
        "runtime_seconds": rt,
    })

    # ============================================================
    # Moving-average baseline
    # ============================================================

    smooth_series = (
        pd.Series(y)
        .rolling(window=5, center=True, min_periods=1)
        .mean()
        .to_numpy()
    )

    y_moving_average = np.interp(x_grid, x, smooth_series)
    mse, rmse, mae = regression_metrics(y_reference, y_moving_average)

    results.append({
        "dataset": "Australia inflation",
        "method": "Moving average baseline",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "turning_points": turning_points(y_moving_average),
        "runtime_seconds": 0.0,
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
    # Plot 1: raw inflation data
    # ============================================================

    plt.figure(figsize=(10, 6))
    plt.scatter(years, y, label="Observed inflation rate", s=35)
    plt.plot(years, y, alpha=0.6)
    plt.title("Australia Inflation Rate Data")
    plt.xlabel("Year")
    plt.ylabel("Inflation rate (% annual)")
    plt.grid(True)
    plt.legend()

    path = PLOTS_DIR / "australia_inflation_raw_data.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 2: full method comparison
    # ============================================================

    y_min = np.nanmin(y) - 5
    y_max = np.nanmax(y) + 5
    y_our_taut_clipped = np.clip(y_our_taut, y_min, y_max)

    plt.figure(figsize=(11, 6))
    plt.scatter(years, y, label="Observed inflation", s=25)

    plt.plot(year_grid, y_our_cubic, label="Our cubic spline")
    plt.plot(year_grid, y_our_taut_clipped, label="Our taut spline tau=10.0")
    plt.plot(year_grid, y_scipy_cubic, "--", label="SciPy CubicSpline")
    plt.plot(year_grid, y_scipy_pchip, "--", label="SciPy PCHIP")
    plt.plot(year_grid, y_smoothing, label="SciPy smoothing spline")
    plt.plot(year_grid, y_poly, label="Polynomial regression degree 4")
    plt.plot(year_grid, y_moving_average, label="Moving average baseline")

    plt.title("Australia Inflation: Interpolation and Smoothing Comparison")
    plt.xlabel("Year")
    plt.ylabel("Inflation rate (% annual)")
    plt.grid(True)
    plt.legend()

    path = PLOTS_DIR / "australia_inflation_method_comparison.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 3: stable smoothing/approximation methods
    # ============================================================

    plt.figure(figsize=(11, 6))
    plt.scatter(years, y, label="Observed inflation", s=25)

    plt.plot(year_grid, y_our_cubic, label="Our cubic spline")
    plt.plot(year_grid, y_scipy_cubic, "--", label="SciPy CubicSpline")
    plt.plot(year_grid, y_scipy_pchip, "--", label="SciPy PCHIP")
    plt.plot(year_grid, y_smoothing, label="SciPy smoothing spline")
    plt.plot(year_grid, y_poly, label="Polynomial regression degree 4")
    plt.plot(year_grid, y_moving_average, label="Moving average baseline")

    plt.title("Australia Inflation: Stable Method Comparison")
    plt.xlabel("Year")
    plt.ylabel("Inflation rate (% annual)")
    plt.grid(True)
    plt.legend()

    path = PLOTS_DIR / "australia_inflation_stable_method_comparison.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 4: zoomed recent years
    # ============================================================

    mask = year_grid >= 1990

    plt.figure(figsize=(11, 6))
    plt.scatter(years[years >= 1990], y[years >= 1990],label="Observed GDP", s=25)

    plt.plot(year_grid[mask], y_our_cubic[mask], label="Our cubic spline")
    plt.plot(year_grid[mask], y_scipy_cubic[mask], "--", label="SciPy CubicSpline")
    plt.plot(year_grid[mask], y_scipy_pchip[mask], "--", label="SciPy PCHIP")
    plt.plot(year_grid[mask], y_smoothing[mask], label="SciPy smoothing spline")
    plt.plot(year_grid[mask], y_poly[mask], label="Polynomial regression degree 4")
    plt.plot(year_grid[mask], y_moving_average[mask], label="Moving average baseline")

    plt.title("Australia Inflation: Zoomed Comparison from 1990 Onwards")
    plt.xlabel("Year")
    plt.ylabel("Inflation rate (% annual)")
    plt.grid(True)
    plt.legend()

    path = PLOTS_DIR / "australia_inflation_zoom_1990_onwards.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 5: MAE comparison
    # ============================================================

    plt.figure(figsize=(10, 6))
    plt.bar(results_df["method"], results_df["mae"])
    plt.xticks(rotation=35, ha="right")
    plt.title("Australia Inflation: MAE Comparison")
    plt.ylabel("Mean Absolute Error")
    plt.grid(axis="y")

    path = PLOTS_DIR / "australia_inflation_mae_comparison.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 6: Turning point comparison
    # ============================================================

    plt.figure(figsize=(10, 6))
    plt.bar(results_df["method"], results_df["turning_points"])
    plt.xticks(rotation=35, ha="right")
    plt.title("Australia Inflation: Turning Points by Method")
    plt.ylabel("Number of turning points")
    plt.grid(axis="y")

    path = PLOTS_DIR / "australia_inflation_turning_points.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 7: Runtime comparison
    # ============================================================

    plt.figure(figsize=(10, 6))
    plt.bar(results_df["method"], results_df["runtime_seconds"])
    plt.xticks(rotation=35, ha="right")
    plt.title("Australia Inflation: Runtime Comparison")
    plt.ylabel("Runtime seconds")
    plt.grid(axis="y")

    path = PLOTS_DIR / "australia_inflation_runtime_comparison.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    print("\nAll plots saved in:", PLOTS_DIR)


if __name__ == "__main__":
    main()