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
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lstspline import CubicSpline


DATA_PATH = PROJECT_ROOT / "data" / "API_NY.GDP.MKTP.CD_DS2_en_csv_v2_126992.csv"

OUT_DIR = PROJECT_ROOT / "benchmarks" / "gdp_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PLOTS_DIR = OUT_DIR / "plots"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

RESULTS_CSV = OUT_DIR / "australia_gdp_benchmark_results.csv"


def load_world_bank_gdp(country_name="Australia"):
    df = pd.read_csv(DATA_PATH, skiprows=4)

    row = df[df["Country Name"] == country_name]

    if row.empty:
        raise ValueError(f"Country not found: {country_name}")

    year_cols = [col for col in df.columns if col.isdigit()]
    values = row[year_cols].iloc[0]

    data = pd.DataFrame({
        "year": year_cols,
        "gdp": values
    })

    data["year"] = data["year"].astype(int)
    data["gdp"] = pd.to_numeric(data["gdp"], errors="coerce")

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


def monotonicity_violations_increasing(y_pred):
    diffs = np.diff(y_pred)
    diffs = diffs[np.isfinite(diffs)]
    return int(np.sum(diffs < -1e-8))


def runtime_seconds(func):
    start = time.perf_counter()
    result = func()
    end = time.perf_counter()
    return result, end - start


def main():
    data = load_world_bank_gdp("Australia")

    years = data["year"].to_numpy(dtype=float)
    gdp = data["gdp"].to_numpy(dtype=float)

    
    gdp_trillion = gdp / 1e12

    # Normalise x for numerical stability
    x = years - years.min()
    y = gdp_trillion

    x_grid = np.linspace(x.min(), x.max(), 800)
    year_grid = x_grid + years.min()

    # Reference curve for metric comparison
    y_reference = np.interp(x_grid, x, y)

    results = []

    # ============================================================
    # Our Cubic Spline
    # ============================================================

    def run_our_cubic():
        model = CubicSpline(x, y, type=1)
        return safe_values(model, x_grid)

    y_our_cubic, rt = runtime_seconds(run_our_cubic)
    mse, rmse, mae = regression_metrics(y_reference, y_our_cubic)

    results.append({
        "dataset": "Australia GDP",
        "method": "Our cubic spline",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "monotonicity_violations": monotonicity_violations_increasing(y_our_cubic),
        "runtime_seconds": rt,
    })

    # ============================================================
    # Our Monotone Spline
    # ============================================================

    def run_our_monotone():
        model = CubicSpline(x, y, type=3)
        return safe_values(model, x_grid)

    y_our_monotone, rt = runtime_seconds(run_our_monotone)
    mse, rmse, mae = regression_metrics(y_reference, y_our_monotone)

    results.append({
        "dataset": "Australia GDP",
        "method": "Our monotone spline",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "monotonicity_violations": monotonicity_violations_increasing(y_our_monotone),
        "runtime_seconds": rt,
    })

    # ============================================================
    # Our Taut Spline - experimental
    # ============================================================

    def run_our_taut():
        model = CubicSpline(x, y, type=2, tau=10.0)
        return safe_values(model, x_grid)

    y_our_taut, rt = runtime_seconds(run_our_taut)
    mse, rmse, mae = regression_metrics(y_reference, y_our_taut)

    results.append({
        "dataset": "Australia GDP",
        "method": "Our taut spline tau=10.0",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "monotonicity_violations": monotonicity_violations_increasing(y_our_taut),
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
        "dataset": "Australia GDP",
        "method": "SciPy CubicSpline",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "monotonicity_violations": monotonicity_violations_increasing(y_scipy_cubic),
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
        "dataset": "Australia GDP",
        "method": "SciPy PCHIP",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "monotonicity_violations": monotonicity_violations_increasing(y_scipy_pchip),
        "runtime_seconds": rt,
    })

    # ============================================================
    # Isotonic Regression + Our Monotone Spline
    # ============================================================

    ir = IsotonicRegression(increasing=True)
    y_iso = ir.fit_transform(x, y)

    def run_iso_monotone():
        model = CubicSpline(x, y_iso, type=3)
        return safe_values(model, x_grid)

    y_iso_monotone, rt = runtime_seconds(run_iso_monotone)
    mse, rmse, mae = regression_metrics(y_reference, y_iso_monotone)

    results.append({
        "dataset": "Australia GDP",
        "method": "Isotonic + our monotone spline",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "monotonicity_violations": monotonicity_violations_increasing(y_iso_monotone),
        "runtime_seconds": rt,
    })

    # ============================================================
    # SciPy Smoothing Spline
    # ============================================================

    def run_smoothing_spline():
        model = UnivariateSpline(x, y, s=len(x) * np.var(y) * 0.15)
        return model(x_grid)

    y_smoothing, rt = runtime_seconds(run_smoothing_spline)
    mse, rmse, mae = regression_metrics(y_reference, y_smoothing)

    results.append({
        "dataset": "Australia GDP",
        "method": "SciPy smoothing spline",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "monotonicity_violations": monotonicity_violations_increasing(y_smoothing),
        "runtime_seconds": rt,
    })

    # ============================================================
    # Polynomial Regression Baseline
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
        "dataset": "Australia GDP",
        "method": "Polynomial regression degree 3",
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "monotonicity_violations": monotonicity_violations_increasing(y_poly),
        "runtime_seconds": rt,
    })

    results_df = pd.DataFrame(results)
    results_df.to_csv(RESULTS_CSV, index=False)

    print("\nBenchmark results:")
    print(results_df)
    print("\nSaved results to:", RESULTS_CSV)

    # ============================================================
    # Plot 1: Raw GDP data
    # ============================================================

    plt.figure(figsize=(10, 6))
    plt.scatter(years, y, label="Observed GDP", s=35)
    plt.title("Australia GDP Data from World Bank")
    plt.xlabel("Year")
    plt.ylabel("GDP (trillion current US$)")
    plt.grid(True)
    plt.legend()

    path = PLOTS_DIR / "australia_gdp_raw_data.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 2: Full method comparison
    # ============================================================

    plt.figure(figsize=(11, 6))
    plt.scatter(years, y, label="Observed GDP", s=25)

    plt.plot(year_grid, y_our_cubic, label="Our cubic spline")
    plt.plot(year_grid, y_our_monotone, label="Our monotone spline")
    plt.plot(year_grid, y_our_taut, label="Our taut spline tau=10")
    plt.plot(year_grid, y_scipy_cubic, "--", label="SciPy CubicSpline")
    plt.plot(year_grid, y_scipy_pchip, "--", label="SciPy PCHIP")
    plt.plot(year_grid, y_iso_monotone, ":", label="Isotonic + our monotone spline")
    plt.plot(year_grid, y_smoothing, label="SciPy smoothing spline")
    plt.plot(year_grid, y_poly, label="Polynomial regression degree 3")

    plt.title("Australia GDP: Spline and Regression Method Comparison")
    plt.xlabel("Year")
    plt.ylabel("GDP (trillion current US$)")
    plt.grid(True)
    plt.legend()

    path = PLOTS_DIR / "australia_gdp_method_comparison.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 3: Stable methods only
    # ============================================================

    plt.figure(figsize=(11, 6))
    plt.scatter(years, y, label="Observed GDP", s=25)

    plt.plot(year_grid, y_our_cubic, label="Our cubic spline")
    plt.plot(year_grid, y_our_monotone, label="Our monotone spline")
    plt.plot(year_grid, y_scipy_cubic, "--", label="SciPy CubicSpline")
    plt.plot(year_grid, y_scipy_pchip, "--", label="SciPy PCHIP")
    plt.plot(year_grid, y_iso_monotone, ":", label="Isotonic + our monotone spline")
    plt.plot(year_grid, y_smoothing, label="SciPy smoothing spline")
    plt.plot(year_grid, y_poly, label="Polynomial regression degree 3")

    plt.title("Australia GDP: Stable Method Comparison")
    plt.xlabel("Year")
    plt.ylabel("GDP (trillion current US$)")
    plt.grid(True)
    plt.legend()

    path = PLOTS_DIR / "australia_gdp_stable_method_comparison.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 4: Zoom 2000 onwards
    # ============================================================

    mask = year_grid >= 2000

    plt.figure(figsize=(11, 6))
    plt.xlim(2000, 2025)
    plt.scatter(years[years >= 2000], y[years >= 2000],label="Observed GDP", s=25)

    plt.plot(year_grid[mask], y_our_cubic[mask], label="Our cubic spline")
    plt.plot(year_grid[mask], y_our_monotone[mask], label="Our monotone spline")
    plt.plot(year_grid[mask], y_scipy_cubic[mask], "--", label="SciPy CubicSpline")
    plt.plot(year_grid[mask], y_scipy_pchip[mask], "--", label="SciPy PCHIP")
    plt.plot(year_grid[mask], y_iso_monotone[mask], ":", label="Isotonic + our monotone spline")
    plt.plot(year_grid[mask], y_smoothing[mask], label="SciPy smoothing spline")
    plt.plot(year_grid[mask], y_poly[mask], label="Polynomial regression degree 3")

    plt.title("Australia GDP: Zoomed Comparison from 2000 Onwards")
    plt.xlabel("Year")
    plt.ylabel("GDP (trillion current US$)")
    plt.grid(True)
    plt.legend()

    path = PLOTS_DIR / "australia_gdp_zoom_2000_onwards.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 5: Monotonicity violations
    # ============================================================

    plt.figure(figsize=(10, 6))
    plt.bar(results_df["method"], results_df["monotonicity_violations"])
    plt.xticks(rotation=35, ha="right")
    plt.title("Australia GDP: Monotonicity Violations by Method")
    plt.ylabel("Number of monotonicity violations")
    plt.grid(axis="y")

    path = PLOTS_DIR / "australia_gdp_monotonicity_violations.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 6: MAE comparison
    # ============================================================

    plt.figure(figsize=(10, 6))
    plt.bar(results_df["method"], results_df["mae"])
    plt.xticks(rotation=35, ha="right")
    plt.title("Australia GDP: MAE Comparison")
    plt.ylabel("Mean Absolute Error")
    plt.grid(axis="y")

    path = PLOTS_DIR / "australia_gdp_mae_comparison.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 7: Runtime comparison
    # ============================================================

    plt.figure(figsize=(10, 6))
    plt.bar(results_df["method"], results_df["runtime_seconds"])
    plt.xticks(rotation=35, ha="right")
    plt.title("Australia GDP: Runtime Comparison")
    plt.ylabel("Runtime seconds")
    plt.grid(axis="y")

    path = PLOTS_DIR / "australia_gdp_runtime_comparison.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    print("\nAll plots saved in:", PLOTS_DIR)


if __name__ == "__main__":
    main()