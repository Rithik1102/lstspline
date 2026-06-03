import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline as SciPyCubicSpline
from scipy.interpolate import PchipInterpolator
from sklearn.isotonic import IsotonicRegression

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lstspline import CubicSpline


# ============================================================
# Paths
# ============================================================

DATA_PATH = PROJECT_ROOT / "data" / "API_SP.POP.TOTL_DS2_en_csv_v2_127039.csv"

OUT_DIR = PROJECT_ROOT / "benchmarks" / "population_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PLOTS_DIR = OUT_DIR / "plots"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

RESULTS_CSV = OUT_DIR / "australia_population_benchmark_results.csv"


# ============================================================
# Helper functions
# ============================================================

def load_world_bank_population(country_name="Australia"):
    """
    Loads World Bank population dataset and extracts one country.
    World Bank CSV has metadata rows before the actual header.
    """

    df = pd.read_csv(DATA_PATH, skiprows=4)

    row = df[df["Country Name"] == country_name]

    if row.empty:
        raise ValueError(f"Country not found: {country_name}")

    year_cols = [col for col in df.columns if col.isdigit()]

    values = row[year_cols].iloc[0]

    data = pd.DataFrame({
        "year": year_cols,
        "population": values
    })

    data["year"] = data["year"].astype(int)
    data["population"] = pd.to_numeric(data["population"], errors="coerce")

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


def mse(y_true, y_pred):
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    return np.mean((y_true[mask] - y_pred[mask]) ** 2)


def mae(y_true, y_pred):
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    return np.mean(np.abs(y_true[mask] - y_pred[mask]))


def monotonicity_violations(y_pred):
    diffs = np.diff(y_pred)
    diffs = diffs[np.isfinite(diffs)]
    return int(np.sum(diffs < -1e-8))


def runtime_seconds(func):
    start = time.perf_counter()
    result = func()
    end = time.perf_counter()
    return result, end - start


# ============================================================
# Main benchmark
# ============================================================

def main():
    country = "Australia"

    data = load_world_bank_population(country)

    x_years = data["year"].to_numpy(dtype=float)
    y_population = data["population"].to_numpy(dtype=float)

    
    x = x_years - x_years.min()
    y = y_population

    x_grid = np.linspace(x.min(), x.max(), 800)
    year_grid = x_grid + x_years.min()


    y_reference = np.interp(x_grid, x, y)

    results = []

    # ------------------------------------------------------------
    # Our Cubic Spline
    # ------------------------------------------------------------
    def run_our_cubic():
        model = CubicSpline(x, y, type=1)
        return safe_values(model, x_grid)

    y_our_cubic, rt = runtime_seconds(run_our_cubic)

    results.append({
        "dataset": "Australia population",
        "method": "Our cubic spline",
        "mse": mse(y_reference, y_our_cubic),
        "mae": mae(y_reference, y_our_cubic),
        "monotonicity_violations": monotonicity_violations(y_our_cubic),
        "runtime_seconds": rt,
    })

    # ------------------------------------------------------------
    # Our Monotone Spline
    # ------------------------------------------------------------
    def run_our_monotone():
        model = CubicSpline(x, y, type=3)
        return safe_values(model, x_grid)

    y_our_monotone, rt = runtime_seconds(run_our_monotone)

    results.append({
        "dataset": "Australia population",
        "method": "Our monotone spline",
        "mse": mse(y_reference, y_our_monotone),
        "mae": mae(y_reference, y_our_monotone),
        "monotonicity_violations": monotonicity_violations(y_our_monotone),
        "runtime_seconds": rt,
    })

    # ------------------------------------------------------------
    # Our Taut Spline - included as experimental
    # ------------------------------------------------------------
    def run_our_taut():
        model = CubicSpline(x, y, type=2, tau=10.0)
        return safe_values(model, x_grid)

    y_our_taut, rt = runtime_seconds(run_our_taut)

    results.append({
        "dataset": "Australia population",
        "method": "Our taut spline tau=10.0",
        "mse": mse(y_reference, y_our_taut),
        "mae": mae(y_reference, y_our_taut),
        "monotonicity_violations": monotonicity_violations(y_our_taut),
        "runtime_seconds": rt,
    })

    # ------------------------------------------------------------
    # SciPy CubicSpline
    # ------------------------------------------------------------
    def run_scipy_cubic():
        model = SciPyCubicSpline(x, y)
        return model(x_grid)

    y_scipy_cubic, rt = runtime_seconds(run_scipy_cubic)

    results.append({
        "dataset": "Australia population",
        "method": "SciPy CubicSpline",
        "mse": mse(y_reference, y_scipy_cubic),
        "mae": mae(y_reference, y_scipy_cubic),
        "monotonicity_violations": monotonicity_violations(y_scipy_cubic),
        "runtime_seconds": rt,
    })

    # ------------------------------------------------------------
    # SciPy PCHIP
    # ------------------------------------------------------------
    def run_scipy_pchip():
        model = PchipInterpolator(x, y)
        return model(x_grid)

    y_scipy_pchip, rt = runtime_seconds(run_scipy_pchip)

    results.append({
        "dataset": "Australia population",
        "method": "SciPy PCHIP",
        "mse": mse(y_reference, y_scipy_pchip),
        "mae": mae(y_reference, y_scipy_pchip),
        "monotonicity_violations": monotonicity_violations(y_scipy_pchip),
        "runtime_seconds": rt,
    })

    # ------------------------------------------------------------
    # Isotonic Regression + Our Monotone Spline
    # ------------------------------------------------------------
    ir = IsotonicRegression(increasing=True)
    y_iso = ir.fit_transform(x, y)

    def run_iso_monotone():
        model = CubicSpline(x, y_iso, type=3)
        return safe_values(model, x_grid)

    y_iso_monotone, rt = runtime_seconds(run_iso_monotone)

    results.append({
        "dataset": "Australia population",
        "method": "Isotonic + our monotone spline",
        "mse": mse(y_reference, y_iso_monotone),
        "mae": mae(y_reference, y_iso_monotone),
        "monotonicity_violations": monotonicity_violations(y_iso_monotone),
        "runtime_seconds": rt,
    })

    # ============================================================
    # Save benchmark results
    # ============================================================

    results_df = pd.DataFrame(results)
    results_df.to_csv(RESULTS_CSV, index=False)

    print("\nBenchmark results:")
    print(results_df)
    print("\nSaved results to:", RESULTS_CSV)

    # ============================================================
    # Plot 1: raw data
    # ============================================================

    plt.figure(figsize=(10, 6))
    plt.scatter(x_years, y_population, label="Observed population data", s=35)
    plt.title("Australia Population Data from World Bank")
    plt.xlabel("Year")
    plt.ylabel("Population")
    plt.grid(True)
    plt.legend()

    path = PLOTS_DIR / "australia_population_raw_data.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 2: full method comparison
    # ============================================================

    plt.figure(figsize=(11, 6))
    plt.scatter(x_years, y_population, label="Observed data", s=25)
    plt.plot(year_grid, y_our_cubic, label="Our cubic spline")
    plt.plot(year_grid, y_our_monotone, label="Our monotone spline")
    plt.plot(year_grid, y_our_taut, label="Our taut spline tau=10.0")
    plt.plot(year_grid, y_scipy_cubic, "--", label="SciPy CubicSpline")
    plt.plot(year_grid, y_scipy_pchip, "--", label="SciPy PCHIP")
    plt.plot(year_grid, y_iso_monotone, ":", label="Isotonic + our monotone spline")

    plt.title("Australia Population: Spline Method Comparison")
    plt.xlabel("Year")
    plt.ylabel("Population")
    plt.grid(True)
    plt.legend()

    path = PLOTS_DIR / "australia_population_spline_comparison.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 3: zoomed recent years
    # ============================================================

    zoom_min_year = 2000
    zoom_max_year = int(x_years.max())

    mask = (year_grid >= zoom_min_year) & (year_grid <= zoom_max_year)

    plt.figure(figsize=(11, 6))
    plt.scatter(x_years[x_years >= 2000], y_population[x_years >= 2000], label="Observed population data", s=25)
    plt.plot(year_grid[mask], y_our_cubic[mask], label="Our cubic spline")
    plt.plot(year_grid[mask], y_our_monotone[mask], label="Our monotone spline")
    plt.plot(year_grid[mask], y_our_taut[mask], label="Our taut spline tau=10.0")
    plt.plot(year_grid[mask], y_scipy_cubic[mask], "--", label="SciPy CubicSpline")
    plt.plot(year_grid[mask], y_scipy_pchip[mask], "--", label="SciPy PCHIP")
    plt.plot(year_grid[mask], y_iso_monotone[mask], ":", label="Isotonic + monotone")

    plt.title("Australia Population: Zoomed Comparison from 2000 Onwards")
    plt.xlabel("Year")
    plt.ylabel("Population")
    plt.grid(True)
    plt.legend()

    path = PLOTS_DIR / "australia_population_zoom_2000_onwards.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    # ============================================================
    # Plot 4: monotonicity comparison
    # ============================================================

    plt.figure(figsize=(10, 6))
    plt.bar(
        results_df["method"],
        results_df["monotonicity_violations"]
    )
    plt.xticks(rotation=35, ha="right")
    plt.title("Australia Population: Monotonicity Violations by Method")
    plt.ylabel("Number of monotonicity violations")
    plt.grid(axis="y")

    path = PLOTS_DIR / "australia_population_monotonicity_violations.png"
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
    plt.title("Australia Population: MAE Comparison by Method")
    plt.ylabel("Mean Absolute Error")
    plt.grid(axis="y")

    path = PLOTS_DIR / "australia_population_mae_comparison.png"
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)

    print("\nAll plots saved in:", PLOTS_DIR)


if __name__ == "__main__":
    main()