import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lstspline import CubicSpline


OUT_DIR = PROJECT_ROOT / "benchmarks" / "taut_tau_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TAU_VALUES = [0.05, 2.0, 10.0, 20.0]


def safe_values(spline, x_grid):
    vals = []
    for xx in x_grid:
        try:
            yy = spline.value(float(xx))
            vals.append(yy if np.isfinite(yy) else np.nan)
        except Exception:
            vals.append(np.nan)
    return np.array(vals)


def plot_taut_tau_experiment(name, x, y, tau_values=TAU_VALUES, zoom=None):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)

    x_grid = np.linspace(x.min(), x.max(), 800)

    plt.figure(figsize=(11, 6))
    plt.scatter(x, y, label="Observed data", s=45, zorder=5)

    for tau in tau_values:
        spline = CubicSpline(x, y, type=2, tau=tau)
        y_grid = safe_values(spline, x_grid)

        plt.plot(
            x_grid,
            y_grid,
            linewidth=2.2,
            label=f"Taut spline tau={tau}"
        )

    plt.title(f"Taut Spline Tau Experiment: {name}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()

    if zoom is not None:
        plt.xlim(zoom[0], zoom[1])

    fname = name.lower().replace(" ", "_").replace("/", "_").replace("(", "").replace(")", "").replace(",", "")
    if zoom is None:
        path = OUT_DIR / f"{fname}_tau_experiment.png"
    else:
        path = OUT_DIR / f"{fname}_tau_experiment_zoom.png"

    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()

    print("Saved:", path)


def load_world_bank_csv(data_path, country_name, value_name):
    df = pd.read_csv(data_path, skiprows=4)
    row = df[df["Country Name"] == country_name]

    if row.empty:
        raise ValueError(f"Country not found: {country_name}")

    year_cols = [col for col in df.columns if col.isdigit()]
    values = row[year_cols].iloc[0]

    data = pd.DataFrame({
        "year": year_cols,
        value_name: values
    })

    data["year"] = data["year"].astype(int)
    data[value_name] = pd.to_numeric(data[value_name], errors="coerce")

    data = data.dropna()
    data = data.sort_values("year")

    return data


def run_synthetic_tau_tests():
    # Dataset 1: plateau and sharp increase
    x1 = np.arange(0, 10)
    y1 = np.array([10, 10, 10, 10, 11, 25, 27, 28, 40, 50], dtype=float)

    plot_taut_tau_experiment(
        "Plateau then sharp increase",
        x1,
        y1
    )

    plot_taut_tau_experiment(
        "Plateau then sharp increase",
        x1,
        y1,
        zoom=(3, 7)
    )

    # Dataset 2: irregular monotone increase
    x2 = np.array([0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15], dtype=float)
    y2 = np.array([10, 10, 10, 10, 10.5, 15, 50, 60, 85, 90, 95], dtype=float)

    plot_taut_tau_experiment(
        "Irregular monotone increase",
        x2,
        y2
    )

    plot_taut_tau_experiment(
        "Irregular monotone increase",
        x2,
        y2,
        zoom=(6, 10)
    )

    # Dataset 3: decreasing monotone data
    y2_decreasing = -y2

    plot_taut_tau_experiment(
        "Decreasing monotone data",
        x2,
        y2_decreasing
    )

    plot_taut_tau_experiment(
        "Decreasing monotone data",
        x2,
        y2_decreasing,
        zoom=(6, 10)
    )

    # Dataset 4: non-monotone max(cos(x), 0)
    x_cos = np.linspace(-3, 3, 25)
    y_cos = np.maximum(np.cos(x_cos), 0)

    plot_taut_tau_experiment(
        "max cos x 0",
        x_cos,
        y_cos
    )

    plot_taut_tau_experiment(
        "max cos x 0",
        x_cos,
        y_cos,
        zoom=(-1.5, 1.5)
    )


def run_population_tau_test():
    data_path = PROJECT_ROOT / "data" / "API_SP.POP.TOTL_DS2_en_csv_v2_127039.csv"

    if not data_path.exists():
        print("Population file not found, skipping:", data_path)
        return

    data = load_world_bank_csv(
        data_path=data_path,
        country_name="Australia",
        value_name="population"
    )

    years = data["year"].to_numpy(dtype=float)
    population = data["population"].to_numpy(dtype=float) / 1e6

    x = years - years.min()
    y = population

    plot_taut_tau_experiment(
        "Australia population",
        x,
        y
    )

    zoom_start = 2000 - years.min()
    zoom_end = years.max() - years.min()

    plot_taut_tau_experiment(
        "Australia population",
        x,
        y,
        zoom=(zoom_start, zoom_end)
    )


def run_gdp_tau_test():
    data_path = PROJECT_ROOT / "data" / "API_NY.GDP.MKTP.CD_DS2_en_csv_v2_126992.csv"

    if not data_path.exists():
        print("GDP file not found, skipping:", data_path)
        return

    data = load_world_bank_csv(
        data_path=data_path,
        country_name="Australia",
        value_name="gdp"
    )

    years = data["year"].to_numpy(dtype=float)
    gdp = data["gdp"].to_numpy(dtype=float) / 1e12

    x = years - years.min()
    y = gdp

    plot_taut_tau_experiment(
        "Australia GDP",
        x,
        y
    )

    zoom_start = 2000 - years.min()
    zoom_end = years.max() - years.min()

    plot_taut_tau_experiment(
        "Australia GDP",
        x,
        y,
        zoom=(zoom_start, zoom_end)
    )


def main():
    print("Running taut spline tau experiments...")

    run_synthetic_tau_tests()
    run_population_tau_test()
    run_gdp_tau_test()

    print("\nDone. Plots saved in:", OUT_DIR)


if __name__ == "__main__":
    main()