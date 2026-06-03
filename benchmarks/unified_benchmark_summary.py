from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[1]

OUT_DIR = PROJECT_ROOT / "benchmarks" / "unified_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PLOTS_DIR = OUT_DIR / "plots"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_CSV = OUT_DIR / "unified_benchmark_summary.csv"
BEST_METHOD_CSV = OUT_DIR / "best_method_summary.csv"


# ============================================================
# Existing benchmark result CSV files
# ============================================================

CSV_FILES = {
    "Population": PROJECT_ROOT / "benchmarks" / "population_outputs" / "australia_population_benchmark_results.csv",
    "GDP": PROJECT_ROOT / "benchmarks" / "gdp_outputs" / "australia_gdp_benchmark_results.csv",
    "Inflation": PROJECT_ROOT / "benchmarks" / "inflation_outputs" / "australia_inflation_benchmark_results.csv",
    "Housing": PROJECT_ROOT / "benchmarks" / "housing_outputs" / "melbourne_housing_benchmark_results.csv",
}


# ============================================================
# Helper functions
# ============================================================

def clean_method_name(method):
    method = str(method)

    replacements = {
        "Our taut spline tau=10.0": "Our taut spline",
        "Polynomial regression degree 3": "Polynomial regression",
        "Polynomial regression degree 4": "Polynomial regression",
        "Isotonic + our monotone spline": "Isotonic + monotone spline",
        "SciPy CubicSpline": "SciPy cubic spline",
        "SciPy PCHIP": "SciPy PCHIP",
        "SciPy smoothing spline": "Smoothing spline",
        "Moving average baseline": "Moving average",
    }

    return replacements.get(method, method)


def load_all_results():
    all_frames = []

    for dataset_name, csv_path in CSV_FILES.items():
        if not csv_path.exists():
            print(f"WARNING: Missing file: {csv_path}")
            continue

        df = pd.read_csv(csv_path)

        df["dataset_group"] = dataset_name
        df["method_clean"] = df["method"].apply(clean_method_name)

        if "rmse" not in df.columns and "mse" in df.columns:
            df["rmse"] = np.sqrt(df["mse"])

        if "mae" not in df.columns:
            raise ValueError(f"Missing MAE column in {csv_path}")

        if "runtime_seconds" not in df.columns:
            df["runtime_seconds"] = np.nan


        max_mae = df["mae"].replace([np.inf, -np.inf], np.nan).max()
        max_rmse = df["rmse"].replace([np.inf, -np.inf], np.nan).max()

        df["relative_mae"] = df["mae"] / max_mae if max_mae and max_mae > 0 else np.nan
        df["relative_rmse"] = df["rmse"] / max_rmse if max_rmse and max_rmse > 0 else np.nan

      
        if "monotonicity_violations" in df.columns:
            df["stability_metric"] = df["monotonicity_violations"]
            df["stability_metric_name"] = "monotonicity violations"

        elif "monotonicity_violations_decreasing" in df.columns:
            df["stability_metric"] = df["monotonicity_violations_decreasing"]
            df["stability_metric_name"] = "decreasing trend violations"

        elif "turning_points" in df.columns:
            df["stability_metric"] = df["turning_points"]
            df["stability_metric_name"] = "turning points"

        else:
            df["stability_metric"] = np.nan
            df["stability_metric_name"] = "not available"

        keep_cols = [
            "dataset_group",
            "dataset",
            "method_clean",
            "method",
            "mse",
            "rmse",
            "mae",
            "relative_mae",
            "relative_rmse",
            "runtime_seconds",
            "stability_metric",
            "stability_metric_name",
        ]

        keep_cols = [c for c in keep_cols if c in df.columns]
        all_frames.append(df[keep_cols])

    if not all_frames:
        raise ValueError("No benchmark CSV files were found.")

    return pd.concat(all_frames, ignore_index=True)


def save_summary_table(df):
    df.to_csv(SUMMARY_CSV, index=False)
    print("Saved unified benchmark summary:", SUMMARY_CSV)

    print("\nUnified benchmark preview:")
    print(df.head(30))


def plot_grouped_bar(df, metric, title, ylabel, filename, exclude_taut=False):
    plot_df = df.copy()

    if exclude_taut:
        plot_df = plot_df[plot_df["method_clean"] != "Our taut spline"]

    pivot = plot_df.pivot_table(
        index="method_clean",
        columns="dataset_group",
        values=metric,
        aggfunc="mean"
    )

    pivot = pivot.sort_index()

    ax = pivot.plot(kind="bar", figsize=(13, 7))

    plt.title(title)
    plt.xlabel("Method")
    plt.ylabel(ylabel)
    plt.xticks(rotation=35, ha="right")
    plt.grid(axis="y")
    plt.legend(title="Dataset")
    plt.tight_layout()

    path = PLOTS_DIR / filename
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()

    print("Saved:", path)


def plot_heatmap(df, metric, title, filename, exclude_taut=False):
    plot_df = df.copy()

    if exclude_taut:
        plot_df = plot_df[plot_df["method_clean"] != "Our taut spline"]

    pivot = plot_df.pivot_table(
        index="dataset_group",
        columns="method_clean",
        values=metric,
        aggfunc="mean"
    )

    fig, ax = plt.subplots(figsize=(13, 5))

    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    data = pivot.to_numpy(dtype=float)

    im = ax.imshow(data, aspect="auto", cmap="Blues")

    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_yticks(np.arange(len(pivot.index)))

    ax.set_xticklabels(pivot.columns, rotation=35, ha="right")
    ax.set_yticklabels(pivot.index)

    plt.colorbar(im, ax=ax, label=metric)

    ax.set_title(title)

    for i in range(len(pivot.index)):
        for j in range(len(pivot.columns)):
            value = data[i, j]
            if np.isfinite(value):
                ax.text(
                    j,
                    i,
                    f"{value:.3g}",
                    ha="center",
                    va="center",
                    fontsize=8
                )

    plt.tight_layout()

    path = PLOTS_DIR / filename
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()

    print("Saved:", path)


def create_best_method_table(df):
    stable_df = df[df["method_clean"] != "Our taut spline"].copy()

    best_mae = (
        stable_df.sort_values("mae")
        .groupby("dataset_group")
        .first()
        .reset_index()
    )

    best_rmse = (
        stable_df.sort_values("rmse")
        .groupby("dataset_group")
        .first()
        .reset_index()
    )

    best_relative_mae = (
        stable_df.sort_values("relative_mae")
        .groupby("dataset_group")
        .first()
        .reset_index()
    )

    best_runtime = (
        stable_df.sort_values("runtime_seconds")
        .groupby("dataset_group")
        .first()
        .reset_index()
    )

    best_summary = pd.DataFrame({
        "dataset": best_mae["dataset_group"],
        "best_by_mae": best_mae["method_clean"],
        "best_mae": best_mae["mae"],
        "best_by_rmse": best_rmse["method_clean"],
        "best_rmse": best_rmse["rmse"],
        "best_by_relative_mae": best_relative_mae["method_clean"],
        "best_relative_mae": best_relative_mae["relative_mae"],
        "fastest_method": best_runtime["method_clean"],
        "fastest_runtime": best_runtime["runtime_seconds"],
    })

    best_summary.to_csv(BEST_METHOD_CSV, index=False)

    print("\nBest method summary:")
    print(best_summary)
    print("Saved:", BEST_METHOD_CSV)


def main():
    df = load_all_results()

    save_summary_table(df)
    create_best_method_table(df)

    # ============================================================
    # Raw metric plots
    # ============================================================

    plot_grouped_bar(
        df,
        metric="mae",
        title="Unified Benchmark: MAE Across All Methods",
        ylabel="Mean Absolute Error",
        filename="unified_mae_all_methods.png",
        exclude_taut=False,
    )

    plot_grouped_bar(
        df,
        metric="rmse",
        title="Unified Benchmark: RMSE Across All Methods",
        ylabel="Root Mean Squared Error",
        filename="unified_rmse_all_methods.png",
        exclude_taut=False,
    )

    plot_grouped_bar(
        df,
        metric="mae",
        title="Unified Benchmark: MAE Across Stable Methods",
        ylabel="Mean Absolute Error",
        filename="unified_mae_stable_methods.png",
        exclude_taut=True,
    )

    plot_grouped_bar(
        df,
        metric="rmse",
        title="Unified Benchmark: RMSE Across Stable Methods",
        ylabel="Root Mean Squared Error",
        filename="unified_rmse_stable_methods.png",
        exclude_taut=True,
    )



    plot_grouped_bar(
        df,
        metric="relative_mae",
        title="Unified Benchmark: Relative MAE Across Stable Methods",
        ylabel="Relative MAE",
        filename="unified_relative_mae_stable_methods.png",
        exclude_taut=True,
    )

    plot_grouped_bar(
        df,
        metric="relative_rmse",
        title="Unified Benchmark: Relative RMSE Across Stable Methods",
        ylabel="Relative RMSE",
        filename="unified_relative_rmse_stable_methods.png",
        exclude_taut=True,
    )

    # ============================================================
    # Runtime plots
    # ============================================================

    plot_grouped_bar(
        df,
        metric="runtime_seconds",
        title="Unified Benchmark: Runtime Across All Methods",
        ylabel="Runtime seconds",
        filename="unified_runtime_all_methods.png",
        exclude_taut=False,
    )

    plot_grouped_bar(
        df,
        metric="runtime_seconds",
        title="Unified Benchmark: Runtime Across Stable Methods",
        ylabel="Runtime seconds",
        filename="unified_runtime_stable_methods.png",
        exclude_taut=True,
    )

    # ============================================================
    # Stability metric
    # ============================================================

    plot_grouped_bar(
        df,
        metric="stability_metric",
        title="Unified Benchmark: Dataset-Specific Stability Metric",
        ylabel="Stability metric value",
        filename="unified_stability_metric.png",
        exclude_taut=True,
    )

    # ============================================================
    # Heatmaps
    # ============================================================

    plot_heatmap(
        df,
        metric="mae",
        title="Unified MAE Heatmap Across Datasets and Methods",
        filename="unified_mae_heatmap_stable_methods.png",
        exclude_taut=True,
    )

    plot_heatmap(
        df,
        metric="rmse",
        title="Unified RMSE Heatmap Across Datasets and Methods",
        filename="unified_rmse_heatmap_stable_methods.png",
        exclude_taut=True,
    )

    plot_heatmap(
        df,
        metric="relative_mae",
        title="Unified Relative MAE Heatmap Across Datasets and Methods",
        filename="unified_relative_mae_heatmap_stable_methods.png",
        exclude_taut=True,
    )

    plot_heatmap(
        df,
        metric="relative_rmse",
        title="Unified Relative RMSE Heatmap Across Datasets and Methods",
        filename="unified_relative_rmse_heatmap_stable_methods.png",
        exclude_taut=True,
    )

    print("\nDone. Unified outputs saved in:", OUT_DIR)


if __name__ == "__main__":
    main()