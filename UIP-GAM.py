#!/usr/bin/env python3

"""
Ecological driver analysis for UIP abundance.

This script performs three complementary analyses:
1. Variance partitioning analysis (VPA)
2. PLS-VIP variable importance analysis
3. GAM-based nonlinear response analysis

Input:
    A metadata table in XLSX or CSV format containing:
    - one response column, e.g., UIP_abundance
    - multiple environmental, climatic, soil, or socioeconomic predictors

Outputs:
    - VPA source data
    - PLS-VIP scores
    - GAM curve and scatter source data
    - PDF and PNG figures

Requirements:
    pandas, numpy, matplotlib, scikit-learn, pygam, matplotlib-venn, openpyxl
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from pygam import LinearGAM, s
from sklearn.cross_decomposition import PLSRegression
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze ecological drivers of UIP abundance using VPA, PLS-VIP, and GAM."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input metadata table in XLSX or CSV format."
    )

    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Directory for output figures and source data."
    )

    parser.add_argument(
        "-y", "--target",
        default="UIP_abundance",
        help="Response variable column name. Default: UIP_abundance"
    )

    parser.add_argument(
        "--drop-cols",
        default="name,longitude,latitude",
        help="Comma-separated columns to exclude from predictors. Default: name,longitude,latitude"
    )

    parser.add_argument(
        "--gam-top-n",
        type=int,
        default=4,
        help="Number of variables shown in GAM plots. Default: 4"
    )

    parser.add_argument(
        "--pls-components",
        type=int,
        default=2,
        help="Number of PLS components. Default: 2"
    )

    return parser.parse_args()


def clean_column_names(columns):
    cleaned = (
        pd.Series(columns)
        .astype(str)
        .str.strip()
        .str.replace(r"[^a-zA-Z0-9_]", "_", regex=True)
        .str.replace(r"_+", "_", regex=True)
        .str.strip("_")
    )
    return cleaned.tolist()


def load_table(input_file):
    input_file = Path(input_file)

    if input_file.suffix.lower() in [".xlsx", ".xls"]:
        df = pd.read_excel(input_file)
    elif input_file.suffix.lower() in [".csv", ".txt", ".tsv"]:
        sep = "\t" if input_file.suffix.lower() in [".tsv", ".txt"] else ","
        df = pd.read_csv(input_file, sep=sep)
    else:
        raise ValueError("Input file must be XLSX, XLS, CSV, TSV, or TXT.")

    df.columns = clean_column_names(df.columns)
    return df


def assign_variable_groups(feature_names):
    socioeconomic_keywords = [
        "pop", "gdp", "gni", "pesticide", "fertilizer",
        "expansion", "industry", "energy", "co2"
    ]

    climate_keywords = [
        "temp", "precip", "evapor", "wind", "solar",
        "pressure", "climate"
    ]

    socioeconomic = []
    climate = []
    soil = []

    for feature in feature_names:
        name = feature.lower()
        if any(keyword in name for keyword in socioeconomic_keywords):
            socioeconomic.append(feature)
        elif any(keyword in name for keyword in climate_keywords):
            climate.append(feature)
        else:
            soil.append(feature)

    return socioeconomic, climate, soil


def adjusted_r2(X, y):
    if X.shape[1] == 0:
        return 0.0

    model = LinearRegression()
    model.fit(X, y)

    r2 = r2_score(y, model.predict(X))
    n, p = X.shape

    if n <= p + 1:
        return max(0.0, r2)

    adj_r2 = 1 - (1 - r2) * (n - 1) / (n - p - 1)
    return max(0.0, adj_r2)


def run_vpa(X_scaled, y_scaled, groups, output_dir):
    socioeconomic, climate, soil = groups

    r1 = adjusted_r2(X_scaled[socioeconomic], y_scaled)
    r2 = adjusted_r2(X_scaled[climate], y_scaled)
    r3 = adjusted_r2(X_scaled[soil], y_scaled)

    r12 = adjusted_r2(X_scaled[socioeconomic + climate], y_scaled)
    r13 = adjusted_r2(X_scaled[socioeconomic + soil], y_scaled)
    r23 = adjusted_r2(X_scaled[climate + soil], y_scaled)

    r123 = adjusted_r2(X_scaled[socioeconomic + climate + soil], y_scaled)

    pure_socioeconomic = r123 - r23
    pure_climate = r123 - r13
    pure_soil = r123 - r12

    shared_socio_climate = r12 - pure_socioeconomic - pure_climate
    shared_socio_soil = r13 - pure_socioeconomic - pure_soil
    shared_climate_soil = r23 - pure_climate - pure_soil

    shared_all = (
        r123
        - pure_socioeconomic
        - pure_climate
        - pure_soil
        - shared_socio_climate
        - shared_socio_soil
        - shared_climate_soil
    )

    fractions = [
        pure_socioeconomic,
        pure_climate,
        shared_socio_climate,
        pure_soil,
        shared_socio_soil,
        shared_climate_soil,
        shared_all,
    ]

    fractions = [max(0.0, value) for value in fractions]
    unexplained = max(0.0, 1.0 - r123)

    source_data = pd.DataFrame({
        "Variance_component": [
            "Pure_socioeconomic",
            "Pure_climate",
            "Shared_socioeconomic_climate",
            "Pure_soil",
            "Shared_socioeconomic_soil",
            "Shared_climate_soil",
            "Shared_all_three",
            "Unexplained",
        ],
        "Explained_percentage": [value * 100 for value in fractions] + [unexplained * 100],
    })

    source_data.to_csv(output_dir / "vpa_source_data.tsv", sep="\t", index=False)

    labels = [f"{value * 100:.1f}%" if value >= 0.001 else "<0.1%" for value in fractions]

    plt.figure(figsize=(8, 7))
    plt.title("Variance Partitioning Analysis of UIP Drivers", fontsize=14, fontweight="bold")

    v = venn3(
        subsets=(1, 1, 1, 1, 1, 1, 1),
        set_labels=("Socioeconomic", "Climate", "Soil properties")
    )

    subset_ids = ["100", "010", "110", "001", "101", "011", "111"]
    for label_text, subset_id in zip(labels, subset_ids):
        label = v.get_label_by_id(subset_id)
        if label:
            label.set_text(label_text)
            label.set_fontsize(12)
            label.set_fontweight("bold")

    plt.text(
        0.5,
        0.05,
        f"Unexplained\n{unexplained * 100:.1f}%",
        fontsize=12,
        fontweight="bold",
        ha="center",
        va="center",
        transform=plt.gca().transAxes,
    )

    plt.tight_layout()
    plt.savefig(output_dir / "vpa_venn_diagram.pdf")
    plt.savefig(output_dir / "vpa_venn_diagram.png", dpi=300)
    plt.close()


def calculate_vip(pls_model):
    scores = pls_model.x_scores_
    weights = pls_model.x_weights_
    loadings = pls_model.y_loadings_

    n_features, n_components = weights.shape
    vip_scores = np.zeros(n_features)

    explained_y = np.diag(scores.T @ scores @ loadings.T @ loadings).reshape(n_components, -1)
    total_explained_y = np.sum(explained_y)

    for i in range(n_features):
        weight = np.array([
            (weights[i, j] / np.linalg.norm(weights[:, j])) ** 2
            for j in range(n_components)
        ])
        vip_scores[i] = float(np.sqrt(n_features * np.dot(explained_y.T, weight) / total_explained_y))

    return vip_scores


def run_pls_vip(X_scaled, y_scaled, n_components, output_dir):
    n_components = min(n_components, X_scaled.shape[1])

    pls = PLSRegression(n_components=n_components)
    pls.fit(X_scaled, y_scaled)

    vip_scores = calculate_vip(pls)

    vip_df = pd.DataFrame({
        "Feature": X_scaled.columns,
        "VIP_score": vip_scores,
    })

    vip_df["Display_name"] = vip_df["Feature"].str.replace("_", " ", regex=False)
    vip_df = vip_df.sort_values("VIP_score", ascending=False)

    vip_df.to_csv(output_dir / "pls_vip_scores.tsv", sep="\t", index=False)

    plot_df = vip_df.sort_values("VIP_score", ascending=True)

    plt.figure(figsize=(9, max(5, len(plot_df) * 0.25)))
    bar_colors = np.where(plot_df["VIP_score"] >= 1.0, "#E74C3C", "#BDBDBD")

    plt.barh(plot_df["Display_name"], plot_df["VIP_score"], color=bar_colors)
    plt.axvline(x=1.0, linestyle="--", linewidth=1.2, color="black", label="VIP = 1.0")
    plt.xlabel("Variable Importance in Projection (VIP score)")
    plt.title("Driver Importance for UIP Abundance")
    plt.legend(loc="lower right")
    plt.tight_layout()

    plt.savefig(output_dir / "pls_vip_scores.pdf")
    plt.savefig(output_dir / "pls_vip_scores.png", dpi=300)
    plt.close()


def select_gam_features(X_imputed, y, top_n):
    corr_df = X_imputed.copy()
    corr_df["target"] = y

    ranked_features = (
        corr_df
        .corr(method="spearman")["target"]
        .abs()
        .drop("target")
        .sort_values(ascending=False)
        .index
        .tolist()
    )

    selected = []

    population_feature = next((col for col in X_imputed.columns if "pop" in col.lower()), None)
    gdp_feature = next((col for col in X_imputed.columns if "gdp" in col.lower()), None)

    if population_feature:
        selected.append(population_feature)

    if gdp_feature and gdp_feature not in selected:
        selected.append(gdp_feature)

    for feature in ranked_features:
        if len(selected) >= top_n:
            break
        if feature not in selected:
            selected.append(feature)

    return selected[:top_n]


def run_gam(X_imputed, y, target_features, output_dir):
    n_features = len(target_features)
    n_cols = 2
    n_rows = int(np.ceil(n_features / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(10, 4.5 * n_rows))
    axes = np.array(axes).reshape(-1)

    curve_data = []
    scatter_data = []

    for i, feature in enumerate(target_features):
        ax = axes[i]

        X_feature = X_imputed[[feature]].values
        gam = LinearGAM(s(0, n_splines=10)).fit(X_feature, y)

        grid = gam.generate_X_grid(term=0)
        partial_dependence, confidence = gam.partial_dependence(term=0, X=grid, width=0.95)

        y_pred = gam.predict(X_feature)
        r2_value = max(0.0, r2_score(y, y_pred))
        p_value = gam.statistics_["p_values"][0]

        partial_residual = y - gam.coef_[-1]

        ax.fill_between(grid[:, 0], confidence[:, 0], confidence[:, 1], alpha=0.2)
        ax.plot(grid[:, 0], partial_dependence, linewidth=2.0)
        ax.scatter(X_feature[:, 0], partial_residual, alpha=0.35, s=20)

        p_text = "P < 0.001" if p_value < 0.001 else f"P = {p_value:.3f}"
        display_name = feature.replace("_", " ")

        ax.set_title(f"UIP vs {display_name}\nR² = {r2_value:.3f}, {p_text}", fontsize=10)
        ax.set_xlabel(display_name)
        ax.set_ylabel("Partial effect on UIP abundance")

        curve_data.append(pd.DataFrame({
            "Feature": feature,
            "X_grid_value": grid[:, 0],
            "GAM_fitted_y": partial_dependence,
            "Confidence_lower_95": confidence[:, 0],
            "Confidence_upper_95": confidence[:, 1],
        }))

        scatter_data.append(pd.DataFrame({
            "Feature": feature,
            "Original_x_value": X_feature[:, 0],
            "Partial_residual_y": np.asarray(partial_residual),
        }))

    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.tight_layout()
    plt.savefig(output_dir / "gam_ecological_responses.pdf")
    plt.savefig(output_dir / "gam_ecological_responses.png", dpi=300)
    plt.close()

    pd.concat(curve_data).to_csv(output_dir / "gam_source_data_curves.tsv", sep="\t", index=False)
    pd.concat(scatter_data).to_csv(output_dir / "gam_source_data_scatters.tsv", sep="\t", index=False)


def main():
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = load_table(args.input)

    if args.target not in df.columns:
        raise ValueError(f"Target column not found: {args.target}")

    drop_cols = [col.strip().lower() for col in args.drop_cols.split(",")]

    excluded_cols = [
        col for col in df.columns
        if col.lower() in drop_cols or col == args.target
    ]

    features = [col for col in df.columns if col not in excluded_cols]

    X = df[features]
    y = df[args.target]

    imputer = SimpleImputer(strategy="median")
    scaler = StandardScaler()

    X_imputed = pd.DataFrame(imputer.fit_transform(X), columns=features)
    X_scaled = pd.DataFrame(scaler.fit_transform(X_imputed), columns=features)
    y_scaled = (y - y.mean()) / y.std()

    groups = assign_variable_groups(features)

    print("Variable groups:")
    print(f"  Socioeconomic: {len(groups[0])}")
    print(f"  Climate:       {len(groups[1])}")
    print(f"  Soil:          {len(groups[2])}")

    print("Running VPA...")
    run_vpa(X_scaled, y_scaled, groups, output_dir)

    print("Running PLS-VIP...")
    run_pls_vip(X_scaled, y_scaled, args.pls_components, output_dir)

    print("Running GAM...")
    target_features = select_gam_features(X_imputed, y, args.gam_top_n)
    run_gam(X_imputed, y, target_features, output_dir)

    print("Analysis completed.")
    print(f"Results saved to: {output_dir}")


if __name__ == "__main__":
    main()