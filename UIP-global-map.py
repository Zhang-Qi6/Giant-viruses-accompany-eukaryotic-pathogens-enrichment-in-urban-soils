#!/usr/bin/env python3

"""
Plot global UIP abundance on an Anthromes background raster.

This script overlays sampling points onto a global Anthromes raster map.
Point size represents UIP abundance, and the background raster represents
human-modified landscape categories.

Requirements:
    pandas
    numpy
    matplotlib
    rasterio
    openpyxl
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import rasterio


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot global UIP abundance on an Anthromes raster background."
    )

    parser.add_argument(
        "-r", "--raster",
        required=True,
        help="Input Anthromes raster file in GeoTIFF format."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input table containing longitude, latitude, and UIP abundance."
    )

    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Output directory for figures and source data."
    )

    parser.add_argument(
        "--lon-col",
        default="Longitude",
        help="Longitude column name. Default: Longitude"
    )

    parser.add_argument(
        "--lat-col",
        default="Latitude",
        help="Latitude column name. Default: Latitude"
    )

    parser.add_argument(
        "--value-col",
        default="UIP",
        help="UIP abundance column name. Default: UIP"
    )

    parser.add_argument(
        "--title",
        default="Global Distribution of Urban Eukaryotic Pathogens",
        help="Figure title."
    )

    parser.add_argument(
        "--scale-factor",
        type=float,
        default=9.0,
        help="Point-size scaling factor. Default: 9.0"
    )

    parser.add_argument(
        "--xlim",
        nargs=2,
        type=float,
        default=[-160, 180],
        help="Longitude range for plotting. Default: -160 180"
    )

    parser.add_argument(
        "--ylim",
        nargs=2,
        type=float,
        default=[-60, 85],
        help="Latitude range for plotting. Default: -60 85"
    )

    return parser.parse_args()


def load_table(input_file):
    input_file = Path(input_file)

    if input_file.suffix.lower() in [".xlsx", ".xls"]:
        return pd.read_excel(input_file)

    if input_file.suffix.lower() in [".csv"]:
        return pd.read_csv(input_file)

    if input_file.suffix.lower() in [".tsv", ".txt"]:
        return pd.read_csv(input_file, sep="\t")

    raise ValueError("Input table must be XLSX, XLS, CSV, TSV, or TXT.")


def format_longitude(x, _):
    if x < 0:
        return f"{abs(int(x))}°W"
    if x > 0:
        return f"{int(x)}°E"
    return "0°"


def format_latitude(y, _):
    if y < 0:
        return f"{abs(int(y))}°S"
    if y > 0:
        return f"{int(y)}°N"
    return "0°"


def classify_anthrome_value(value, nodata=None):
    if nodata is not None and value == nodata:
        return "Unknown"
    if pd.isna(value) or value == 0:
        return "Unknown"
    if value < 30:
        return "High modification"
    if value < 50:
        return "Medium modification"
    return "Low modification"


def extract_raster_values(df, raster_path, lon_col, lat_col):
    coordinates = list(zip(df[lon_col], df[lat_col]))

    with rasterio.open(raster_path) as src:
        nodata = src.nodata
        values = [sample[0] for sample in src.sample(coordinates)]

    df["Anthrome_value"] = values
    df["Anthrome_class"] = [
        classify_anthrome_value(value, nodata=nodata)
        for value in values
    ]

    return df


def add_size_legend(ax, values, scale_factor):
    max_value = np.nanmax(values)

    legend_values = [10, 100, 500, 1500]
    legend_values = [v for v in legend_values if v <= max_value]

    if len(legend_values) == 0:
        legend_values = [max_value]
    elif max_value > legend_values[-1]:
        legend_values.append(int(max_value))

    handles = []
    labels = []

    for value in legend_values:
        handle = ax.scatter(
            [],
            [],
            s=np.sqrt(value) * scale_factor,
            color="#E74C3C",
            alpha=0.75,
            edgecolors="white",
            linewidth=0.5,
        )
        handles.append(handle)
        labels.append(str(value))

    legend = ax.legend(
        handles,
        labels,
        title="UIP abundance\nMedian TPM",
        loc="lower left",
        bbox_to_anchor=(0.03, 0.05),
        frameon=True,
        facecolor="white",
        framealpha=0.9,
        edgecolor="#BDBDBD",
        labelspacing=1.5,
        borderpad=1.2,
    )

    legend.get_title().set_fontweight("bold")


def plot_map(df, raster_path, output_dir, args):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["axes.unicode_minus"] = False

    fig, ax = plt.subplots(figsize=(16, 9), facecolor="white")
    ax.set_facecolor("#DCE6F1")

    with rasterio.open(raster_path) as src:
        raster_data = src.read(1).astype(float)

        if src.nodata is not None:
            raster_data[raster_data == src.nodata] = np.nan

        raster_data[raster_data == 0] = np.nan

        extent = [
            src.bounds.left,
            src.bounds.right,
            src.bounds.bottom,
            src.bounds.top,
        ]

        background = ax.imshow(
            raster_data,
            extent=extent,
            cmap="Greys_r",
            alpha=0.6,
            vmin=10,
            vmax=65,
            zorder=1,
            aspect="auto",
        )

    point_sizes = np.sqrt(df[args.value_col]) * args.scale_factor

    ax.scatter(
        df[args.lon_col],
        df[args.lat_col],
        s=point_sizes,
        color="#E74C3C",
        alpha=0.75,
        edgecolors="white",
        linewidth=0.5,
        zorder=5,
    )

    ax.set_xlim(args.xlim)
    ax.set_ylim(args.ylim)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(60))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(30))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_longitude))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(format_latitude))

    ax.tick_params(axis="both", labelsize=11, colors="#34495E", length=5, width=1.2)

    for spine in ax.spines.values():
        spine.set_color("#BDBDBD")
        spine.set_linewidth(1.2)

    ax.grid(
        True,
        color="white",
        linestyle="-",
        linewidth=0.8,
        alpha=0.5,
        zorder=2,
    )

    colorbar = plt.colorbar(
        background,
        ax=ax,
        orientation="horizontal",
        shrink=0.52,
        pad=0.06,
        aspect=35,
    )

    colorbar.set_label(
        "Background gradient: human modification Anthromes",
        fontsize=12,
        fontweight="bold",
        color="#2C3E50",
        labelpad=10,
    )

    colorbar.set_ticks([15, 35, 60])
    colorbar.set_ticklabels([
        "High modification",
        "Medium modification",
        "Low modification",
    ])
    colorbar.ax.tick_params(labelsize=10, colors="#34495E")
    colorbar.outline.set_edgecolor("#BDBDBD")

    add_size_legend(ax, df[args.value_col].values, args.scale_factor)

    fig.suptitle(
        args.title,
        fontsize=18,
        fontweight="bold",
        y=0.95,
        color="#2C3E50",
    )

    plt.tight_layout()

    fig.savefig(output_dir / "global_uip_anthrome_map.pdf")
    fig.savefig(output_dir / "global_uip_anthrome_map.png", dpi=300)
    plt.close(fig)


def main():
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = load_table(args.input)

    required_cols = [args.lon_col, args.lat_col, args.value_col]
    missing_cols = [col for col in required_cols if col not in df.columns]

    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    df = df.dropna(subset=required_cols).copy()
    df = df.sort_values(args.value_col, ascending=False)

    df = extract_raster_values(
        df=df,
        raster_path=args.raster,
        lon_col=args.lon_col,
        lat_col=args.lat_col,
    )

    source_data_path = output_dir / "global_uip_map_source_data.tsv"
    df.to_csv(source_data_path, sep="\t", index=False)

    plot_map(
        df=df,
        raster_path=args.raster,
        output_dir=output_dir,
        args=args,
    )

    print("Global UIP Anthromes map completed.")
    print(f"Valid sampling points: {len(df)}")
    print(f"Source data: {source_data_path}")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()