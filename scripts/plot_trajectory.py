#!/usr/bin/env python3
"""Visualisasi lintasan peluru dari results/trajectory.csv.

Membaca keluaran program Fortran (kolom t, x, y, vx, vy), lalu menggambar
lintasan dengan hambatan udara beserta lintasan vakum sebagai pembanding.
Grafik disimpan ke results/trajectory.png.

Pemakaian (dari folder mana pun):
    python scripts/plot_trajectory.py
    python scripts/plot_trajectory.py --input results/trajectory.csv \
        --output results/trajectory.png --no-vacuum
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # backend tanpa jendela, aman untuk terminal/CI

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = ROOT / "results" / "trajectory.csv"
DEFAULT_OUTPUT = ROOT / "results" / "trajectory.png"

# Samakan dengan nilai g pada src/ballistics.f90
G = 9.80665  # m/s^2

# Palet
COLOR_DRAG = "#2a78d6"     # lintasan dengan hambatan udara
COLOR_VACUUM = "#898781"   # lintasan vakum (acuan)
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
SURFACE = "#fcfcfb"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot lintasan peluru 2D.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT,
                        help="berkas CSV lintasan (bawaan: results/trajectory.csv)")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT,
                        help="berkas PNG keluaran (bawaan: results/trajectory.png)")
    parser.add_argument("--no-vacuum", action="store_true",
                        help="jangan gambar lintasan vakum pembanding")
    parser.add_argument("--dpi", type=int, default=150, help="resolusi PNG")
    return parser.parse_args()


def load_trajectory(path: Path) -> pd.DataFrame:
    if not path.is_file():
        sys.exit(f"Berkas tidak ditemukan: {path}\n"
                 "Jalankan simulasi Fortran terlebih dahulu (./ballistics).")
    df = pd.read_csv(path, skipinitialspace=True)
    missing = {"t", "x", "y", "vx", "vy"} - set(df.columns)
    if missing:
        sys.exit(f"Kolom hilang di {path}: {', '.join(sorted(missing))}")
    if len(df) < 2:
        sys.exit(f"Data lintasan terlalu sedikit di {path}")
    return df


def vacuum_trajectory(vx0: float, vy0: float, y0: float, n: int = 400):
    """Parabola vakum dari kondisi awal yang sama, hingga y = 0."""
    t_flight = (vy0 + np.sqrt(vy0**2 + 2.0 * G * y0)) / G
    t = np.linspace(0.0, t_flight, n)
    return vx0 * t, y0 + vy0 * t - 0.5 * G * t**2


def style_axes(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(BASELINE)
    ax.tick_params(colors=INK_MUTED, labelsize=9, length=0, pad=6)
    ax.grid(True, axis="y", color=GRID, linewidth=0.8)
    ax.grid(True, axis="x", color=GRID, linewidth=0.8, alpha=0.6)
    ax.set_axisbelow(True)


def main() -> None:
    args = parse_args()
    df = load_trajectory(args.input)

    x, y, t = df["x"].to_numpy(), df["y"].to_numpy(), df["t"].to_numpy()
    i_apex = int(np.argmax(y))
    x_range, h_max, t_flight = x[-1], y[i_apex], t[-1]

    fig, ax = plt.subplots(figsize=(10, 5.2), facecolor=SURFACE)
    style_axes(ax)

    if not args.no_vacuum:
        xv, yv = vacuum_trajectory(df["vx"].iloc[0], df["vy"].iloc[0], y[0])
        ax.plot(xv, yv, color=COLOR_VACUUM, linewidth=1.6, linestyle=(0, (5, 4)),
                label="Tanpa hambatan (vakum)", zorder=2)
        ax.annotate("vakum", xy=(xv[len(xv) // 2], yv.max()),
                    xytext=(0, 6), textcoords="offset points",
                    ha="center", va="bottom", fontsize=9, color=INK_SECONDARY)

    ax.plot(x, y, color=COLOR_DRAG, linewidth=2.2,
            label="Dengan hambatan udara (RK4)", zorder=3)

    # Titik puncak dan titik tumbukan, dengan cincin warna permukaan
    ax.scatter([x[i_apex], x_range], [h_max, y[-1]], s=64, color=COLOR_DRAG,
               edgecolors=SURFACE, linewidths=2, zorder=4)
    ax.annotate(f"H maks = {h_max:.2f} m", xy=(x[i_apex], h_max),
                xytext=(0, 10), textcoords="offset points",
                ha="center", va="bottom", fontsize=10, color=INK_PRIMARY)
    ax.annotate(f"R = {x_range:.2f} m\nT = {t_flight:.2f} s", xy=(x_range, 0),
                xytext=(10, 6), textcoords="offset points",
                ha="left", va="bottom", fontsize=10, color=INK_PRIMARY)

    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.margins(x=0.04, y=0.18)
    ax.set_xlabel("Jarak horizontal x [m]", color=INK_SECONDARY, fontsize=10)
    ax.set_ylabel("Ketinggian y [m]", color=INK_SECONDARY, fontsize=10)
    ax.set_title("Lintasan peluru 2D dengan hambatan udara kuadratik",
                 loc="left", fontsize=13, color=INK_PRIMARY, pad=14)

    legend = ax.legend(loc="upper right", frameon=False, fontsize=9)
    for text in legend.get_texts():
        text.set_color(INK_SECONDARY)

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, facecolor=SURFACE)
    plt.close(fig)

    print(f"Grafik disimpan ke {args.output}")
    print(f"  R = {x_range:.4f} m | H = {h_max:.4f} m | T = {t_flight:.4f} s")


if __name__ == "__main__":
    main()
