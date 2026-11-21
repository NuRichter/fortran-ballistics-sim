#!/usr/bin/env python3
"""Plot the projectile trajectory from results/trajectory.csv.

Reads the output of the Fortran program (columns t, x, y, vx, vy) and draws
the trajectory with air drag together with the vacuum trajectory for
comparison. The figure is saved to results/trajectory.png.

Usage (from any folder):
    python scripts/plot_trajectory.py
    python scripts/plot_trajectory.py --input results/trajectory.csv \
        --output results/trajectory.png --no-vacuum
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # non-interactive backend, safe for terminals and CI

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = ROOT / "results" / "trajectory.csv"
DEFAULT_OUTPUT = ROOT / "results" / "trajectory.png"

# Keep in sync with the value of g in src/ballistics.f90
G = 9.80665  # m/s^2

# Palette
COLOR_DRAG = "#2a78d6"     # trajectory with air drag
COLOR_VACUUM = "#898781"   # vacuum trajectory (reference)
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
SURFACE = "#fcfcfb"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot a 2D projectile trajectory.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT,
                        help="trajectory CSV file (default: results/trajectory.csv)")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT,
                        help="output PNG file (default: results/trajectory.png)")
    parser.add_argument("--no-vacuum", action="store_true",
                        help="do not draw the vacuum comparison trajectory")
    parser.add_argument("--dpi", type=int, default=150, help="PNG resolution")
    return parser.parse_args()


def load_trajectory(path: Path) -> pd.DataFrame:
    if not path.is_file():
        sys.exit(f"File not found: {path}\n"
                 "Run the Fortran simulation first (./ballistics).")
    df = pd.read_csv(path, skipinitialspace=True)
    missing = {"t", "x", "y", "vx", "vy"} - set(df.columns)
    if missing:
        sys.exit(f"Missing columns in {path}: {', '.join(sorted(missing))}")
    if len(df) < 2:
        sys.exit(f"Not enough trajectory data in {path}")
    return df


def vacuum_trajectory(vx0: float, vy0: float, y0: float, n: int = 400):
    """Vacuum parabola from the same initial conditions, down to y = 0."""
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
                label="No drag (vacuum)", zorder=2)
        ax.annotate("vacuum", xy=(xv[len(xv) // 2], yv.max()),
                    xytext=(0, 6), textcoords="offset points",
                    ha="center", va="bottom", fontsize=9, color=INK_SECONDARY)

    ax.plot(x, y, color=COLOR_DRAG, linewidth=2.2,
            label="With air drag (RK4)", zorder=3)

    # Apex and impact point, with a ring in the surface color
    ax.scatter([x[i_apex], x_range], [h_max, y[-1]], s=64, color=COLOR_DRAG,
               edgecolors=SURFACE, linewidths=2, zorder=4)
    ax.annotate(f"H max = {h_max:.2f} m", xy=(x[i_apex], h_max),
                xytext=(0, 10), textcoords="offset points",
                ha="center", va="bottom", fontsize=10, color=INK_PRIMARY)
    ax.annotate(f"R = {x_range:.2f} m\nT = {t_flight:.2f} s", xy=(x_range, 0),
                xytext=(10, 6), textcoords="offset points",
                ha="left", va="bottom", fontsize=10, color=INK_PRIMARY)

    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.margins(x=0.04, y=0.18)
    ax.set_xlabel("Horizontal distance x [m]", color=INK_SECONDARY, fontsize=10)
    ax.set_ylabel("Height y [m]", color=INK_SECONDARY, fontsize=10)
    ax.set_title("2D projectile trajectory with quadratic air drag",
                 loc="left", fontsize=13, color=INK_PRIMARY, pad=14)

    legend = ax.legend(loc="upper right", frameon=False, fontsize=9)
    for text in legend.get_texts():
        text.set_color(INK_SECONDARY)

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, facecolor=SURFACE)
    plt.close(fig)

    print(f"Figure saved to {args.output}")
    print(f"  R = {x_range:.4f} m | H = {h_max:.4f} m | T = {t_flight:.4f} s")


if __name__ == "__main__":
    main()
