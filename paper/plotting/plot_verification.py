#!/usr/bin/env python3
"""Create quadrature, adaptive-branch, and splitting verification figures.

By default, read CSV output from ``paper/verification/build/`` and write PDFs
to ``paper/plotting/`` in the QMOMKit repository.

Example:
    python3 paper/plotting/plot_verification.py
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "qmomkit-paper-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "qmomkit-paper-cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext, LogLocator, NullFormatter


FAMILY_STYLE = {
    "legendre_unit_weight_minus1_to_1": ("Gauss-Legendre", "#0072B2", "o"),
    "laguerre_alpha0_exp_minus_x": ("Gauss-Laguerre", "#D55E00", "s"),
    "hermite_probabilists_normalized_standard_normal": ("Gauss-Hermite", "#009E73", "^"),
}
SPLIT_STYLE = {
    "advsrc": (r"$\Phi_{\mathrm{src}}^{\Delta t} \circ \Phi_{\mathrm{adv}}^{\Delta t}$", "#0072B2", "o"),
    "srcadv": (r"$\Phi_{\mathrm{adv}}^{\Delta t} \circ \Phi_{\mathrm{src}}^{\Delta t}$", "#D55E00", "s"),
    "strang": (r"$\Phi_{\mathrm{adv}}^{\Delta t/2} \circ \Phi_{\mathrm{src}}^{\Delta t} \circ \Phi_{\mathrm{adv}}^{\Delta t/2}$", "#009E73", "^"),
}
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIRECTORY = Path(__file__).resolve().parent
VERIFICATION_BUILD_DIRECTORY = REPOSITORY_ROOT / "paper" / "verification" / "build"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def numeric(rows: Iterable[dict[str, str]], key: str) -> list[float]:
    return [float(row[key]) for row in rows]


def positive_or_nan(value: float) -> float:
    """Leave exact zeros unplotted on a logarithmic axis without changing data."""
    return value if value > 0.0 else math.nan


def set_paper_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9,
            "legend.fontsize": 7.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "lines.linewidth": 1.3,
            "lines.markersize": 4.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save(figure: plt.Figure, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, bbox_inches="tight")
    plt.close(figure)
    print(output)


def plot_quadrature(rows: list[dict[str, str]], output: Path) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(8.05, 2.25), sharex=True)
    quantities = (
        ("max_scaled_node_error", "Abscissa error"),
        ("max_scaled_weight_error", "Weight error"),
        ("max_scaled_moment_error", "Moment error"),
    )

    for family, (label, color, marker) in FAMILY_STYLE.items():
        family_rows = sorted((row for row in rows if row["family"] == family), key=lambda row: int(row["node_count"]))
        counts = [int(row["node_count"]) for row in family_rows]
        for axis, (key, title) in zip(axes, quantities):
            values = [positive_or_nan(float(row[key])) for row in family_rows]
            axis.plot(counts, values, color=color, marker=marker, label=label)
            axis.set_yscale("log")
            axis.set_title(title)
            axis.grid(True, which="both", alpha=0.25)
            axis.set_xticks(range(1, 11))
            axis.set_xlabel("Node count, $N$")

    axes[0].set_ylabel("Maximum scaled error, $E$")
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.98))
    figure.subplots_adjust(wspace=0.38, top=0.76)
    save(figure, output)


def plot_adaptive(rows: list[dict[str, str]], output: Path) -> None:
    rows = sorted(rows, key=lambda row: float(row["time"]))
    times = numeric(rows, "time")
    number_density = numeric(rows, "total_mass")
    mean = numeric(rows, "mean_velocity")
    active_nodes = [int(row["active_nodes"]) for row in rows]
    number_density_error = [abs(value - number_density[0]) for value in number_density]
    # Relaxation toward the equilibrium velocity 0.5 at rate 5, as in the study.
    exact_mean = [0.5 + (mean[0] - 0.5) * math.exp(-5.0 * time) for time in times]
    mean_error = [abs(value - reference) for value, reference in zip(mean, exact_mean)]

    figure, axes = plt.subplots(1, 2, figsize=(6, 2), sharex=True)
    axes[0].step(times, active_nodes, where="post", color="#0072B2")
    axes[0].set_ylim(0.7, 5.3)
    axes[0].set_yticks((1, 2, 3, 4, 5))
    axes[0].set_ylabel("Active velocity nodes")
    axes[0].set_xlabel("Time [s]")

    axes[1].plot(times, [positive_or_nan(value) for value in number_density_error], color="#0072B2", label="$M_{0,0}$")
    axes[1].plot(times, [positive_or_nan(value) for value in mean_error], color="#D55E00", label="$M_{0,1}/M_{0,0}$")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Absolute error")
    axes[1].set_xlabel("Time [s]")
    axes[1].legend(loc="lower left", frameon=False, fontsize=8)

    for axis in axes:
        axis.axvspan(4.0, 4.05, color="0.55", alpha=0.5, label="moment-source interval")
        axis.grid(True, which="both", alpha=0.25)
    figure.subplots_adjust(wspace=0.38)
    save(figure, output)


def grouped_rows(rows: list[dict[str, str]], key: str) -> dict[str, list[dict[str, str]]]:
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[row[key]].append(row)
    for group in groups.values():
        group.sort(key=lambda row: float(row["outer_step_size"]), reverse=True)
    return groups


def plot_splitting_convergence(rows: list[dict[str, str]], output: Path) -> None:
    figure, axis = plt.subplots(figsize=(3.5, 2.7))

    split_groups = grouped_rows(rows, "split_mode")
    for mode in ("advsrc", "srcadv", "strang"):
        rows = split_groups[mode]
        label, color, marker = SPLIT_STYLE[mode]
        steps = numeric(rows, "outer_step_size")
        errors = numeric(rows, "max_scaled_moment_error")
        axis.loglog(steps, errors, color=color, marker=marker, label=label)

    # Reference slopes for the expected convergence orders.
    x0, x1 = 2.0e-2, 4.0e-2
    x_mid = (x0 * x1)**0.5
    # First-order reference line
    y0 = 2.4e-3
    y1 = y0 * (x1 / x0)
    axis.loglog([x0, x1], [y0, y1], color="black", linewidth=0.8)
    y_mid = (y0 * y1)**0.5
    axis.annotate("1", xy=(x_mid, y_mid), xytext=(0, 9), textcoords="offset points", ha="center", va="top", fontsize=8)
    # Second-order reference line
    y0 = 1.4e-6
    y1 = y0 * (x1 / x0) ** 2
    axis.loglog([x0, x1], [y0, y1], color="black", linewidth=0.8)
    y_mid = (y0 * y1)**0.5
    axis.annotate("2", xy=(x_mid, y_mid), xytext=(0, 4), textcoords="offset points", ha="center", va="bottom", fontsize=8)

    axis.set_xlabel("$\\Delta t$ [s]")
    axis.set_ylabel("$E_{\\Delta t}$")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(frameon=False, loc="lower right", bbox_to_anchor=(0.52, 0.3), fontsize=8)
    axis.xaxis.set_major_locator(LogLocator(base=10, numticks=3))
    axis.xaxis.set_major_formatter(LogFormatterMathtext(base=10))
    axis.xaxis.set_minor_formatter(NullFormatter())
    save(figure, output)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create quadrature, adaptive-branch, and splitting figures from QMOMKit CSV output.")
    parser.add_argument(
        "--quadrature",
        type=Path,
        default=VERIFICATION_BUILD_DIRECTORY / "quadrature_1d.csv",
        help="Path to quadrature_1d.csv (default: paper/verification/build/quadrature_1d.csv).",
    )
    parser.add_argument(
        "--adaptive",
        type=Path,
        default=VERIFICATION_BUILD_DIRECTORY / "adaptive_branch_quadrature.csv",
        help="Path to adaptive_branch_quadrature.csv (default: paper/verification/build/adaptive_branch_quadrature.csv).",
    )
    parser.add_argument(
        "--splitting",
        type=Path,
        default=VERIFICATION_BUILD_DIRECTORY / "splitting_convergence.csv",
        help="Path to splitting_convergence.csv (default: paper/verification/build/splitting_convergence.csv).",
    )
    parser.add_argument("--output-dir", type=Path, default=SCRIPT_DIRECTORY, help="Directory for vector-PDF figures.")
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    set_paper_style()
    plot_quadrature(read_csv(arguments.quadrature), arguments.output_dir / "quadrature_reconstruction.pdf")
    plot_adaptive(read_csv(arguments.adaptive), arguments.output_dir / "adaptive_branch_quadrature.pdf")
    plot_splitting_convergence(read_csv(arguments.splitting), arguments.output_dir / "splitting_convergence.pdf")


if __name__ == "__main__":
    main()
