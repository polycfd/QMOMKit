#!/usr/bin/env python3
"""Create the application-study figures from the QMOMKit paper-study output.

By default, read the study output from ``paper/examples/build/`` and write PDFs
to ``paper/plotting/`` in the QMOMKit repository.

Example:
    python3 paper/plotting/plot_examples.py
"""

from __future__ import annotations

import argparse
import math
import os
import tempfile
from collections.abc import Callable
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "qmomkit-paper-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "qmomkit-paper-cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


STYLE = {
    "monte_carlo": ("Monte Carlo", "#B3B3B3", 2.6, 1),
    "qmomkit": ("QMOMKit", "#0072B2", 1.15, 2),
}
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIRECTORY = Path(__file__).resolve().parent

Row = dict[str, float]
Comparison = list[dict[str, float]]


def read_history(path: Path) -> list[Row]:
    """Read a QMOMKit history file, one dictionary per row keyed by column name."""
    with path.open(encoding="utf-8") as handle:
        # QMOMKit centers the "# time" label inside a padded column, so the
        # comment marker is not at the start of the line.
        header = handle.readline().strip().lstrip("#").split()
        return [dict(zip(header, (float(value) for value in line.split()))) for line in handle if line.split()]


def compare(
    qmomkit_rows: list[Row], monte_carlo_rows: list[Row], quantities: dict[str, Callable[[Row], float]]
) -> dict[str, Comparison]:
    """Evaluate each compared quantity on the QMOMKit and Monte Carlo rows of the same output time."""
    return {
        name: [
            {"time": qmomkit_row["time"], "qmomkit": evaluate(qmomkit_row), "monte_carlo": evaluate(monte_carlo_row)}
            for qmomkit_row, monte_carlo_row in zip(qmomkit_rows, monte_carlo_rows)
        ]
        for name, evaluate in quantities.items()
    }


def column(name: str) -> Callable[[Row], float]:
    return lambda row: row[name]


def ratio(numerator: str, denominator: str, factor: float = 1.0) -> Callable[[Row], float]:
    return lambda row: factor * row[numerator] / row[denominator]


def cavitation_data(results_dir: Path) -> dict[str, Comparison]:
    """Number density and the population-normalized first, mixed and second moments in (R, U)."""
    return compare(
        read_history(results_dir / "moment_history.txt"),
        read_history(results_dir / "moment_history_mc.txt"),
        {
            "M_00": column("M[0,0]"),
            "M_10/M_00": ratio("M[1,0]", "M[0,0]"),
            "M_01/M_00": ratio("M[0,1]", "M[0,0]"),
            "M_11/M_00": ratio("M[1,1]", "M[0,0]"),
            "M_20/M_00": ratio("M[2,0]", "M[0,0]"),
            "M_02/M_00": ratio("M[0,2]", "M[0,0]"),
        },
    )


def condensing_droplet_data(results_dir: Path) -> dict[str, Comparison]:
    """Droplet quantities in (r, w, T). Radii are stored in meters; the micrometer factors are presentation only."""

    def radius_variance(row: Row) -> float:
        return 1.0e12 * (row["M[2,0,0]"] / row["M[0,0,0]"] - (row["M[1,0,0]"] / row["M[0,0,0]"]) ** 2)

    return compare(
        read_history(results_dir / "moment_history.txt"),
        read_history(results_dir / "moment_history_mc.txt"),
        {
            "M_000": column("M[0,0,0]"),
            "M_300/M_200": ratio("M[3,0,0]", "M[2,0,0]", 1.0e6),
            "M_010/M_000": ratio("M[0,1,0]", "M[0,0,0]"),
            "M_001/M_000": ratio("M[0,0,1]", "M[0,0,0]"),
            "sigma_r^2": radius_variance,
            "M_101/M_000": ratio("M[1,0,1]", "M[0,0,0]", 1.0e6),
        },
    )


def spherical_breakup_quadrature_4d_data(results_dir: Path) -> dict[str, Comparison]:
    """Tracked moments, and the untracked moments the configuration adds as generalized columns, in (d, u, v, w)."""
    return compare(
        read_history(results_dir / "moment_history.txt"),
        read_history(results_dir / "moment_history_mc.txt"),
        {
            "M_0000": column("M[0,0,0,0]"),
            "M_1000/M_0000": ratio("M[1,0,0,0]", "M[0,0,0,0]", 1.0e3),
            "M_3000": column("M[3,0,0,0]"),
            "M_3100": column("G[3,1,0,0]"),
            "M_3010": column("G[3,0,1,0]"),
            "M_3001": column("G[3,0,0,1]"),
        },
    )


def spherical_breakup_moments_5d_data(results_dir: Path) -> dict[str, Comparison]:
    """Tracked moments, and the generalized moments of volume order 4/3 the configuration adds, in (V, T, u, v, w)."""
    return compare(
        read_history(results_dir / "moment_history.txt"),
        read_history(results_dir / "moment_history_mc.txt"),
        {
            "M_00000": column("M[0,0,0,0,0]"),
            "M_10000": column("M[1,0,0,0,0]"),
            "M_10100": column("M[1,0,1,0,0]"),
            "M_11000": column("M[1,1,0,0,0]"),
            "M_4over3_0000": column("G[1.3333333333333333,0,0,0,0]"),
            "M_4over3_1000": column("G[1.3333333333333333,1,0,0,0]"),
        },
    )


def set_paper_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9,
            "legend.fontsize": 8,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "lines.linewidth": 1.35,
            "lines.markersize": 4.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def plot_comparison(axis: plt.Axes, rows: Comparison) -> None:
    time = [row["time"] for row in rows]
    for run in ("monte_carlo", "qmomkit"):
        label, color, linewidth, zorder = STYLE[run]
        axis.plot(time, [row[run] for row in rows], color=color, linewidth=linewidth, zorder=zorder, label=label)


def plot_panels(panels: tuple[tuple[plt.Axes, Comparison, str | None, str], ...]) -> None:
    for axis, rows, title, y_label in panels:
        plot_comparison(axis, rows)
        if title is not None:
            axis.set_title(title)
        axis.set_ylabel(y_label)
        axis.grid(True, alpha=0.25)


def add_legend(figure: plt.Figure, axis: plt.Axes, ncol: int, top: float) -> None:
    handles, labels = axis.get_legend_handles_labels()
    figure.legend(handles, labels, loc="upper center", ncol=ncol, frameon=False, bbox_to_anchor=(0.5, top))


def save_figure(figure: plt.Figure, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, bbox_inches="tight")
    plt.close(figure)
    print(output)


def plot_cavitation(data: dict[str, Comparison], output: Path) -> None:
    figure, axes = plt.subplots(2, 3, figsize=(7.05, 4.45), sharex=True)
    plot_panels(
        (
            (axes[0, 0], data["M_00"], "Number density", r"$M_{0,0}$"),
            (axes[0, 1], data["M_10/M_00"], "Mean bubble radius", r"$M_{1,0}$"),
            (axes[0, 2], data["M_01/M_00"], "Mean radial velocity", r"$M_{0,1}$"),
            (axes[1, 0], data["M_11/M_00"], "Radius-velocity moment", r"$M_{1,1}$"),
            (axes[1, 1], data["M_20/M_00"], "Second radius moment", r"$M_{2,0}$"),
            (axes[1, 2], data["M_02/M_00"], "Second velocity moment", r"$M_{0,2}$"),
        ),
    )
    for axis in axes[1, :]:
        axis.set_xlabel("Time")
    add_legend(figure, axes[0, 0], 2, 0.91)
    figure.subplots_adjust(left=0.10, right=0.985, bottom=0.12, top=0.8, hspace=0.28, wspace=0.43)
    save_figure(figure, output)


def plot_condensing_droplets(data: dict[str, Comparison], output: Path) -> None:
    figure, axes = plt.subplots(2, 3, figsize=(7.05, 3.7), sharex=True)
    # The shaded band is the nucleation interval.
    for axis in axes.flat:
        axis.axvspan(0.20, 0.40, color="#E69F00", alpha=0.12, zorder=0)
    plot_panels(
        (
            (axes[0, 0], data["M_000"], "Total number density", r"$M_{0,0,0}$ [-]"),
            (axes[0, 1], data["M_300/M_200"], "Characteristic radius", r"$M_{3,0,0}/M_{2,0,0}$ [$\mu$m]"),
            (axes[0, 2], data["M_010/M_000"], "Mean vertical velocity", r"$M_{0,1,0}/M_{0,0,0}$ [m/s]"),
            (axes[1, 0], data["M_001/M_000"], "Mean temperature", r"$M_{0,0,1}/M_{0,0,0}$ [K]"),
            (axes[1, 2], data["sigma_r^2"], "Radius variance", r"$\sigma_r^2$ [$\mu$m$^2$]"),
            (axes[1, 1], data["M_101/M_000"], "Radius-weighted mean temperature", r"$M_{1,0,1}/M_{0,0,0}$ [$\mu$m K]"),
        ),
    )
    for axis in axes[1, :]:
        axis.set_xlabel("Time [s]")
    add_legend(figure, axes[0, 0], 2, 1)
    figure.subplots_adjust(left=0.11, right=0.985, bottom=0.09, top=0.88, hspace=0.25, wspace=0.45)
    save_figure(figure, output)


def plot_spherical_breakup_quadrature_4d(data: dict[str, Comparison], output: Path) -> None:
    figure, axes = plt.subplots(2, 3, figsize=(7.05, 4.3), sharex=True)
    plot_panels(
        (
            (axes[0, 0], data["M_0000"], None, r"$M_{0,0,0,0}$ [m$^{-3}$]"),
            (axes[0, 1], data["M_1000/M_0000"], None, r"$M_{1,0,0,0}/M_{0,0,0,0}$ [mm]"),
            (axes[0, 2], data["M_3000"], None, r"$M_{3,0,0,0}$ [-]"),
            (axes[1, 0], data["M_3100"], None, r"$M_{3,1,0,0}$ [m/s]"),
            (axes[1, 1], data["M_3010"], None, r"$M_{3,0,1,0}$ [m/s]"),
            (axes[1, 2], data["M_3001"], None, r"$M_{3,0,0,1}$ [m/s]"),
        ),
    )
    reference_volume_moment = data["M_3000"][0]["qmomkit"]
    axes[0, 2].set_ylim(0.9 * reference_volume_moment, 1.1 * reference_volume_moment)
    for axis in axes[1, :]:
        axis.set_xlabel("Time [s]")
    add_legend(figure, axes[0, 0], 2, 0.88)
    figure.subplots_adjust(left=0.10, right=0.985, bottom=0.12, top=0.8, hspace=0.2, wspace=0.43)
    save_figure(figure, output)


def plot_spherical_breakup_moments_5d(data: dict[str, Comparison], output: Path) -> None:
    figure, axes = plt.subplots(2, 3, figsize=(7.05, 4.3), sharex=True)
    plot_panels(
        (
            (axes[0, 0], data["M_00000"], None, r"$M_{0,0,0,0,0}$ [m$^{-3}$]"),
            (axes[0, 1], data["M_10000"], None, r"$M_{1,0,0,0,0}$ [-]"),
            (axes[0, 2], data["M_10100"], None, r"$M_{1,0,1,0,0}$ [m/s]"),
            (axes[1, 0], data["M_11000"], None, r"$M_{1,1,0,0,0}$ [K]"),
            (axes[1, 1], data["M_4over3_0000"], None, r"$M_{4/3,0,0,0,0}$ [m]"),
            (axes[1, 2], data["M_4over3_1000"], None, r"$M_{4/3,1,0,0,0}$ [m K]"),
        ),
    )
    reference_volume_moment = data["M_10000"][0]["qmomkit"]
    axes[0, 1].set_ylim(0.9 * reference_volume_moment, 1.1 * reference_volume_moment)
    for axis in axes[1, :]:
        axis.set_xlabel("Time [s]")
    add_legend(figure, axes[0, 0], 2, 0.9)
    figure.subplots_adjust(left=0.10, right=0.985, bottom=0.12, top=0.8, hspace=0.2, wspace=0.43)
    save_figure(figure, output)


def read_statistics_history(path: Path) -> list[tuple[float, ...]]:
    """Read the rows "time M0 d32 diameter_variance" of an aerosol history."""
    with path.open(encoding="utf-8") as handle:
        return [tuple(float(value) for value in line.split()) for line in handle if line and not line.startswith("#")]


def plot_aerosol_evaporation(results_dir: Path, output: Path) -> None:
    mc_history = read_statistics_history(results_dir / "moment_history_mc.txt")
    qmom_histories = {
        node_count: read_statistics_history(results_dir / f"moment_history_qmom_n{node_count}.txt") for node_count in (2, 4, 6, 8, 10)
    }
    with (results_dir / "initial_mc_samples.txt").open(encoding="utf-8") as handle:
        initial_diameters = [float(line) for line in handle if line and not line.startswith("#")]

    figure, axes = plt.subplots(2, 2, figsize=(6.0, 4.5))

    # Number fraction per decade of the retained initial sample, on 60
    # logarithmic bins over the modeled diameter support [0.8, 1000] um.
    histogram_axis = axes[0, 0]
    log_minimum = math.log10(0.8)
    log_bin_width = (math.log10(1000.0) - log_minimum) / 60.0
    histogram_bins = [10.0 ** (log_minimum + bin_index * log_bin_width) for bin_index in range(61)]
    bin_counts = [0] * 60
    for diameter in initial_diameters:
        bin_counts[min(59, max(0, int((math.log10(diameter) - log_minimum) / log_bin_width)))] += 1
    nonzero_bins = [
        (
            math.sqrt(histogram_bins[index] * histogram_bins[index + 1]),
            bin_counts[index] / (len(initial_diameters) * log_bin_width),
            histogram_bins[index + 1] - histogram_bins[index],
        )
        for index in range(60)
        if bin_counts[index] > 0
    ]
    histogram_axis.bar(
        [bin_center for bin_center, _, _ in nonzero_bins],
        [fraction for _, fraction, _ in nonzero_bins],
        width=[bin_width for _, _, bin_width in nonzero_bins],
        align="center",
        color="#D9D9D9",
        edgecolor="#8C8C8C",
        linewidth=0.35,
        zorder=1,
    )
    histogram_axis.set_xscale("log")
    histogram_axis.set_yscale("log")
    histogram_axis.set_xlim(0.8, 1000.0)
    histogram_axis.set_ylim(bottom=0.7 * min(fraction for _, fraction, _ in nonzero_bins))
    histogram_axis.set_xlabel(r"$d$ [$\mu$m]")
    histogram_axis.set_ylabel("Normalized number distribution")
    histogram_axis.grid(True, which="both", alpha=0.25)

    qmom_styles = {
        2: ("#4D4D4D", ":", 1.20),
        4: ("#6A6A6A", "--", 1.20),
        6: ("#0072B2", "-.", 1.20),
        8: ("#005A9C", (0, (4, 1, 1, 1)), 1.20),
        10: ("#003B6F", "-", 1.35),
    }
    panels = (
        (axes[0, 1], 1, r"$M_0$ [$10^5$ m$^{-3}$]", 1.0e-5),
        (axes[1, 0], 2, r"$d_{32}$ [$\mu$m]", 1.0),
        (axes[1, 1], 3, r"$\sigma_d^2$ [$\mu$m$^2$]", 1.0),
    )
    for axis, statistic, y_label, scale in panels:
        axis.plot(
            [row[0] for row in mc_history],
            [scale * row[statistic] for row in mc_history],
            color="#B3B3B3",
            linewidth=2.6,
            zorder=1,
            label="Monte Carlo",
        )
        for node_count, history in qmom_histories.items():
            color, line_style, line_width = qmom_styles[node_count]
            axis.plot(
                [row[0] for row in history],
                [scale * row[statistic] for row in history],
                color=color,
                linestyle=line_style,
                linewidth=line_width,
                zorder=2,
                label=f"QMOMKit, $N={node_count}$",
            )
        axis.set_xlim(0.0, 20.0)
        axis.set_ylabel(y_label)
        if statistic == 1:
            axis.set_ylim(2.20, 2.26)
            axis.ticklabel_format(axis="y", style="plain", useOffset=False)
        axis.grid(True, alpha=0.25)

    axes[0, 1].set_xlabel("Time [s]")
    axes[1, 0].set_xlabel("Time [s]")
    axes[1, 1].set_xlabel("Time [s]")
    add_legend(figure, axes[0, 1], 3, 0.95)
    figure.subplots_adjust(left=0.12, right=0.985, bottom=0.13, top=0.84, wspace=0.4, hspace=0.34)
    save_figure(figure, output)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create paper application-study figures from QMOMKit output.")
    parser.add_argument(
        "--examples-build-dir",
        type=Path,
        default=REPOSITORY_ROOT / "paper" / "examples" / "build",
        help="Paper examples build directory containing the per-case output directories.",
    )
    parser.add_argument("--output-dir", type=Path, default=SCRIPT_DIRECTORY, help="Directory for vector-PDF figures.")
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    build = arguments.examples_build_dir
    figures = arguments.output_dir
    set_paper_style()

    plot_aerosol_evaporation(build / "1d_aerosols_evaporation" / "output", figures / "aerosol_evaporation_moments.pdf")
    plot_cavitation(cavitation_data(build / "2d_bubbles_cavitation" / "output"), figures / "cavitation_moments.pdf")
    plot_condensing_droplets(
        condensing_droplet_data(build / "3d_droplets_condensation-nucleation" / "output"), figures / "condensing_droplets_moments.pdf"
    )
    plot_spherical_breakup_quadrature_4d(
        spherical_breakup_quadrature_4d_data(build / "4d_particles_breakup_recoil" / "output"),
        figures / "cqmom_spherical_breakup_quadrature_4d_moments.pdf",
    )
    plot_spherical_breakup_moments_5d(
        spherical_breakup_moments_5d_data(build / "5d_particles_breakup_heat-transfer" / "output"),
        figures / "cqmom_spherical_breakup_moments_5d_moments.pdf",
    )


if __name__ == "__main__":
    main()
