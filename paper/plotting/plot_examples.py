#!/usr/bin/env python3
"""Create the application-example figures from QMOMKit paper-case output."""

from __future__ import annotations

import argparse
import math
import os
import tempfile
from collections.abc import Callable
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "qmomkit-paper-matplotlib"),
)
os.environ.setdefault(
    "XDG_CACHE_HOME",
    str(Path(tempfile.gettempdir()) / "qmomkit-paper-cache"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


STYLE = {
    "monte_carlo": ("Monte Carlo", "#B3B3B3", 2.6, 1),
    "qmomkit": ("QMOMKit", "#0072B2", 1.15, 2),
}


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIRECTORY = Path(__file__).resolve().parent


def read_moment_rows(
    path: Path,
    labels: list[str],
) -> list[dict[str, float]]:
    """Read selected raw moments from a native history, one row per output time."""

    history = read_selected_moment_history(
        path,
        {label: label for label in labels},
    )
    columns = [history[label] for label in labels]

    return [
        {
            "time": entries[0]["time"],
            **{
                label: entry["value"]
                for label, entry in zip(labels, entries)
            },
        }
        for entries in zip(*columns)
    ]


def compare_native_histories(
    results_dir: Path,
    labels: list[str],
    quantities: dict[str, Callable[[dict[str, float]], float]],
) -> dict[str, list[dict[str, float]]]:
    """Derive comparison quantities from the QMOMKit and Monte Carlo histories.

    Both runs write QMOMKit-native moment histories at the same output times, so
    the compared quantities are derived here rather than written to a separate
    comparison file by the example itself.
    """

    qmomkit_rows = read_moment_rows(
        results_dir / "moment_history.txt",
        labels,
    )
    monte_carlo_rows = read_moment_rows(
        results_dir / "moment_history_mc.txt",
        labels,
    )

    if len(qmomkit_rows) != len(monte_carlo_rows):
        raise ValueError(
            f"{results_dir} QMOMKit and Monte Carlo histories have different output-time counts."
        )

    grouped: dict[str, list[dict[str, float]]] = {
        name: []
        for name in quantities
    }

    for qmomkit_row, monte_carlo_row in zip(qmomkit_rows, monte_carlo_rows):
        if not math.isclose(
            qmomkit_row["time"],
            monte_carlo_row["time"],
            rel_tol=0.0,
            abs_tol=1.0e-14,
        ):
            raise ValueError(
                f"{results_dir} QMOMKit and Monte Carlo histories have different output times."
            )

        for name, evaluate in quantities.items():
            grouped[name].append(
                {
                    "time": qmomkit_row["time"],
                    "qmomkit": evaluate(qmomkit_row),
                    "monte_carlo": evaluate(monte_carlo_row),
                }
            )

    return grouped


def read_selected_moment_history(
    path: Path,
    quantities: dict[str, str],
) -> dict[str, list[dict[str, float]]]:
    with path.open(encoding="utf-8") as handle:
        # QMOMKit centres the "# time" label inside a padded column, so the
        # comment marker is not at the start of the line.
        header = handle.readline().strip().lstrip("#").split()

        if not header or header[0] != "time":
            raise ValueError(
                f"{path} does not contain a moment-history header."
            )

        positions = {
            name: header.index(source_name)
            for name, source_name in quantities.items()
        }

        selected: dict[str, list[dict[str, float]]] = {
            name: []
            for name in quantities
        }

        for line in handle:
            values = line.split()
            if not values:
                continue

            time = float(values[0])

            for name, position in positions.items():
                selected[name].append(
                    {
                        "time": time,
                        "value": float(values[position]),
                    }
                )

    return selected


def read_terminal_quadrature_moments(
    path: Path,
    quantities: dict[str, tuple[float, ...]],
) -> dict[str, list[dict[str, float]]]:
    """Evaluate selected moments from a terminal-quadrature history."""

    with path.open(encoding="utf-8") as handle:
        # QMOMKit centres the "# time" label inside a padded column, so the
        # comment marker is not at the start of the line.
        header = handle.readline().strip().lstrip("#").split()

        if not header or header[0] != "time":
            raise ValueError(
                f"{path} does not contain a quadrature-history header."
            )

        weight_position = header.index("weight")

        dimension_count = len(
            next(iter(quantities.values()))
        )

        coordinate_positions = [
            header.index(f"x{dimension}")
            for dimension in range(dimension_count)
        ]

        totals: dict[float, dict[str, float]] = {}

        for line in handle:
            values = line.split()
            if not values:
                continue

            time = float(values[0])
            weight = float(values[weight_position])

            coordinates = [
                float(values[position])
                for position in coordinate_positions
            ]

            totals.setdefault(
                time,
                {
                    name: 0.0
                    for name in quantities
                },
            )

            for name, orders in quantities.items():
                moment = weight

                for coordinate, order in zip(
                    coordinates,
                    orders,
                ):
                    moment *= coordinate**order

                totals[time][name] += moment

    return {
        name: [
            {
                "time": time,
                "value": totals[time][name],
            }
            for time in sorted(totals)
        ]
        for name in quantities
    }


def combine_moment_histories(
    qmomkit: dict[str, list[dict[str, float]]],
    monte_carlo: dict[str, list[dict[str, float]]],
) -> dict[str, list[dict[str, float]]]:
    if qmomkit.keys() != monte_carlo.keys():
        raise ValueError(
            "QMOMKit and Monte Carlo histories do not contain "
            "the same selected moments."
        )

    combined: dict[str, list[dict[str, float]]] = {}

    for name, qmomkit_rows in qmomkit.items():
        monte_carlo_rows = monte_carlo[name]

        if len(qmomkit_rows) != len(monte_carlo_rows):
            raise ValueError(
                "QMOMKit and Monte Carlo histories have different "
                f"lengths for {name}."
            )

        rows: list[dict[str, float]] = []

        for qmomkit_row, monte_carlo_row in zip(
            qmomkit_rows,
            monte_carlo_rows,
        ):
            if not math.isclose(
                qmomkit_row["time"],
                monte_carlo_row["time"],
                rel_tol=0.0,
                abs_tol=1.0e-14,
            ):
                raise ValueError(
                    "QMOMKit and Monte Carlo histories have different "
                    f"output times for {name}."
                )

            rows.append(
                {
                    "time": qmomkit_row["time"],
                    "qmomkit": qmomkit_row["value"],
                    "monte_carlo": monte_carlo_row["value"],
                }
            )

        combined[name] = rows

    return combined


def spherical_breakup_quadrature_4d_data(
    results_dir: Path,
) -> dict[str, list[dict[str, float]]]:
    """Construct the 4D comparison directly from retained raw histories."""

    # Transported CQMOM moments.
    tracked = read_selected_moment_history(
        results_dir / "moment_history.txt",
        {
            "M_0000": "M[0,0,0,0]",
            "M_1000": "M[1,0,0,0]",
            "M_3000": "M[3,0,0,0]",
        },
    )

    # Untracked moments evaluated from the retained terminal quadrature.
    reconstructed = read_terminal_quadrature_moments(
        results_dir / "quadrature_history.txt",
        {
            "M_3100": (3, 1, 0, 0),
            "M_3010": (3, 0, 1, 0),
            "M_3001": (3, 0, 0, 1),
        },
    )

    # Monte Carlo contains all six selected moments directly.
    monte_carlo = read_selected_moment_history(
        results_dir / "moment_history_mc.txt",
        {
            "M_0000": "M[0,0,0,0]",
            "M_1000": "M[1,0,0,0]",
            "M_3000": "M[3,0,0,0]",
            "M_3100": "M[3,1,0,0]",
            "M_3010": "M[3,0,1,0]",
            "M_3001": "M[3,0,0,1]",
        },
    )

    qmomkit = {
        **tracked,
        **reconstructed,
    }

    return combine_moment_histories(
        qmomkit,
        monte_carlo,
    )


def spherical_breakup_moments_5d_data(
    results_dir: Path,
) -> dict[str, list[dict[str, float]]]:
    """Construct the 5D closure comparison from retained raw histories."""

    tracked = read_selected_moment_history(
        results_dir / "moment_history.txt",
        {
            "M_00000": "M[0,0,0,0,0]",
            "M_10000": "M[1,0,0,0,0]",
            "M_10100": "M[1,0,1,0,0]",
            "M_11000": "M[1,1,0,0,0]",
        },
    )

    generalized = read_terminal_quadrature_moments(
        results_dir / "quadrature_history.txt",
        {
            "M_4over3_0000": (4.0 / 3.0, 0.0, 0.0, 0.0, 0.0),
            "M_4over3_1000": (4.0 / 3.0, 1.0, 0.0, 0.0, 0.0),
        },
    )

    monte_carlo = {
        **read_selected_moment_history(
            results_dir / "moment_history_mc.txt",
            {
                "M_00000": "M[0,0,0,0,0]",
                "M_10000": "M[1,0,0,0,0]",
                "M_10100": "M[1,0,1,0,0]",
                "M_11000": "M[1,1,0,0,0]",
            },
        ),
        **read_selected_moment_history(
            results_dir / "generalized_moment_history_mc.txt",
            {
                "M_4over3_0000": "M[4/3,0,0,0,0]",
                "M_4over3_1000": "M[4/3,1,0,0,0]",
            },
        ),
    }

    return combine_moment_histories(
        {
            **tracked,
            **generalized,
        },
        monte_carlo,
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


def plot_comparison(
    axis: plt.Axes,
    rows: list[dict[str, float]],
) -> None:
    time = [
        row["time"]
        for row in rows
    ]

    qmomkit = [
        row["qmomkit"]
        for row in rows
    ]

    monte_carlo = [
        row["monte_carlo"]
        for row in rows
    ]

    axis.plot(
        time,
        monte_carlo,
        color=STYLE["monte_carlo"][1],
        linewidth=STYLE["monte_carlo"][2],
        zorder=STYLE["monte_carlo"][3],
        label=STYLE["monte_carlo"][0],
    )

    axis.plot(
        time,
        qmomkit,
        color=STYLE["qmomkit"][1],
        linewidth=STYLE["qmomkit"][2],
        zorder=STYLE["qmomkit"][3],
        label=STYLE["qmomkit"][0],
    )


def format_axis(
    axis: plt.Axes,
    title: str | None,
    y_label: str,
) -> None:
    if title is not None:
        axis.set_title(title)

    axis.set_ylabel(y_label)
    axis.grid(True, alpha=0.25)


def save_figure(
    figure: plt.Figure,
    output: Path,
) -> None:
    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure.savefig(
        output,
        bbox_inches="tight",
    )

    plt.close(figure)
    print(output)


def plot_cavitation(
    data: dict[str, list[dict[str, float]]],
    output: Path,
) -> None:
    figure, axes = plt.subplots(
        2,
        3,
        figsize=(7.05, 4.45),
        sharex=True,
    )

    plot_comparison(
        axes[0, 0],
        data["M_00"],
    )
    format_axis(
        axes[0, 0],
        "Number density",
        r"$M_{0,0}$",
    )

    plot_comparison(
        axes[0, 1],
        data["M_10/M_00"],
    )
    format_axis(
        axes[0, 1],
        "Mean bubble radius",
        r"$M_{1,0}$",
    )

    plot_comparison(
        axes[0, 2],
        data["M_01/M_00"],
    )
    format_axis(
        axes[0, 2],
        "Mean radial velocity",
        r"$M_{0,1}$",
    )

    plot_comparison(
        axes[1, 0],
        data["M_11/M_00"],
    )
    format_axis(
        axes[1, 0],
        "Radius-velocity moment",
        r"$M_{1,1}$",
    )

    plot_comparison(
        axes[1, 1],
        data["M_20/M_00"],
    )
    format_axis(
        axes[1, 1],
        "Second radius moment",
        r"$M_{2,0}$",
    )

    plot_comparison(
        axes[1, 2],
        data["M_02/M_00"],
    )
    format_axis(
        axes[1, 2],
        "Second velocity moment",
        r"$M_{0,2}$",
    )

    for axis in axes[1, :]:
        axis.set_xlabel("Time")

    handles, labels = axes[0, 0].get_legend_handles_labels()

    figure.legend(
        handles,
        labels,
        loc="upper center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 0.91),
    )

    figure.subplots_adjust(
        left=0.10,
        right=0.985,
        bottom=0.12,
        top=0.8,
        hspace=0.28,
        wspace=0.43,
    )

    save_figure(
        figure,
        output,
    )


def plot_condensing_droplets(
    data: dict[str, list[dict[str, float]]],
    output: Path,
) -> None:
    figure, axes = plt.subplots(
        2,
        3,
        figsize=(7.05, 3.7),
        sharex=True,
    )

    panels = (
        (
            axes[0, 0],
            data["M_000"],
            "Total number density",
            r"$M_{0,0,0}$ [-]",
        ),
        (
            axes[0, 1],
            data["M_300/M_200"],
            "Characteristic radius",
            r"$M_{3,0,0}/M_{2,0,0}$ [$\mu$m]",
        ),
        (
            axes[0, 2],
            data["M_010/M_000"],
            "Mean vertical velocity",
            r"$M_{0,1,0}/M_{0,0,0}$ [m/s]",
        ),
        (
            axes[1, 0],
            data["M_001/M_000"],
            "Mean temperature",
            r"$M_{0,0,1}/M_{0,0,0}$ [K]",
        ),
        (
            axes[1, 2],
            data["sigma_r^2"],
            "Radius variance",
            r"$\sigma_r^2$ [$\mu$m$^2$]",
        ),
        (
            axes[1, 1],
            data["M_101/M_000"],
            "Radius-weighted mean temperature",
            r"$M_{1,0,1}/M_{0,0,0}$ [$\mu$m K]",
        ),
    )

    for axis, rows, title, y_label in panels:
        axis.axvspan(
            0.20,
            0.40,
            color="#E69F00",
            alpha=0.12,
            zorder=0,
        )

        plot_comparison(
            axis,
            rows,
        )

        format_axis(
            axis,
            title,
            y_label,
        )

    for axis in axes[1, :]:
        axis.set_xlabel("Time [s]")

    handles, labels = axes[0, 0].get_legend_handles_labels()

    figure.legend(
        handles,
        labels,
        loc="upper center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 1),
    )

    figure.subplots_adjust(
        left=0.11,
        right=0.985,
        bottom=0.09,
        top=0.88,
        hspace=0.25,
        wspace=0.45,
    )

    save_figure(
        figure,
        output,
    )


# Phase-space order is (r, w, T); these are the raw moments the compared
# droplet quantities are built from.
DROPLET_MOMENTS = [
    "M[0,0,0]",
    "M[3,0,0]",
    "M[2,0,0]",
    "M[0,1,0]",
    "M[0,0,1]",
    "M[1,0,0]",
    "M[1,0,1]",
]


def condensing_droplet_data(
    results_dir: Path,
) -> dict[str, list[dict[str, float]]]:
    """Present SI-radius droplet data in micrometres for readability."""

    def mean_radius(row: dict[str, float]) -> float:
        return row["M[1,0,0]"] / row["M[0,0,0]"]

    def mean_temperature(row: dict[str, float]) -> float:
        return row["M[0,0,1]"] / row["M[0,0,0]"]

    # Radii are stored in metres; the micrometre factors are presentation only.
    return compare_native_histories(
        results_dir,
        DROPLET_MOMENTS,
        {
            "M_000": lambda row: row["M[0,0,0]"],
            "M_300/M_200": lambda row: 1.0e6 * row["M[3,0,0]"] / row["M[2,0,0]"],
            "M_010/M_000": lambda row: row["M[0,1,0]"] / row["M[0,0,0]"],
            "M_001/M_000": mean_temperature,
            "sigma_r^2": lambda row: 1.0e12 * (row["M[2,0,0]"] / row["M[0,0,0]"] - mean_radius(row) ** 2),
            "M_101/M_000": lambda row: 1.0e6 * row["M[1,0,1]"] / row["M[0,0,0]"],
        },
    )


# Phase-space order is (R, U); the cavitation figure compares the zeroth raw
# moment and the population-normalised first, mixed, and second moments.
CAVITATION_MOMENTS = [
    "M[0,0]",
    "M[1,0]",
    "M[0,1]",
    "M[1,1]",
    "M[2,0]",
    "M[0,2]",
]


def cavitation_data(
    results_dir: Path,
) -> dict[str, list[dict[str, float]]]:
    """Population-normalised cavitation moments compared against Monte Carlo."""

    return compare_native_histories(
        results_dir,
        CAVITATION_MOMENTS,
        {
            "M_00": lambda row: row["M[0,0]"],
            "M_10/M_00": lambda row: row["M[1,0]"] / row["M[0,0]"],
            "M_01/M_00": lambda row: row["M[0,1]"] / row["M[0,0]"],
            "M_11/M_00": lambda row: row["M[1,1]"] / row["M[0,0]"],
            "M_20/M_00": lambda row: row["M[2,0]"] / row["M[0,0]"],
            "M_02/M_00": lambda row: row["M[0,2]"] / row["M[0,0]"],
        },
    )


def plot_non_spherical_particles(
    data: dict[str, list[dict[str, float]]],
    output: Path,
) -> None:
    figure, axes = plt.subplots(
        3,
        3,
        figsize=(7.05, 5.35),
        sharex=True,
    )

    plot_comparison(
        axes[0, 0],
        data["M_00000"],
    )
    format_axis(
        axes[0, 0],
        "Number density",
        r"$M_{0,0,0,0,0}$ [-]",
    )

    plot_comparison(
        axes[0, 1],
        data["M_30000/M_20000"],
    )
    format_axis(
        axes[0, 1],
        "Sauter mean diameter",
        r"$d_{32}$ [m]",
    )

    plot_comparison(
        axes[0, 2],
        data["M_01000/M_00000"],
    )
    format_axis(
        axes[0, 2],
        "Mean aspect ratio",
        r"$<\chi>$ [-]",
    )

    plot_comparison(
        axes[1, 0],
        data["M_30100/M_00000"],
    )
    format_axis(
        axes[1, 0],
        "Size-weighted $u$ moment",
        r"$q_x$ [m$^4$/s]",
    )

    plot_comparison(
        axes[1, 1],
        data["M_30010/M_00000"],
    )
    format_axis(
        axes[1, 1],
        "Size-weighted $v$ moment",
        r"$q_y$ [m$^4$/s]",
    )

    plot_comparison(
        axes[1, 2],
        data["M_30001/M_00000"],
    )
    format_axis(
        axes[1, 2],
        "Size-weighted $w$ moment",
        r"$q_z$ [m$^4$/s]",
    )

    plot_comparison(
        axes[2, 0],
        data["Theta_x"],
    )
    format_axis(
        axes[2, 0],
        "Velocity variance in $x$",
        r"$\Theta_x$ [m$^2$/s$^{2}$]",
    )

    plot_comparison(
        axes[2, 1],
        data["Theta_y"],
    )
    format_axis(
        axes[2, 1],
        "Velocity variance in $y$",
        r"$\Theta_y$ [m$^2$/s$^{2}$]",
    )

    plot_comparison(
        axes[2, 2],
        data["Theta_z"],
    )
    format_axis(
        axes[2, 2],
        "Velocity variance in $z$",
        r"$\Theta_z$ [m$^2$/s$^{2}$]",
    )

    for axis in axes[2, :]:
        axis.set_xlabel("Time [s]")

    handles, labels = axes[0, 0].get_legend_handles_labels()

    figure.legend(
        handles,
        labels,
        loc="upper center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 0.95),
    )

    figure.subplots_adjust(
        left=0.09,
        right=0.985,
        bottom=0.10,
        top=0.86,
        hspace=0.475,
        wspace=0.5,
    )

    save_figure(
        figure,
        output,
    )


def plot_spherical_breakup_quadrature_4d(
    data: dict[str, list[dict[str, float]]],
    output: Path,
) -> None:
    figure, axes = plt.subplots(
        2,
        3,
        figsize=(7.05, 4.3),
        sharex=True,
    )

    def ratio_rows(
        numerator: str,
        denominator: str,
        factor: float = 1.0,
    ) -> list[dict[str, float]]:
        return [
            {
                "time": numerator_row["time"],
                "qmomkit": (
                    factor
                    * numerator_row["qmomkit"]
                    / denominator_row["qmomkit"]
                ),
                "monte_carlo": (
                    factor
                    * numerator_row["monte_carlo"]
                    / denominator_row["monte_carlo"]
                ),
            }
            for numerator_row, denominator_row in zip(
                data[numerator],
                data[denominator],
            )
        ]

    panels = (
        (
            axes[0, 0],
            data["M_0000"],
            r"$M_{0,0,0,0}$ [m$^{-3}$]",
        ),
        (
            axes[0, 1],
            ratio_rows(
                "M_1000",
                "M_0000",
                1.0e3,
            ),
            r"$M_{1,0,0,0}/M_{0,0,0,0}$ [mm]",
        ),
        (
            axes[0, 2],
            data["M_3000"],
            r"$M_{3,0,0,0}$ [-]",
        ),
        (
            axes[1, 0],
            data["M_3100"],
            r"$M_{3,1,0,0}$ [m/s]",
        ),
        (
            axes[1, 1],
            data["M_3010"],
            r"$M_{3,0,1,0}$ [m/s]",
        ),
        (
            axes[1, 2],
            data["M_3001"],
            r"$M_{3,0,0,1}$ [m/s]",
        ),
    )

    for axis, rows, y_label in panels:
        plot_comparison(
            axis,
            rows,
        )

        format_axis(
            axis,
            None,
            y_label,
        )

    reference_volume_moment = data["M_3000"][0]["qmomkit"]

    axes[0, 2].set_ylim(
        0.9 * reference_volume_moment,
        1.1 * reference_volume_moment,
    )

    for axis in axes[1, :]:
        axis.set_xlabel("Time [s]")

    handles, labels = axes[0, 0].get_legend_handles_labels()

    figure.legend(
        handles,
        labels,
        loc="upper center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 0.88),
    )

    figure.subplots_adjust(
        left=0.10,
        right=0.985,
        bottom=0.12,
        top=0.8,
        hspace=0.2,
        wspace=0.43,
    )

    save_figure(
        figure,
        output,
    )


def plot_spherical_breakup_moments_5d(
    data: dict[str, list[dict[str, float]]],
    output: Path,
) -> None:
    figure, axes = plt.subplots(
        2,
        3,
        figsize=(7.05, 4.3),
        sharex=True,
    )

    panels = (
        (
            axes[0, 0],
            data["M_00000"],
            r"$M_{0,0,0,0,0}$ [m$^{-3}$]",
        ),
        (
            axes[0, 1],
            data["M_10000"],
            r"$M_{1,0,0,0,0}$ [-]",
        ),
        (
            axes[0, 2],
            data["M_10100"],
            r"$M_{1,0,1,0,0}$ [m/s]",
        ),
        (
            axes[1, 0],
            data["M_11000"],
            r"$M_{1,1,0,0,0}$ [K]",
        ),
        (
            axes[1, 1],
            data["M_4over3_0000"],
            r"$M_{4/3,0,0,0,0}$ [m]",
        ),
        (
            axes[1, 2],
            data["M_4over3_1000"],
            r"$M_{4/3,1,0,0,0}$ [m K]",
        ),
    )

    for axis, rows, y_label in panels:
        plot_comparison(
            axis,
            rows,
        )

        format_axis(
            axis,
            None,
            y_label,
        )

    reference_volume_moment = data["M_10000"][0]["qmomkit"]

    axes[0, 1].set_ylim(
        0.9 * reference_volume_moment,
        1.1 * reference_volume_moment,
    )

    for axis in axes[1, :]:
        axis.set_xlabel("Time [s]")

    handles, labels = axes[0, 0].get_legend_handles_labels()

    figure.legend(
        handles,
        labels,
        loc="upper center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 0.9),
    )

    figure.subplots_adjust(
        left=0.10,
        right=0.985,
        bottom=0.12,
        top=0.8,
        hspace=0.2,
        wspace=0.43,
    )

    save_figure(
        figure,
        output,
    )


def plot_aerosol_evaporation(
    results_dir: Path,
    output: Path,
) -> None:
    def read_history(
        filename: str,
    ) -> list[tuple[float, float, float, float]]:
        rows: list[
            tuple[float, float, float, float]
        ] = []

        with (
            results_dir / filename
        ).open(encoding="utf-8") as handle:
            for line in handle:
                if line and not line.startswith("#"):
                    (
                        time,
                        number_density,
                        d32,
                        variance,
                    ) = (
                        float(value)
                        for value in line.split()
                    )

                    rows.append(
                        (
                            time,
                            number_density,
                            d32,
                            variance,
                        )
                    )

        return rows

    mc_history = read_history(
        "moment_history_mc.txt"
    )

    qmom_histories = {
        2: read_history("moment_history_qmom_n2.txt"),
        4: read_history("moment_history_qmom_n4.txt"),
        6: read_history("moment_history_qmom_n6.txt"),
        8: read_history("moment_history_qmom_n8.txt"),
        10: read_history("moment_history_qmom_n10.txt"),
    }

    initial_diameters: list[float] = []

    with (
        results_dir / "initial_mc_samples.txt"
    ).open(encoding="utf-8") as handle:
        for line in handle:
            if line and not line.startswith("#"):
                initial_diameters.append(
                    float(line)
                )

    figure, axes = plt.subplots(
        2,
        2,
        figsize=(6.0, 4.5),
    )

    histogram_axis = axes[0, 0]

    log_minimum = math.log10(0.8)
    log_bin_width = (
        math.log10(1000.0) - log_minimum
    ) / 60.0

    histogram_bins = [
        10.0
        ** (
            log_minimum
            + bin_index * log_bin_width
        )
        for bin_index in range(61)
    ]

    bin_counts = [0] * 60

    for diameter in initial_diameters:
        bin_index = min(
            59,
            max(
                0,
                int(
                    (
                        math.log10(diameter)
                        - log_minimum
                    )
                    / log_bin_width
                ),
            ),
        )

        bin_counts[bin_index] += 1

    number_fraction_per_decade = [
        count
        / (
            len(initial_diameters)
            * log_bin_width
        )
        for count in bin_counts
    ]

    if not math.isclose(
        sum(
            fraction * log_bin_width
            for fraction
            in number_fraction_per_decade
        ),
        1.0,
        rel_tol=1.0e-12,
        abs_tol=1.0e-12,
    ):
        raise ValueError(
            "The retained-population histogram "
            "is not normalized."
        )

    nonzero_bins = [
        (
            math.sqrt(
                histogram_bins[index]
                * histogram_bins[index + 1]
            ),
            number_fraction_per_decade[index],
            histogram_bins[index + 1]
            - histogram_bins[index],
        )
        for index in range(60)
        if bin_counts[index] > 0
    ]

    histogram_axis.bar(
        [
            bin_center
            for bin_center, _, _
            in nonzero_bins
        ],
        [
            fraction
            for _, fraction, _
            in nonzero_bins
        ],
        width=[
            bin_width
            for _, _, bin_width
            in nonzero_bins
        ],
        align="center",
        color="#D9D9D9",
        edgecolor="#8C8C8C",
        linewidth=0.35,
        zorder=1,
    )

    histogram_axis.set_xscale("log")
    histogram_axis.set_yscale("log")
    histogram_axis.set_xlim(
        0.8,
        1000.0,
    )

    histogram_axis.set_ylim(
        bottom=0.7
        * min(
            fraction
            for _, fraction, _
            in nonzero_bins
        )
    )

    histogram_axis.set_xlabel(
        r"$d$ [$\mu$m]"
    )

    histogram_axis.set_ylabel(
        "Normalized number distribution"
    )

    histogram_axis.grid(
        True,
        which="both",
        alpha=0.25,
    )

    qmom_styles = {
        2: ("#4D4D4D", ":", 1.20),
        4: ("#6A6A6A", "--", 1.20),
        6: ("#0072B2", "-.", 1.20),
        8: (
            "#005A9C",
            (0, (4, 1, 1, 1)),
            1.20,
        ),
        10: ("#003B6F", "-", 1.35),
    }

    panels = (
        (
            axes[0, 1],
            1,
            r"$M_0$ [$10^5$ m$^{-3}$]",
            1.0e-5,
        ),
        (
            axes[1, 0],
            2,
            r"$d_{32}$ [$\mu$m]",
            1.0,
        ),
        (
            axes[1, 1],
            3,
            r"$\sigma_d^2$ [$\mu$m$^2$]",
            1.0,
        ),
    )

    for axis, column, y_label, scale in panels:
        axis.plot(
            [
                row[0]
                for row in mc_history
            ],
            [
                scale * row[column]
                for row in mc_history
            ],
            color="#B3B3B3",
            linewidth=2.6,
            zorder=1,
            label="Monte Carlo",
        )

        for node_count, history in qmom_histories.items():
            (
                color,
                line_style,
                line_width,
            ) = qmom_styles[node_count]

            axis.plot(
                [
                    row[0]
                    for row in history
                ],
                [
                    scale * row[column]
                    for row in history
                ],
                color=color,
                linestyle=line_style,
                linewidth=line_width,
                zorder=2,
                label=(
                    f"QMOMKit, $N={node_count}$"
                ),
            )

        axis.set_xlim(
            0.0,
            20.0,
        )

        axis.set_ylabel(
            y_label
        )

        if column == 1:
            axis.set_ylim(
                2.20,
                2.26,
            )

            axis.ticklabel_format(
                axis="y",
                style="plain",
                useOffset=False,
            )

        axis.grid(
            True,
            alpha=0.25,
        )

    axes[0, 1].set_xlabel("Time [s]")
    axes[1, 0].set_xlabel("Time [s]")
    axes[1, 1].set_xlabel("Time [s]")

    history_handles, history_labels = (
        axes[0, 1].get_legend_handles_labels()
    )

    figure.legend(
        history_handles,
        history_labels,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.95),
    )

    figure.subplots_adjust(
        left=0.12,
        right=0.985,
        bottom=0.13,
        top=0.84,
        wspace=0.4,
        hspace=0.34,
    )

    save_figure(
        figure,
        output,
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create paper application-example figures from QMOMKit output.")
    parser.add_argument(
        "--examples-build-dir",
        type=Path,
        default=REPOSITORY_ROOT / "paper" / "examples" / "build",
        help="Paper examples build directory containing the per-case output directories.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SCRIPT_DIRECTORY,
        help="Directory for vector-PDF figures.",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    examples_build_directory = arguments.examples_build_dir
    figure_directory = arguments.output_dir
    set_paper_style()

    plot_aerosol_evaporation(
        examples_build_directory / "1d_aerosols_evaporation" / "output",
        figure_directory / "aerosol_evaporation_moments.pdf",
    )

    plot_cavitation(
        cavitation_data(
            examples_build_directory / "2d_bubbles_cavitation" / "output",
        ),
        figure_directory / "cavitation_moments.pdf",
    )

    plot_condensing_droplets(
        condensing_droplet_data(
            examples_build_directory / "3d_droplets_condensation-nucleation" / "output",
        ),
        figure_directory / "condensing_droplets_moments.pdf",
    )

    breakup_quadrature_4d_dir = examples_build_directory / "4d_particles_breakup_recoil" / "output"

    plot_spherical_breakup_quadrature_4d(
        spherical_breakup_quadrature_4d_data(
            breakup_quadrature_4d_dir
        ),
        figure_directory
        / "cqmom_spherical_breakup_quadrature_4d_moments.pdf",
    )

    breakup_moments_5d_dir = examples_build_directory / "5d_particles_breakup_heat-transfer" / "output"

    breakup_moments_5d = spherical_breakup_moments_5d_data(
        breakup_moments_5d_dir
    )

    plot_spherical_breakup_moments_5d(
        breakup_moments_5d,
        figure_directory
        / "cqmom_spherical_breakup_moments_5d_moments.pdf",
    )


if __name__ == "__main__":
    main()
