#!/usr/bin/env python3
"""Plot raw moment histories from QMOMKit MomentWriter output."""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "qmomkit-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "qmomkit-cache"))

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_qmomkit_table(path: Path) -> tuple[list[str], np.ndarray]:
    """Read a whitespace-delimited QMOMKit table with a '# time ...' header."""

    labels = None
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("#"):
                candidate = stripped[1:].split()
                if candidate and candidate[0] == "time":
                    labels = candidate
                    break

    if labels is None:
        raise ValueError(f"{path} does not contain a '# time ...' header")

    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] != len(labels):
        raise ValueError(f"{path} has {data.shape[1]} data columns but {len(labels)} header labels")
    return labels, data


def choose_moments(labels: list[str], requested: list[str] | None) -> list[str]:
    """Return requested moment labels, or all non-time columns."""

    moments = labels[1:]
    if not requested:
        return moments

    missing = [name for name in requested if name not in moments]
    if missing:
        raise ValueError("missing moment columns: " + ", ".join(missing))
    return requested


def plot_moments(path: Path, output: Path, requested: list[str] | None) -> None:
    """Plot selected moment columns against time."""

    output.parent.mkdir(parents=True, exist_ok=True)
    labels, data = read_qmomkit_table(path)
    moments = choose_moments(labels, requested)
    time = data[:, 0]

    fig, ax = plt.subplots(figsize=(8, 5))
    for moment in moments:
        ax.plot(time, data[:, labels.index(moment)], marker="o", markersize=3, label=moment)

    ax.set_title(path.name)
    ax.set_xlabel("time")
    ax.set_ylabel("raw moment")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small")
    fig.tight_layout()
    fig.savefig(output, dpi=160)


def main() -> None:
    """Command-line entry point."""

    parser = argparse.ArgumentParser(description="Plot QMOMKit MomentWriter output.")
    parser.add_argument("moment_file", type=Path, help="whitespace-delimited MomentWriter text file")
    parser.add_argument("--moments", nargs="+", help="moment labels to plot, for example 'M[0]' 'M[2]'")
    parser.add_argument("--output", type=Path, help="figure path, default: <moment_file stem>_moments.png")
    args = parser.parse_args()

    output = args.output or args.moment_file.with_name(args.moment_file.stem + "_moments.png")
    plot_moments(args.moment_file, output, args.moments)
    print(output)


if __name__ == "__main__":
    main()
