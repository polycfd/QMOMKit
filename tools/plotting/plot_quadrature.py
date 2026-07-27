#!/usr/bin/env python3
"""Plot basic quadrature histories from QMOMKit QuadratureWriter output."""

from __future__ import annotations

import argparse
import os
import re
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
            fields = line.strip()[1:].split() if line.lstrip().startswith("#") else []
            if fields and fields[0] == "time":
                labels = fields
                break
    if labels is None:
        raise ValueError(f"{path} does not contain a '# time ...' header")

    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] != len(labels):
        raise ValueError(f"{path} has {data.shape[1]} data columns but {len(labels)} header labels")
    return labels, data


def qmom_columns(labels: list[str]) -> tuple[list[str], list[str]] | None:
    columns = labels[1:]
    weights = columns[0::2]
    abscissas = columns[1::2]
    if len(columns) % 2 == 0 and all(w == f"w{i}" and x == f"x{i}" for i, (w, x) in enumerate(zip(weights, abscissas))):
        return weights, abscissas
    return None


def cqmom_columns(labels: list[str]) -> tuple[list[int], list[int], dict[tuple[int, int], str]]:
    nodes: list[int] = []
    dims: set[int] = set()
    xs: dict[tuple[int, int], str] = {}
    for label in labels[1:]:
        if match := re.fullmatch(r"w(\d+)", label):
            nodes.append(int(match.group(1)))
        elif match := re.fullmatch(r"x(\d+)_(\d+)", label):
            dim, node = int(match.group(1)), int(match.group(2))
            dims.add(dim)
            xs[(node, dim)] = label
    if not nodes or not dims:
        raise ValueError("quadrature file is neither standard QMOM nor CQMOM output")
    return nodes, sorted(dims), xs


def plot_lines(ax, labels: list[str], data: np.ndarray, names: list[str], ylabel: str) -> None:
    for name in names:
        ax.plot(data[:, 0], data[:, labels.index(name)], marker="o", markersize=3, label=name)
    ax.set_xlabel("time")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small")


def save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_qmom(labels: list[str], data: np.ndarray, output: Path) -> None:
    weights, abscissas = qmom_columns(labels) or ([], [])
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    plot_lines(axes[0], labels, data, weights, "weight")
    plot_lines(axes[1], labels, data, abscissas, "abscissa")
    fig.suptitle("QMOM quadrature history")
    save(fig, output)


def plot_cqmom_history(labels: list[str], data: np.ndarray, nodes: list[int], dims: list[int], xs: dict[tuple[int, int], str], path: Path) -> None:
    shown_dims = dims[:2]
    fig, axes = plt.subplots(1 + len(shown_dims), 1, figsize=(8, 3 + 2.5 * len(shown_dims)), sharex=True)
    axes = np.atleast_1d(axes)
    plot_lines(axes[0], labels, data, [f"w{node}" for node in nodes], "weight")
    axes[0].set_title("CQMOM quadrature history")
    for axis, dim in zip(axes[1:], shown_dims):
        plot_lines(axis, labels, data, [xs[(node, dim)] for node in nodes], f"x{dim}")
    save(fig, path)


def plot_cqmom_final(labels: list[str], data: np.ndarray, nodes: list[int], dims: list[int], xs: dict[tuple[int, int], str], path: Path) -> None:
    final = data[-1]
    weights = np.array([final[labels.index(f"w{node}")] for node in nodes])
    sizes = 80 + 220 * weights / max(float(np.max(weights)), 1.0e-30)
    fig, ax = plt.subplots(figsize=(6, 5))

    if len(dims) == 1:
        xvals = [final[labels.index(xs[(node, dims[0])])] for node in nodes]
        ax.scatter(xvals, weights, s=sizes, c=weights, cmap="viridis", edgecolor="black")
        ax.set_xlabel(f"x{dims[0]}")
        ax.set_ylabel("weight")
    else:
        first, second = dims[:2]
        xvals = [final[labels.index(xs[(node, first)])] for node in nodes]
        yvals = [final[labels.index(xs[(node, second)])] for node in nodes]
        ax.scatter(xvals, yvals, s=sizes, c=weights, cmap="viridis", edgecolor="black")
        ax.set_xlabel(f"x{first}")
        ax.set_ylabel(f"x{second}")
    ax.set_title("CQMOM final terminal nodes")
    ax.grid(True, alpha=0.3)
    save(fig, path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot QMOMKit QuadratureWriter output.")
    parser.add_argument("quadrature_file", type=Path, help="whitespace-delimited QuadratureWriter text file")
    parser.add_argument("--output", type=Path, help="figure path or prefix, default: input file stem")
    args = parser.parse_args()

    labels, data = read_qmomkit_table(args.quadrature_file)
    output = args.output or args.quadrature_file.with_suffix("")
    if qmom_columns(labels):
        output = output.with_suffix(".png")
        plot_qmom(labels, data, output)
        print(output)
    else:
        nodes, dims, xs = cqmom_columns(labels)
        history = output.with_name(output.name + "_history.png")
        final = output.with_name(output.name + "_final.png")
        plot_cqmom_history(labels, data, nodes, dims, xs, history)
        plot_cqmom_final(labels, data, nodes, dims, xs, final)
        print(history)
        print(final)


if __name__ == "__main__":
    main()
