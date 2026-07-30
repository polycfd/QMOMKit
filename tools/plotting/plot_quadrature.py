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


def is_qmom_labels(labels: list[str]) -> bool:
    return labels == ["time", "node_index", "weight", "x0"]


def cqmom_columns(labels: list[str]) -> tuple[list[str], list[str]]:
    if len(labels) < 6 or labels[:2] != ["time", "node_index"]:
        raise ValueError("quadrature file is neither standard QMOM nor CQMOM output")

    cursor = 2
    path_labels: list[str] = []
    while cursor < len(labels) and re.fullmatch(r"path_\d+", labels[cursor]):
        path_labels.append(labels[cursor])
        cursor += 1

    if not path_labels or cursor >= len(labels) or labels[cursor] != "weight":
        raise ValueError("quadrature file is neither standard QMOM nor CQMOM output")
    cursor += 1

    x_labels = labels[cursor:]
    if not x_labels or not all(re.fullmatch(r"x\d+", label) for label in x_labels):
        raise ValueError("quadrature file is neither standard QMOM nor CQMOM output")
    return path_labels, x_labels


def rows_by_node(data: np.ndarray) -> list[tuple[int, np.ndarray]]:
    node_column = 1
    nodes = sorted({int(node) for node in data[:, node_column]})
    grouped: list[tuple[int, np.ndarray]] = []
    for node in nodes:
        rows = data[data[:, node_column] == float(node)]
        rows = rows[np.argsort(rows[:, 0], kind="stable")]
        grouped.append((node, rows))
    return grouped


def plot_node_lines(ax, data: np.ndarray, value_column: int, ylabel: str) -> None:
    for node, rows in rows_by_node(data):
        ax.plot(rows[:, 0], rows[:, value_column], marker="o", markersize=3, label=f"node {node}")
    ax.set_xlabel("time")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small")


def save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_qmom(data: np.ndarray, output: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    plot_node_lines(axes[0], data, 2, "weight")
    plot_node_lines(axes[1], data, 3, "x0")
    fig.suptitle("QMOM quadrature history")
    save(fig, output)


def plot_cqmom_history(labels: list[str], data: np.ndarray, x_labels: list[str], path: Path) -> None:
    shown_x = x_labels[:2]
    fig, axes = plt.subplots(1 + len(shown_x), 1, figsize=(8, 3 + 2.5 * len(shown_x)), sharex=True)
    axes = np.atleast_1d(axes)
    plot_node_lines(axes[0], data, labels.index("weight"), "weight")
    axes[0].set_title("CQMOM quadrature history")
    for axis, x_label in zip(axes[1:], shown_x):
        plot_node_lines(axis, data, labels.index(x_label), x_label)
    save(fig, path)


def plot_cqmom_final(labels: list[str], data: np.ndarray, x_labels: list[str], path: Path) -> None:
    final_time = np.max(data[:, 0])
    final_rows = data[data[:, 0] == final_time]
    final_rows = final_rows[np.argsort(final_rows[:, 1], kind="stable")]
    weights = final_rows[:, labels.index("weight")]
    sizes = 80 + 220 * weights / max(float(np.max(weights)), 1.0e-30)
    fig, ax = plt.subplots(figsize=(6, 5))

    if len(x_labels) == 1:
        xvals = final_rows[:, labels.index(x_labels[0])]
        ax.scatter(xvals, weights, s=sizes, c=weights, cmap="viridis", edgecolor="black")
        ax.set_xlabel(x_labels[0])
        ax.set_ylabel("weight")
    else:
        first, second = x_labels[:2]
        xvals = final_rows[:, labels.index(first)]
        yvals = final_rows[:, labels.index(second)]
        ax.scatter(xvals, yvals, s=sizes, c=weights, cmap="viridis", edgecolor="black")
        ax.set_xlabel(first)
        ax.set_ylabel(second)
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
    if is_qmom_labels(labels):
        output = output.with_suffix(".png")
        plot_qmom(data, output)
        print(output)
    else:
        _path_labels, x_labels = cqmom_columns(labels)
        history = output.with_name(output.name + "_history.png")
        final = output.with_name(output.name + "_final.png")
        plot_cqmom_history(labels, data, x_labels, history)
        plot_cqmom_final(labels, data, x_labels, final)
        print(history)
        print(final)


if __name__ == "__main__":
    main()
