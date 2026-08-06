#!/usr/bin/env python3
"""Plot QMOMKit paper Runge--Kutta and splitting convergence data."""

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


def reference_line(axis: plt.Axes, x: np.ndarray, y: np.ndarray, order: int, label: str) -> None:
    anchor = y[0]
    axis.loglog(x, anchor * (x / x[0]) ** order, "--", color="0.45", label=label)


def plot_rk(path: Path, output: Path) -> None:
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    data = np.atleast_1d(data)
    figure, axis = plt.subplots(figsize=(7.2, 5.0))
    orders = {"ssp_rk3": 3, "rk23": 3, "rk45": 5}
    for scheme, order in orders.items():
        rows = data[data["scheme"] == scheme]
        axis.loglog(rows["step_size"], rows["absolute_error"], "o-", label=scheme)
        reference_line(axis, rows["step_size"], rows["absolute_error"], order, f"order {order}")
    axis.set_xlabel("fixed step size")
    axis.set_ylabel("final absolute error")
    axis.grid(True, which="both", alpha=0.3)
    axis.legend(fontsize="small", ncol=2)
    figure.tight_layout()
    figure.savefig(output)
    plt.close(figure)


def plot_splitting(path: Path, output: Path) -> None:
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    data = np.atleast_1d(data)
    figure, axis = plt.subplots(figsize=(7.2, 5.0))
    orders = {"advsrc": 1, "srcadv": 1, "strang": 2}
    for mode, order in orders.items():
        rows = data[data["split_mode"] == mode]
        axis.loglog(rows["outer_step_size"], rows["max_scaled_moment_error"], "o-", label=mode)
        reference_line(axis, rows["outer_step_size"], rows["max_scaled_moment_error"], order, f"order {order}")
    axis.set_xlabel("outer splitting step size")
    axis.set_ylabel("maximum scaled moment error")
    axis.grid(True, which="both", alpha=0.3)
    axis.legend(fontsize="small", ncol=2)
    figure.tight_layout()
    figure.savefig(output)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot RK and operator-splitting convergence studies.")
    parser.add_argument("--rk-input", type=Path, default=Path("rk_convergence.csv"))
    parser.add_argument("--splitting-input", type=Path, default=Path("splitting_convergence.csv"))
    parser.add_argument("--rk-output", type=Path, default=Path("rk_convergence.pdf"))
    parser.add_argument("--splitting-output", type=Path, default=Path("splitting_convergence.pdf"))
    args = parser.parse_args()

    args.rk_output.parent.mkdir(parents=True, exist_ok=True)
    args.splitting_output.parent.mkdir(parents=True, exist_ok=True)
    plot_rk(args.rk_input, args.rk_output)
    plot_splitting(args.splitting_input, args.splitting_output)
    print(args.rk_output)
    print(args.splitting_output)


if __name__ == "__main__":
    main()
