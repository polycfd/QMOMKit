#!/usr/bin/env python3
"""Plot the adaptive CQMOM branch-count verification history."""

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


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot adaptive CQMOM branch counts and velocity variance.")
    parser.add_argument("csv", nargs="?", type=Path, default=Path("adaptive_branch_quadrature.csv"))
    parser.add_argument("--output", type=Path, default=Path("adaptive_branch_quadrature.pdf"))
    parser.add_argument("--pulse-start", type=float, default=4.0)
    parser.add_argument("--pulse-end", type=float, default=4.05)
    args = parser.parse_args()

    data = np.genfromtxt(args.csv, delimiter=",", names=True)
    if data.size == 0:
        raise ValueError(f"{args.csv} does not contain data rows")

    data = np.atleast_1d(data)
    figure, axes = plt.subplots(2, 1, figsize=(7.2, 6.0), sharex=True)
    axes[0].step(data["time"], data["branch_0_active_nodes"], where="post", label="size branch 0")
    axes[0].step(data["time"], data["branch_1_active_nodes"], where="post", label="size branch 1")
    axes[0].set_ylabel("active velocity nodes")
    axes[0].set_ylim(0.7, 3.3)
    axes[0].set_yticks([1, 2, 3])
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].plot(data["time"], data["velocity_variance"], color="black", label="velocity variance")
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("velocity variance")
    axes[1].grid(True, alpha=0.3)

    for axis in axes:
        axis.axvspan(args.pulse_start, args.pulse_end, color="tab:orange", alpha=0.2, label="reinjection" if axis is axes[0] else None)

    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output)
    plt.close(figure)
    print(args.output)


if __name__ == "__main__":
    main()
