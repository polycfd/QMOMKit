#!/usr/bin/env python3
"""Plot the adaptive CQMOM paper case."""

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


data = np.atleast_1d(np.genfromtxt("adaptive_branch_quadrature.csv", delimiter=",", names=True))
if data.size == 0:
    raise ValueError("adaptive_branch_quadrature.csv does not contain data rows")

figure, axes = plt.subplots(2, 1, figsize=(7.2, 6.0), sharex=True)
branch_indices = np.unique(data["branch_index"])
for branch_index in branch_indices:
    branch_rows = data[data["branch_index"] == branch_index]
    axes[0].step(branch_rows["time"], branch_rows["active_nodes"], where="post", label=f"size branch {int(branch_index)}")
axes[0].set_ylabel("active velocity nodes")
axes[0].grid(True, alpha=0.3)
axes[0].legend()

global_rows = data[data["branch_index"] == branch_indices[0]]
axes[1].plot(global_rows["time"], global_rows["velocity_variance"], color="black")
axes[1].set_xlabel("time")
axes[1].set_ylabel("velocity variance")
axes[1].grid(True, alpha=0.3)

for axis in axes:
    axis.axvspan(4.0, 4.05, color="tab:orange", alpha=0.2)

figure.tight_layout()
figure.savefig("adaptive_branch_quadrature.pdf")
plt.close(figure)
