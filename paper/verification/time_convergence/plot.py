#!/usr/bin/env python3
"""Plot the fixed QMOMKit time-convergence paper studies."""

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def reference_line(axis, x, y, order, label):
    axis.loglog(x, y[0] * (x / x[0]) ** order, "--", color="0.45", label=label)


rk_data = np.atleast_1d(np.genfromtxt("rk_convergence.csv", delimiter=",", names=True, dtype=None, encoding="utf-8"))
rk_figure, rk_axis = plt.subplots(figsize=(7.2, 5.0))
for scheme, order in {"ssp_rk3": 3, "rk23": 3, "rk45": 5}.items():
    rows = rk_data[rk_data["scheme"] == scheme]
    rk_axis.loglog(rows["step_size"], rows["absolute_error"], "o-", label=scheme)
    reference_line(rk_axis, rows["step_size"], rows["absolute_error"], order, f"order {order}")
rk_axis.set_xlabel("fixed step size")
rk_axis.set_ylabel("final absolute error")
rk_axis.grid(True, which="both", alpha=0.3)
rk_axis.legend(fontsize="small", ncol=2)
rk_figure.tight_layout()
rk_figure.savefig("rk_convergence.pdf")
plt.close(rk_figure)

splitting_data = np.atleast_1d(np.genfromtxt("splitting_convergence.csv", delimiter=",", names=True, dtype=None, encoding="utf-8"))
splitting_figure, splitting_axis = plt.subplots(figsize=(7.2, 5.0))
for mode, order in {"advsrc": 1, "srcadv": 1, "strang": 2}.items():
    rows = splitting_data[splitting_data["split_mode"] == mode]
    splitting_axis.loglog(rows["outer_step_size"], rows["max_scaled_moment_error"], "o-", label=mode)
    reference_line(splitting_axis, rows["outer_step_size"], rows["max_scaled_moment_error"], order, f"order {order}")
splitting_axis.set_xlabel("outer splitting step size")
splitting_axis.set_ylabel("maximum scaled moment error")
splitting_axis.grid(True, which="both", alpha=0.3)
splitting_axis.legend(fontsize="small", ncol=2)
splitting_figure.tight_layout()
splitting_figure.savefig("splitting_convergence.pdf")
plt.close(splitting_figure)
