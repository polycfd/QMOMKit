# CQMOM 2D minimal example

This example demonstrates a small two-dimensional CQMOM phase-space evolution
workflow using QMOMKit's typed public C++ API. All case definition, callbacks,
time controls, and output setup are in `main.cpp`.

## What it demonstrates

The example follows this state-centric workflow:

```text
CqmomDenseCase -> CqmomCase -> initial_state() -> CqmomState
                                                   |
                                                   v
                                  advance_cqmom_phase_space(...)
                                                   |
                                                   v
                                           state.quadrature()
```

The two phase-space coordinates are generic coordinates `x0` and `x1`. The
case uses two nodes in each coordinate and conditioning order `{0, 1}`. It is
constructed from a dense `4 x 4` raw-moment tensor, with rows indexed by the
order in `x0` and columns by the order in `x1`:

$$
\begin{bmatrix}
1 & 0 & 0.25 & 0 \\
2 & 0 & 0.5 & 0 \\
4.25 & 0 & 1.0625 & 0 \\
9.5 & 0 & 2.375 & 0
\end{bmatrix}.
$$

## Model

The phase-space velocity callback is

$$
\frac{\mathrm{d}x_0}{\mathrm{d}t} = 1-x_0, \qquad
\frac{\mathrm{d}x_1}{\mathrm{d}t} = -x_1.
$$

The moment source callback applies binary-fragmentation-style rates along
`x0`:

$$
\frac{\mathrm{d}M_{\boldsymbol{\alpha}}}{\mathrm{d}t}
= 0.5\left(2^{1-\alpha_0}-1\right)M_{\boldsymbol{\alpha}}.
$$

The boundary contribution adds

$$
0.1\,(0.1)^{\alpha_0}
$$

to each supplied moment rate. The run starts at `t = 0`, ends at `t = 1`,
and advances with `dt = 0.01`. Strang splitting is used with SSP-RK3 node
advection. The initial state and subsequent states are written at an output
interval of `0.1`.

## Build and run

From this example directory, with QMOMKit installed in `../../install`:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./cqmom_2d_minimal_example
```

## Output files

On success the executable prints one line confirming where the results went, and
writes:

- `output/moment_history.txt` — raw-moment history;
- `output/quadrature_history.txt` — current terminal CQMOM weights and coordinates.

The moment history uses QMOMKit's state-aware moment writer, and the
quadrature history is written from `state.quadrature()`. No generalized
moments are requested.
