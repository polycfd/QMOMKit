# 1D QMOM breakage and advection

This example evolves a three-node `QmomState` with
`advance_qmom_phase_space(...)`, which combines two callbacks by operator
splitting: a phase-space velocity that moves the quadrature nodes, and a
breakage source that returns moment rates. The source evaluates its rates over
the nodes of the stage quadrature with `context.for_each_node(...)`, which
closes moments of non-integer order that the state does not track.

## Model

The phase-space coordinate is the particle volume $V$ and $n(V)$ is the number
density, with raw moments $m_k = \int V^k\,n(V)\,\mathrm{d}V$. Volume and time
$t$ are dimensionless. Volumes relax toward $V = 1$,

```math
\frac{\mathrm{d}V}{\mathrm{d}t} = 1 - V,
```

and particles break into two of half the volume, $V \to V/2$, at the rate

```math
a(V) = a_\mathrm{ref}\left(\frac{V}{V_\mathrm{ref}}\right)^{p},
\qquad a_\mathrm{ref} = 0.5,\quad V_\mathrm{ref} = 1,\quad p = \tfrac{3}{2}.
```

Breakage conserves the total volume $m_1$, and the moments change as

```math
\frac{\mathrm{d}m_k}{\mathrm{d}t}
= \left(2^{1-k}-1\right)\int a(V)\,V^k\,n(V)\,\mathrm{d}V
\approx \left(2^{1-k}-1\right)\sum_i w_i\,a(V_i)\,V_i^k,
```

where $w_i$ and $V_i$ are the weights and abscissas of the quadrature that
QMOMKit reconstructs from the moments of each Runge–Kutta stage. The exact
right-hand side needs $m_{k+3/2}$, which the state does not track. With $p = 0$
the rate is constant and the source is closed in the moments, the model of
[qmom_1d_breakage_advection_with-mc](../qmom_1d_breakage_advection_with-mc/README.md).
$a_\mathrm{ref}$, $V_\mathrm{ref}$ and $p$ are set in
[physics.cpp](./physics.cpp), which defines both callbacks.

[configuration.qmom](./configuration.qmom) defines the initial number density
$0.5\,\mathcal{N}(1, 0.1^2) + 0.5\,\mathcal{N}(3, 0.2^2)$ and three nodes
(`node_counts = 3`). For $n$ nodes, `maximum_order` must be at least $2n-1$,
or omitted to use that minimum; here it is 5. The file has no
`[time_integration]` section, so each step uses the default Strang splitting:
half a step of node advection, the source step, half a step of node
advection, as described in [Evolving in
time](../../docs/user/5_evolving_in_time.md). [main.cpp](./main.cpp) advances
from $t = 0$ to $1$ with steps of $0.01$.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./qmom_1d_breakage_advection_example
```

## Output

The example prints the initial and final number density $m_0$ and mean volume
$m_1/m_0$, and the final quadrature. It writes the moments and the quadrature
every $0.1$ to `build/output/moments.txt` and `build/output/quadrature.txt`, as
set in the `[output]` section of `configuration.qmom`.

## What to expect

Breakage raises $m_0$. Large particles break faster: the mode near $V = 3$
initially breaks about five times as fast as the mode near $V = 1$. Breakage
and advection both shrink the large particles, so the growth of $m_0$ slows: it
gains about 0.14 in the first $0.1$ and about 0.09 in the last. By $t = 1$,
$m_0$ has grown from 1 to about 2.04, against about 1.65 for the constant rate
$a_\mathrm{ref}$, and the mean volume has fallen from 2 to about 0.87. In
`moments.txt`, $m_5$ falls from about 127 to about 2.8. The third node of the
final quadrature, near $V = 1.7$, carries the least weight.
