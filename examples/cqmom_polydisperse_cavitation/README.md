# CQMOM polydisperse cavitation example

## Purpose and workflows

This standalone example evolves a spatially uniform bubble population in the
two-dimensional cavitation phase space

$$
\boldsymbol{\xi} = (R,U),
$$

where $R$ is the instantaneous dimensionless bubble radius and
$U=\dot R$ is the dimensionless radial velocity. 

The example contains two executables:

- `cqmom_polydisperse_cavitation_example`: initializes CQMOM from empirical raw
  moments of a deterministic finite sampled population and advances the CQMOM
  quadrature through the continuous cavitation velocity.
- `cqmom_polydisperse_cavitation_mc_example`: evolves the same deterministic
  finite sampled population directly with the Monte Carlo companion utilities.

Both workflows use the same `configuration.qmom`, distribution parameters,
sample count, initialization seed, coordinate order, and continuous
phase-space velocity. The Monte Carlo executable is a direct sampled-population
companion for this example; it is not a separate physical model and no
automatic CQMOM/MC comparison is performed.

The CQMOM workflow is:

```text
configuration.qmom
        |
        v
fixed-seed sampled initialization
        |
        v
empirical raw moments
        |
        v
CQMOM case
        |
        v
compute_cqmom(...) -> initial quadrature
        |
        v
continuous CQMOM phase-space evolution
```

The direct Monte Carlo workflow is:

```text
configuration.qmom
        |
        v
fixed-seed sampled initialization
        |
        v
finite sampled population
        |
        v
continuous sample evolution
```

## Phase-space definition and raw moments

The configured phase-space order is:

```text
(R, U)
```

The same order is used by `configuration.qmom`, the sampled ensemble, the
phase-space velocity callback, CQMOM moment exponents, quadrature abscissas,
and output files. Coordinate names are example metadata; QMOMKit itself treats
them as generic phase-space coordinates.

For a number density function $f(R,U)$, the raw unnormalized moments are

$$
\begin{aligned}
M_{i,j}
&=
\int R^i U^j f(R,U)\,\mathrm{d}R\,\mathrm{d}U.
\end{aligned}
$$

Here $i$ is the order in radius and $j$ is the order in radial velocity.
The zeroth moment $M_{0,0}$ is the total number density and is not normalized
away by QMOMKit.

The finite sampled initialization estimates the same raw moments as

$$
\begin{aligned}
M_{i,j}^{(N)}
&=
\frac{n_0}{N}
\sum_{p=1}^{N}
R_p^i U_p^j,
\end{aligned}
$$

where $N$ is the sample count, $n_0$ is the configured total number density,
$(R_p,U_p)$ is sample $p$, and $n_0/N$ is the constant statistical weight
carried by each sample. The empirical moments $M_{i,j}^{(N)}$ initialize the
CQMOM case, while the same samples $(R_p,U_p)$ initialize the direct Monte
Carlo companion.

## Cavitation model

The continuous callback in `physics.cpp` advances each point in phase space
with

$$
\frac{\mathrm{d}R}{\mathrm{d}t}=U,
$$

and

$$
\begin{aligned}
\frac{\mathrm{d}U}{\mathrm{d}t}
&=
\frac{1}{R}
\left[
\left(\frac{R_o}{R}\right)^{3\gamma}
-\frac{1}{C_p}
-\frac{2}{\mathrm{We}\,R_o}
\left(
\frac{R_o}{R}
-
\left(\frac{R_o}{R}\right)^{3\gamma}
\right)
-\frac{3}{2}U^2
-\frac{4}{\mathrm{Re}}\frac{U}{R}
\right].
\end{aligned}
$$

The symbols are:

- $R$: instantaneous dimensionless bubble radius;
- $U$: dimensionless radial velocity;
- $t$: dimensionless time used by the example driver;
- $R_o$: reference radius;
- $C_p$: pressure coefficient;
- $\gamma$: heat-capacity ratio;
- $\mathrm{Re}$: Reynolds number parameter;
- $\mathrm{We}$: Weber number parameter.

The implemented constants are compiled into the example physics callback:

| Symbol | Code name | Value | Role |
|---|---|---:|---|
| $R_o$ | `reference_radius` | `1.0` | reference radius in the radius ratio |
| $C_p$ | `pressure_coefficient` | `0.3` | pressure coefficient in the forcing term |
| $\gamma$ | `heat_capacity_ratio` | `1.4` | exponent coefficient for compression |
| $\mathrm{Re}$ | `reynolds_number` | `10.0` | viscous damping parameter |
| $\mathrm{We}$ | `weber_number` | `1.39` | surface-tension parameter |

These are example-local physics parameters, not entries in `configuration.qmom`
and not generic QMOMKit settings.

## Initialization shared by CQMOM and Monte Carlo

The supplied initialization is an independent product distribution with total
number density

$$
n_0 = 1.
$$

The configured marginals are:

- $R$ is log-normal with `mu_log = -0.02` and `sigma_log = 0.2`.
- $U$ is normal with `mean = 0.0` and `standard_deviation = 0.2`.

For the implemented log-normal parameterization,

$$
\begin{aligned}
\mathbb{E}[R]
&=
\exp\left(
\mu_{\log}
+
\frac{\sigma_{\log}^2}{2}
\right)
=
\exp\left(-0.02+\frac{0.2^2}{2}\right)
=
\exp(0)
=
1.
\end{aligned}
$$

Both executables draw

```text
sample_count = 10000
initialization_seed = 12345
```

through the shared sampling helper. The CQMOM executable computes empirical raw
moments from that finite sample and builds a sampled-moment CQMOM case. The
direct Monte Carlo executable evolves the sampled ensemble itself. Analytical
moments of the configured distributions are not used as the production CQMOM
initial condition in this example.

## Numerical workflow and configuration

The numerical CQMOM hierarchy is selected in `configuration.qmom`:

```text
[phase_space]
names = R, U
node_counts = 2, 2
conditioning_order = 0, 1

[initial]
type = independent_product
number_density = 1.0
maximum_orders = 3, 3
```

The `node_counts = 2, 2` setting gives two nodes in $R$, two conditional
nodes in $U$, and four terminal CQMOM nodes in the full default topology.
The `conditioning_order = 0, 1` setting reconstructs $R$ first and $U$
second. The `maximum_orders = 3, 3` setting requests sampled tensor moments
with orders $0,\ldots,3$ in each configured phase-space dimension. The driver
then uses the hierarchy-driven CQMOM moment layout for reconstruction and
output; it does not impose a paper-specific restricted moment list.

The driver time controls are compiled local constants:

```text
initial_time = 0.0
final_time = 13.9
time_step = 0.01
```

The output settings are read from:

```text
[output]
directory = output
write_moments = true
write_quadrature = true
interval = 0.1
```

`configuration.qmom` does not contain an explicit `[time_integration]` section.
The CQMOM executable loads the default transport options and advances the
accepted CQMOM state with `qmomkit::advance_cqmom_continuous(...)`; the default
node-advection scheme is `ssp_rk3`. During continuous CQMOM phase-space
advection, quadrature abscissas move and weights are preserved over the
continuous advancement; raw moments change because the abscissas change.

The direct Monte Carlo executable advances every sample with
`qmomkit::monte_carlo::evolve_ensemble(...)` and the same cavitation velocity
callback. It constructs default `qmomkit::AdaptiveRkOptions`, whose default
tolerances are `relative_tolerance = 1.0e-6` and
`absolute_tolerance = 1.0e-9` and whose initial step defaults to `1.0e-6 s`;
the example does not override the initial, minimum, or maximum step for the
Monte Carlo path.

## Build

This example is a standalone installed-package consumer. It requires both
QMOMKit and the optional QMOMKitMonteCarlo companion package.

From the repository root, install both packages into the shared
repository-local `install/` prefix:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"

cmake --build build
cmake --install build

cmake -S tools/monte_carlo -B tools/monte_carlo/build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install" \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"

cmake --build tools/monte_carlo/build
cmake --install tools/monte_carlo/build
```

Then configure and build this example from its directory:

```bash
cd examples/cqmom_polydisperse_cavitation

cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
```

This builds:

- `cqmom_polydisperse_cavitation_example`
- `cqmom_polydisperse_cavitation_mc_example`

CMake copies `configuration.qmom` into the example build directory.

## Run the CQMOM workflow

Run from the example build directory so generated output stays there:

```bash
cd build
./cqmom_polydisperse_cavitation_example
```

The executable prints sample metadata, time controls, selected initial and
final raw moments, final terminal quadrature nodes and weights, and output
paths.

## Run the direct Monte Carlo workflow

From the same build directory, run:

```bash
./cqmom_polydisperse_cavitation_mc_example
```

The executable prints sample metadata, time controls, selected initial and
final raw moments, final sample count, and its output path.

## Output files

With the default output configuration, the CQMOM executable writes:

```text
output/initial_moments.txt
output/moment_history.txt
output/initial_quadrature.txt
output/quadrature_history.txt
```

`initial_moments.txt` contains raw moments reconstructed from the initial CQMOM
quadrature built from the sampled initialization. `moment_history.txt` contains
the deterministic CQMOM raw-moment history at the configured output interval.
Both files are whitespace-delimited `MomentWriter` tables with a commented
header containing `time` followed by moment labels such as `M[0,0]`.

`initial_quadrature.txt` contains the initial terminal CQMOM nodes.
`quadrature_history.txt` contains the terminal CQMOM quadrature history. These
`QuadratureWriter` tables include time, public node index, conditioning path
entries, terminal weight, and abscissas ordered as $R$ and $U$.

The direct Monte Carlo executable writes:

```text
output/moment_history_mc.txt
```

`moment_history_mc.txt` contains Monte Carlo estimated raw moments computed
from the evolving finite sample with the constant per-sample weight $n_0/N$.
It uses the same phase-space moment ordering as the sampled moment tensor.

## Plotting

The repository provides generic plotting helpers for QMOMKit-native text
output. From the repository root, after running the example, plot moment
histories with:

```bash
python tools/plotting/plot_moments.py examples/cqmom_polydisperse_cavitation/build/output/moment_history.txt
python tools/plotting/plot_moments.py examples/cqmom_polydisperse_cavitation/build/output/moment_history_mc.txt
```

Selected moment columns can be requested by label:

```bash
python tools/plotting/plot_moments.py examples/cqmom_polydisperse_cavitation/build/output/moment_history.txt --moments 'M[0,0]' 'M[1,0]' 'M[0,1]'
```

Plot the CQMOM quadrature history with:

```bash
python tools/plotting/plot_quadrature.py examples/cqmom_polydisperse_cavitation/build/output/quadrature_history.txt
```

The plotting scripts read the existing whitespace-delimited output files and
write PNG files next to the input file.

## Expected qualitative behavior

The pressure coefficient in the implemented ODE excites radius and
radial-velocity evolution from the initially polydisperse distribution.
Surface-tension and viscous terms, represented by the Weber and Reynolds
parameters, influence the oscillation frequency and damping. CQMOM represents
the evolving population with a finite conditional quadrature, while the direct
Monte Carlo companion evolves the finite initial sample directly.

The CQMOM and Monte Carlo initial moments agree up to sampled-moment and
reconstruction precision because they come from the same deterministic finite
sample. Later differences reflect the finite CQMOM representation, finite
Monte Carlo sample, and numerical integration choices; exact agreement during
evolution is not expected.

## Interpretation and limitations

This example is a compact continuous phase-space cavitation workflow, not a
complete cavitating-flow solver. It does not include physical-space advection,
mesh coupling, gas-liquid feedback, collapse removal, breakup, coalescence,
discontinuous sources, or CHyQMOM. The implemented constants are fixed
example-local physics parameters, and the `.qmom` file describes numerical
configuration rather than application physics.

The Monte Carlo executable is useful for inspecting the direct evolution of
the same finite sampled population used for CQMOM initialization. It is not an
independent analytical reference, and the default `2 x 2` CQMOM topology is not
claimed to be converged or optimal.
