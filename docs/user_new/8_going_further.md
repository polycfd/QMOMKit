# Going further

## Examples

The examples under [`examples/`](../../examples/) are small standalone programs,
built against the installed package as described in
[Installation](1_installation.md). Each has a README with its model, commands,
and expected output. In order of increasing scope:

| Example | Shows |
| --- | --- |
| **Moments and quadrature** | |
| [`qmom_explicit_moments`](../../examples/qmom_explicit_moments/README.md) | Reconstructing a quadrature from given moments. |
| [`qmom_distribution_initialization`](../../examples/qmom_distribution_initialization/README.md) | Moments from a log-normal distribution. |
| [`qmom_mixture_initialization`](../../examples/qmom_mixture_initialization/README.md) | Moments from a mixture with absolute component densities. |
| [`realizability_check`](../../examples/realizability_check/README.md) | Checking moments against a support. |
| [`cqmom_independent_product_initialization`](../../examples/cqmom_independent_product_initialization/README.md) | A CQMOM case from a configuration file, and the conditioning order. |
| **Evolving in time** | |
| [`qmom_1d_minimal`](../../examples/qmom_1d_minimal/README.md) | The smallest QMOM time loop, with a velocity and a fragmentation source. |
| [`qmom_1d_aerosol_evaporation`](../../examples/qmom_1d_aerosol_evaporation/README.md) | Node motion alone: evaporating aerosol. |
| [`qmom_1d_breakage_advection`](../../examples/qmom_1d_breakage_advection/README.md) | Node motion and a size-dependent breakage source written over the nodes. |
| [`cqmom_2d_minimal`](../../examples/cqmom_2d_minimal/README.md) | The CQMOM time loop with velocity, source, and boundary contribution. |
| [`cqmom_2d_adaptive_node_counts`](../../examples/cqmom_2d_adaptive_node_counts/README.md) | Branches losing and regaining nodes during a run. |
| **Coupling** | |
| [`qmom_state_candidate`](../../examples/qmom_state_candidate/README.md) | Creating states from moments your code changed. |
| [`qmom_1d_explicit_strang_splitting`](../../examples/qmom_1d_explicit_strang_splitting/README.md) | Calling the stages of a step yourself. |
| [`cqmom_2d_batch_host_data`](../../examples/cqmom_2d_batch_host_data/README.md) | Many states at once, with host data per state. |
| **With a Monte Carlo reference** | |
| [`qmom_1d_breakage_advection_with-mc`](../../examples/qmom_1d_breakage_advection_with-mc/README.md) | QMOM breakage and advection beside a particle simulation. |
| [`cqmom_2d_breakage_advection_with-mc`](../../examples/cqmom_2d_breakage_advection_with-mc/README.md) | The same in two dimensions. |
| **Physical models** (with a Monte Carlo reference) | |
| [`cqmom_2d_bubbles_cavitation`](../../examples/cqmom_2d_bubbles_cavitation/README.md) | Nonlinear node motion of a cavitating bubble population. |
| [`cqmom_3d_droplets_condensation-nucleation`](../../examples/cqmom_3d_droplets_condensation-nucleation/README.md) | Condensation and nucleation as a boundary contribution. |
| [`cqmom_4d_particles_breakage_recoil`](../../examples/cqmom_4d_particles_breakage_recoil/README.md) | A breakage source evaluated over the nodes of the stage quadrature. |
| [`cqmom_5d_particles_breakage_heat-transfer`](../../examples/cqmom_5d_particles_breakage_heat-transfer/README.md) | A source needing fractional moments. |
| [`cqmom_5d_non-spherical-particles`](../../examples/cqmom_5d_non-spherical-particles/README.md) | Bounds on a coordinate, and time-varying host data. |

## Monte Carlo companion

The optional `QMOMKitMonteCarlo` package samples a finite population from a
distribution, estimates its moments, evolves the sampled particles with the same
velocity, and compares moment histories. Starting QMOMKit and the particle
simulation from the same sample removes the initial sampling difference, so
that what remains comes from the evolution. A Monte Carlo result carries sampling
noise: it is a reference, not an exact solution. See
[tools/monte_carlo/README.md](../../tools/monte_carlo/README.md).

## Citing QMOMKit

If QMOMKit contributes to your work, please cite it as given in the
[repository README](../../README.md#how-to-cite).

## Developer documentation

The [developer documentation](../developer/README.md) describes the algorithms,
numerical conventions, and design of QMOMKit, and how to test and extend it.

## Coupling to host codes

Examples of couplings of QMOMKit to host simulation codes are maintained in
[QMOMKit-Couplings](https://github.com/polycfd/QMOMKit-Couplings).
