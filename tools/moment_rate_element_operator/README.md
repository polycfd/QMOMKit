# QMOM Spatial-Coupling Demonstrator

This standalone downstream example teaches two different ways a host spatial
discretization can couple raw-moment fields to QMOMKit. It is a manufactured
quadratic reference-element example, not a CG, DG, or FVM solver: it has no
mesh, numerical flux, physical-space derivative, physical-space transport,
physical-element mapping, or boundary-condition machinery.

The common host-only spatial mapping is:

```text
host element moment coefficients
        <->
local raw moments at spatial quadrature points
```

The host stores coefficients for raw-moment fields. It interpolates those
fields before every local QMOMKit call and projects local values back afterward.
QMOMKit never receives basis functions, coefficient vectors, a mass matrix, or
a spatial quadrature rule.

The manufactured field is a positive spatial scaling of one realizable
two-node QMOM moment sequence. The host never stores, transports, interpolates,
or projects QMOM nodes or weights; quadrature is an internal temporary
representation used only by QMOMKit's local phase-space operations.

## Shared timing and phase-space model

Both routes start from exactly the same moment coefficients and use the same
constant phase-space velocity, initial time, final time, fixed `dt`, node count,
and SSP-RK3 integrator. The example deliberately removes time-integration
order as a variable: the difference is where the integration occurs relative to
the host projection.

```text
Route A: QMOMKit rate -> projection -> host SSP-RK3
Route B: QMOMKit local SSP-RK3 evolution -> projection
```

The routes are different coupling strategies. The example does not claim that
their final coefficient vectors must be identical, because projection and local
phase-space evolution occur in a different order.

## Route A: host-integrated moment rates

```text
for each host time step:

    host SSP-RK3 stage
        element coefficients
        -> interpolate local raw moments
        -> QMOMKit instantaneous moment rates
        -> project rates
        -> coefficient RHS

    repeat with fresh local rates for all SSP-RK3 stages
```

The program calls `qmomkit::evaluate_qmom_advection_moment_rates(...)` at every
spatial quadrature point and at every SSP-RK3 stage.

> QMOMKit supplies local instantaneous phase-space moment rates; the host owns
> SSP-RK3 time integration.

A real CG/DG/FVM application would replace this manufactured coefficient update
with its native time-integration and physical-space machinery.

## Route B: QMOMKit-integrated finite evolution

```text
for each phase-space substep:

    element coefficients
    -> interpolate local raw moments

    for each spatial quadrature point:
        raw moments
        -> reconstruct local QMOM state
        -> QMOMKit SSP-RK3 evolution over dt
        -> evolved raw moments

    evolved local raw moments
    -> project
    -> updated element coefficients
```

The program constructs a temporary `QmomState` for each local raw-moment set,
then calls `qmomkit::advance_qmom_phase_space(...)`. Its node-advection option
is explicitly `NodeAdvectionScheme::ssp_rk3`; the source callback is explicitly
zero, so the manufactured model contains only phase-space advection.

> QMOMKit owns the local SSP-RK3 phase-space advance; the host owns the spatial
> representation and decides where this phase-space step appears in its outer
> operator splitting.

This example contains no physical-space transport. A real host solver supplies
that evolution independently and can place the finite local phase-space step
inside its own operator-splitting sequence.

## Build and run

Install QMOMKit first, then configure this standalone consumer against that
install prefix:

```sh
cmake -S tools/moment_rate_element_operator \
  -B tools/moment_rate_element_operator/build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$PWD/install"
cmake --build tools/moment_rate_element_operator/build
tools/moment_rate_element_operator/build/moment_rate_element_operator
```

The executable prints the shared initial coefficients and the final coefficient
vectors produced by Routes A and B. Their agreement is not a success criterion.
