# Galerkin Element Coupling Demonstrator

This is a manufactured Galerkin element coupling demonstrator. It uses one
small reference element to show how raw-moment fields represented by spatial
basis coefficients can call QMOMKit's local phase-space APIs.

It is not a complete CG or DG solver, an FVM solver, a general finite-element
framework, or a physical-space transport example. There is no mesh, global
assembly, numerical flux, physical-space derivative, physical-element mapping,
or boundary-condition machinery.

## The general coupling principle

QMOMKit does not prescribe a physical-space field representation. In general,
the coupling boundary is simply:

```text
host spatial representation
-> local raw-moment state
-> QMOMKit local phase-space operation
-> host spatial representation
```

This demonstrator shows the Galerkin version of that boundary:

```text
spatial basis coefficients
-> interpolation to spatial quadrature points
-> local QMOMKit calls
-> mass-matrix projection
-> spatial basis coefficients
```

That mapping is directly relevant to CG and DG implementations. It is not a
generic CG/DG/FVM procedure. A typical cell-centred finite-volume method
already stores one local moment state per cell, so its connection to QMOMKit neither 
require basis interpolation nor mass-matrix projection. Higher-order finite-volume 
reconstruction or quadrature is host-specific and is not illustrated here.

## Galerkin representation in this example

Conceptually, `coefficient[i][k]` is the coefficient of raw moment order `k`
in spatial basis function `i`. QMOMKit instead works on one local raw-moment
set at a time. The Galerkin host therefore evaluates its moment fields at each
spatial quadrature point; every such point supplies an independent local
phase-space problem.

After QMOMKit returns local rates or evolved local moments, the host maps those
values back to its basis. The weighted projection forms the basis inner
products, and the inverse element mass matrix converts them to coefficient
values.

The spatial quadrature points here are Gauss points of the Galerkin reference
element. They are unrelated to QMOM quadrature nodes and weights, which
represent the distribution in phase space. The Galerkin solution stores only
moment fields in physical space; it never spatially stores, transports,
interpolates, or projects QMOM nodes or weights.

The manufactured field is a positive spatial scaling of one realizable,
two-node one-dimensional QMOM moment sequence.

## Shared model and timing

Both routes start from the same coefficient field and use the same phase-space
velocity, node count, time interval, fixed `dt`, and SSP-RK3 scheme. This keeps
the time integrator fixed so the example isolates the coupling strategy:

```text
Route A: local rate -> spatial projection -> host SSP-RK3
Route B: local SSP-RK3 evolution -> spatial projection
```

The routes are different discrete coupling strategies. The example does not
use agreement of their final coefficients as a validation criterion.

## Route A: instantaneous rates, host SSP-RK3

```text
for every SSP-RK3 stage:

    spatial coefficients
    -> interpolate local raw moments
    -> QMOMKit instantaneous moment rates
    -> project local rates
    -> coefficient-space RHS
    -> host SSP-RK3 stage update
```

`qmomkit::evaluate_qmom_advection_moment_rates(...)` supplies the instantaneous
local phase-space rate at every spatial quadrature point. The Galerkin host
projects those rates into its basis and performs the SSP-RK3 time integration.
Each RK stage needs a fresh interpolate -> rate -> project operation because
its stage coefficients define different local moment states.

## Route B: finite local evolution

```text
for every time step:

    spatial coefficients
    -> interpolate local raw moments

    at every spatial quadrature point:
        raw moments
        -> construct local QMOM state
        -> QMOMKit SSP-RK3 evolution over dt
        -> evolved raw moments

    evolved local raw moments
    -> project back to the spatial basis
    -> updated coefficients
```

For each spatial quadrature point, the program constructs a `QmomState` and
calls `qmomkit::advance_qmom_phase_space(...)` with
`NodeAdvectionScheme::ssp_rk3`. The local QMOM state contains its corresponding
QMOM quadrature representation for the phase-space calculation. After the
substep, the Galerkin host receives evolved raw moments and projects only those
moments back into its spatial representation.

QMOMKit performs the finite local phase-space advance. The Galerkin host
chooses where that substep appears in an outer operator-splitting sequence.

## Coupling and performance trade-offs

Route A crosses the Galerkin/QMOMKit boundary three times per SSP-RK3 step:
each stage needs interpolation, local moment-rate evaluation, and projection.
That gives the host a conventional semidiscrete RHS, which can participate
directly in an application-managed method-of-lines integrator alongside other
operators. Its cost is repeating the Galerkin mapping and local reconstruction
or rate evaluation at every stage.

Route B maps into local states once around the full local phase-space substep:

```text
interpolate -> QMOMKit finite SSP-RK3 evolution -> project
```

QMOMKit executes the local SSP-RK3 stages before returning evolved moments.
This can avoid repeated Galerkin interpolation and projection, and can reduce
repeated reconstruction work. It is especially natural when physical-space and
phase-space evolution are already operator-split. The trade-off is that the
host uses a split finite evolution instead of adding the phase-space term to a
single global semidiscrete RHS.

## Build and run

Install QMOMKit first, then configure this standalone consumer against that
install prefix:

```sh
cmake -S tools/galerkin_element_coupling \
  -B tools/galerkin_element_coupling/build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$PWD/install"
cmake --build tools/galerkin_element_coupling/build
tools/galerkin_element_coupling/build/galerkin_element_coupling
```
