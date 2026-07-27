# Output Design

This page documents the current QMOMKit design for output handling. It is a
developer reference for maintainers, contributors, and AI-assisted development
workflows.

For moment, tensor, phase-space, node, and weight conventions, see
`numerical_conventions.md`. For configuration parsing and validation, see
`configuration_design.md`.

## Role Of Output

QMOMKit output converts numerical state into reproducible external
representations. The output layer records moment and quadrature histories in a
stable, deterministic order so examples, validation workflows, and downstream
applications can compare results across runs.

Output is separate from numerical algorithms:

- QMOM and CQMOM reconstruction do not write files.
- Node advection, moment-source updates, and transport integrators do not write
  files.
- Writers do not schedule simulation steps, choose output cadence, or own
  physical model logic.

Applications and drivers decide when to construct writers and when to call
`write(...)`. This keeps the core numerical algorithms reusable in contexts
that do not use filesystem output.

## Output Architecture

The public output API lives in:

```text
include/qmomkit/output.hpp
```

The public concepts are:

```text
OutputOptions
MomentWriter
QuadratureWriter
```

`OutputOptions` is an application-controlled policy object. It stores:

- output directory;
- whether moment history output is requested;
- whether quadrature history output is requested;
- an optional output interval.

The numerical algorithms do not consume `OutputOptions` directly. A driver may
parse it from configuration, create directories, construct writers, and schedule
calls according to its own simulation loop.

`MomentWriter` writes raw moment histories. Construction opens exactly one text
file, writes one header row, and each `write(...)` call appends one time row.
The writer supports one-dimensional `MomentSet` output, dense
multidimensional `RawMomentTensor` compatibility output, advanced lean
`CqmomMomentData` output, and requested CQMOM moments reconstructed from
`CqmomQuadrature`.

`QuadratureWriter` writes quadrature histories. Construction opens exactly one
text file, writes one header row, and each `write(...)` call appends one time
row. The writer supports one-dimensional `Quadrature1D` or `QmomResult` output
and multidimensional `CqmomQuadrature` output.

Both writers write to caller-provided paths. The public writer API does not
choose application-level filenames. Examples may use names such as
`moments.txt` or `quadrature.txt`, but those names are driver choices rather
than writer-owned conventions.

Writers validate structural consistency before writing rows. Examples include
finite time values, finite output values, matching moment counts, matching
tensor extents, matching quadrature node counts, and matching CQMOM dimension
counts. Invalid rows are reported through exceptions rather than silently
truncated or reshaped.

## Moment History Output

Moment history output records native integer-order raw moments of the number
density function at explicit times supplied by the caller. Generalized moments
with real exponents are caller-evaluated post-processing quantities from
quadrature and are not standard `MomentWriter` output.

The text representation is space-separated. The first column is time. Remaining
columns are moment values:

```text
# time M[0] M[1] M[2] ...
```

Numeric values are currently written with `std::setprecision(17)`. This is a
fixed writer behavior, not a user-configurable option. It preserves
high-precision floating-point text representation for reproducible output in
ordinary validation workflows.

For one-dimensional QMOM moment output, the all-moment writer uses ascending raw
moment order:

```text
M[0], M[1], M[2], ...
```

Selected one-dimensional moments use the exact order supplied to
`MomentWriter::for_qmom(path, selected_orders)`. Selection is writer-local
metadata. It does not modify `MomentSet`, `QmomCase`, or any numerical
algorithm.

For dense CQMOM compatibility output, all-moment output uses the row-major
`RawMomentTensor` storage convention. Multi-indices are enumerated in
phase-space dimension order with the last index varying fastest:

```text
# time M[0,0] M[0,1] M[1,0] M[1,1]
```

Selected dense CQMOM moments use the exact `MultiIndex` order supplied to the
writer and are read from the supplied dense tensor.

Advanced lean CQMOM moment output is constructed from a `CqmomMomentLayout`
and writes `CqmomMomentData` rows. Columns use `layout.required_indices` order,
and values are read directly from `CqmomMomentData::values`; the writer does
not require dense tensor extents or `RawMomentTensor` storage.

Requested CQMOM moment output is constructed from a `RawMomentRequest` and
writes rows from `CqmomQuadrature`. Each requested raw moment is reconstructed
directly from the terminal quadrature nodes:

```text
m_alpha = sum_j W_j product_r xi_(r,j)^alpha_r
```

Requested output preserves request order and duplicate requested indices. It
does not create a dense `RawMomentTensor`.

## Quadrature History Output

Quadrature history output records quadrature weights and abscissas at explicit
times supplied by the caller.

For one-dimensional QMOM quadrature, each node is written as a
weight/abscissa pair in the quadrature node order:

```text
# time w0 x0 w1 x1 w2 x2
```

The writer accepts `Quadrature1D` and `QmomResult`. Callers are responsible for
checking numerical status before writing a result object. Output helpers do not
reinterpret a failed reconstruction as valid numerical data.

For CQMOM quadrature, each terminal node is written as its accumulated terminal
weight followed by abscissas in phase-space dimension order:

```text
# time w0 x0_0 x1_0 w1 x0_1 x1_1
```

Here `xD_N` means the abscissa in phase-space dimension `D` for terminal node
`N`. The terminal node order is the order exposed by
`CqmomQuadrature::weights()` and `CqmomQuadrature::abscissas()`. Fixed
rectangular quadratures use the public fixed hierarchy ordering. Adaptive
ragged quadratures use the deterministic terminal ordering represented by the
public quadrature object.

`QuadratureWriter` always writes the full quadrature. Node selection, node
filtering, restart files, visualization formats, checkpoints, parallel I/O, and
logging are outside the current output API.

## Configuration Interaction

Configuration can describe optional output settings:

```text
[output]
directory = output
write_moments = true
write_quadrature = true
interval = 0.1
```

Parsing this section creates an `OutputOptions` object. It does not modify
`QmomCase`, `CqmomDenseCase`, transport options, time integrators, moment
sources, or velocity callbacks.

If no output section is present, parsed output options use the default
application-controlled behavior:

```text
directory = .
write_moments = false
write_quadrature = false
interval = unset
```

The optional interval is cadence metadata for the driver. QMOMKit does not own
an output scheduler. A driver that uses the interval must decide how to align it
with its timestep sequence, final time, retry behavior, and validation output.

Configuration controls whether a driver may choose to create writers. It does
not define output data semantics. Moment and quadrature column ordering remain
the responsibility of the writer contracts and the numerical conventions.
Likewise, configuration may provide an output directory, but it does not choose
the filenames passed to `MomentWriter` or `QuadratureWriter`.

## Extension Principles

Output extensions should preserve reproducibility and keep file I/O separate
from numerical algorithms.

Stable rules for extending output are:

- preserve existing column ordering for current text outputs;
- keep raw moments in the documented `MomentSet` and `RawMomentTensor` order;
- keep quadrature weights paired with their corresponding abscissas;
- preserve phase-space dimension order for multidimensional abscissas;
- make writer selection metadata local to the writer, not to request or result
  objects;
- validate structural mismatches before writing rows;
- avoid hidden global output state or implicit scheduling inside algorithms;
- keep configuration parsing separate from writer construction and simulation
  control.

New output formats should document their data contract before becoming a public
API. If a format is intended only for one example or validation script, keep it
local to that consumer instead of adding it to the core package.
