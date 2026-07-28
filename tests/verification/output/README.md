# Output Verification

This directory verifies the public QMOMKit text output workflow. The output
layer is intentionally application-controlled: numerical algorithms compute
QMOM or CQMOM results, `OutputOptions` records whether a driver should write
moments or quadrature and where files should be placed, and `MomentWriter` or
`QuadratureWriter` write exactly the caller-provided paths.

`output_workflows.cpp` is the canonical output verification case. It uses a
two-node analytical QMOM configuration and the independent 2D analytical CQMOM
configuration, writes moment and quadrature files, and reads those files back
as whitespace-delimited text.

The case documents and verifies:

- enabled output creates the expected driver-chosen files;
- disabled `OutputOptions` do not create files in the representative driver
  workflow;
- headers are commented and describe the public column ordering;
- data rows have the same number of columns as the headers;
- moment values and quadrature values match the computed public result or an
  independent analytical product quadrature;
- CQMOM requested moment output preserves requested multi-index order;
- CQMOM quadrature output writes one row per active terminal node, including the
  local terminal index, recursive path, weight, and phase-space coordinates.

The writer API does not own filenames, output cadence, restart formats, or
visualization formats. Those remain application responsibilities.
