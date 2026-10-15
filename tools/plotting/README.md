# QMOMKit Plotting Helpers

These are lightweight matplotlib helpers for quick visualization of
QMOMKit-native text output. They are intended as small starting points that
researchers can read, copy, and adapt.

The scripts consume files written by `qmomkit::MomentWriter` and
`qmomkit::QuadratureWriter`. Those files are whitespace-delimited text tables
with commented headers, not CSV files:

```text
# time M[0] M[1] ...
# time node_index weight x0
# time node_index path_0 ... weight x0 ...
```

## Requirements

Only `numpy` and `matplotlib` are required.

## Moment Histories

Plot the moments of a `MomentWriter` file. A file with at most ten moment
columns is plotted in full; a larger one shows the raw moments of total order 0
and 1, the number density and the first moments, unless `--moments` selects
others:

```bash
python3 tools/plotting/plot_moments.py output/moments.txt
```

Plot selected moment columns:

```bash
python3 tools/plotting/plot_moments.py output/moments.txt --moments 'M[0]' 'M[2]'
```

Tables whose columns are not `M[...]` or `G[...]` moments are refused. The
figure is written next to the input as `<stem>_moments.png`, or to `--output`.

## Quadrature Histories

Plot a QMOM or CQMOM `QuadratureWriter` file; the script tells them apart by
the header:

```bash
python3 tools/plotting/plot_quadrature.py output/quadrature_history.txt
```

A QMOM file gives `<stem>.png` next to the input, a CQMOM file
`<stem>_history.png` and `<stem>_final.png`; `--output` sets the path or
prefix instead.

For QMOM output, each row is one active node at one output time. Adaptive QMOM
files may contain different active node counts at different times, and inactive
nodes are not padded. For CQMOM output, each row is one terminal node with its
recursive path. The script plots weights and abscissas by `node_index` and
adds a simple final-time CQMOM node visualization. When the number of nodes
changes between output times, an index can refer to a different node at the
next time, so the points are then shown without connecting lines.
Multidimensional CQMOM output is shown in the first two dimensions.
