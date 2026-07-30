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

Plot all moments from a `MomentWriter` file:

```bash
python tools/plotting/plot_moments.py output/moments.txt
```

Plot selected moment columns:

```bash
python tools/plotting/plot_moments.py output/moments.txt --moments 'M[0]' 'M[2]'
```

## Quadrature Histories

Plot a one-dimensional QMOM `QuadratureWriter` file:

```bash
python tools/plotting/plot_quadrature.py output/quadrature.txt
```

Plot a CQMOM quadrature history:

```bash
python tools/plotting/plot_quadrature.py output/quadrature_history.txt
```

For QMOM output, each row is one active node at one output time. Adaptive QMOM
files may contain different active node counts at different times, and inactive
nodes are not padded. For CQMOM output, each row is one terminal node with its
recursive path. The script plots weights and abscissas by public node index and
adds a simple final-time CQMOM node visualization. Multidimensional CQMOM output
is shown in the first two dimensions.

These scripts deliberately avoid a general projection framework, trajectory
analysis, or specialized interpretation of CQMOM nodes. They read the column
labels and numeric rows directly from the native output files.
