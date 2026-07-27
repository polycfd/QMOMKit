# QMOMKit Plotting Helpers

These are lightweight matplotlib helpers for quick visualization of
QMOMKit-native text output. They are intended as small starting points that
researchers can read, copy, and adapt.

The scripts consume files written by `qmomkit::MomentWriter` and
`qmomkit::QuadratureWriter`. Those files are whitespace-delimited text tables
with commented headers, not CSV files:

```text
# time M[0] M[1] ...
# time w0 x0 w1 x1 ...
# time w0 x0_0 x1_0 w1 x0_1 x1_1 ...
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

For QMOM output, the script plots weights versus time and abscissas versus
time. For CQMOM output, it plots terminal weights versus time and a simple
final-time node visualization. Multidimensional CQMOM output is shown in the
first two dimensions.

These scripts deliberately avoid a general projection framework, trajectory
analysis, or specialized interpretation of CQMOM nodes. They read the column
labels and numeric rows directly from the native output files.
