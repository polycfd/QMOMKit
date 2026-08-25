# 4D spherical-breakup manuscript example

This is the manuscript version of
`examples/cqmom_spherical_breakup_quadrature_4d`. It is an independent copy
with identical phase space `(d,u,v,w)`, particle and breakup/recoil models,
sampled initialization, seeds, and CQMOM evolution settings. It is a QMOMKit
capability demonstration of a source-term closure strategy, not a physical
validation case for a specific breakup model.

Diameter-dependent recoil creates a nonlinear daughter-state mapping and mixed
diameter-velocity moments outside the transported finite moment set. Terminal
quadrature therefore directly evaluates the nonlinear population
transformation.

The physical model and CQMOM evolution settings match the user example. For a
controlled comparison, the paper MC reference intentionally reuses the
configured adaptive integration settings.

Build this case through the parent paper-example project. Its one executable
creates the deterministic sampled population once, initializes CQMOM from its
raw tensor, and evolves that same population as the MC reference:

```bash
cmake -S paper/examples -B <build-dir> -DCMAKE_PREFIX_PATH=<install-prefix>
cmake --build <build-dir>
cd <build-dir>/cqmom_spherical_breakup_quadrature_4d
../qmomkit_paper_cqmom_spherical_breakup_quadrature_4d
```

`output/moment_comparison.csv` uses the paper comparison columns
`time,quantity,qmomkit,monte_carlo`. It stores the raw moments
`M[0,0,0,0]`, `M[3,0,0,0]`, `M[1,0,0,0]`, `M[2,0,0,0]`, `M[0,1,0,0]`,
`M[0,2,0,0]`, `M[2,1,0,0]`, `M[2,2,0,0]`, and `M[6,0,0,0]` for later
post-processing. No normalization is applied. Native raw-moment and
quadrature histories are also written according to `configuration.qmom`.

The copied configuration controls CQMOM. Its diameter marginal is a trimodal
lognormal mixture: component number densities `0.40`, `0.50`, and `0.10`,
median diameters `25`, `100`, and `400` micrometres, and the original
log-space standard deviation `0.12`. The component number densities are normalized by
the independent-product initialization before applying the global number
density of one. The case uses 10,000 initialization samples with seed `123`,
breakup seed `321`, and end time `0.1`.
