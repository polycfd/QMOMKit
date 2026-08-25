# 5D spherical-breakup manuscript example

This is the manuscript version of
`examples/cqmom_spherical_breakup_moments_5d`. It is an independent copy with
identical phase space `(V,T,u,v,w)`, particle and breakup/temperature models,
sampled initialization, seeds, and CQMOM evolution settings. It is a QMOMKit
capability demonstration of a source-term closure strategy, not a physical
validation case for a specific breakup model.

The physical model and CQMOM evolution settings match the user example. For a
controlled comparison, the paper MC reference intentionally reuses the
configured adaptive integration settings.

Build this case through the parent paper-example project. Its one executable
creates the deterministic sampled population once, initializes CQMOM from its
raw tensor, and evolves that same population as the MC reference:

```bash
cmake -S paper/examples -B <build-dir> -DCMAKE_PREFIX_PATH=<install-prefix>
cmake --build <build-dir>
cd <build-dir>/cqmom_spherical_breakup_moments_5d
../qmomkit_paper_cqmom_spherical_breakup_moments_5d
```

`output/moment_comparison.csv` uses the paper comparison columns
`time,quantity,qmomkit,monte_carlo`. It stores the raw moments
`M[0,0,0,0,0]`, `M[1,0,0,0,0]`, `M[0,1,0,0,0]`, `M[0,2,0,0,0]`,
`M[1,1,0,0,0]`, `M[1/3,0,0,0,0]`, `M[4/3,0,0,0,0]`,
`M[1/3,1,0,0,0]`, and `M[4/3,1,0,0,0]`. The fractional CQMOM values
are evaluated from the terminal quadrature through QMOMKit generalized
moments. No normalization is applied. Native raw-moment and quadrature
histories are also written according to `configuration.qmom`. Fractional
moments use mathematical notation here; their CSV identifiers use names such
as `M_1over3_*` and `M_4over3_*`.

`output/generalized_moment_history_mc.txt` contains the two direct-ensemble
Monte Carlo histories required by the 5D breakup closure:
`M[4/3,0,0,0,0]` and `M[4/3,1,0,0,0]`. The latter uses the native kK
temperature coordinate.

The copied configuration controls CQMOM. The case uses 10,000 initialization
samples with seed `123`, breakup seed `321`, and end time `0.1`.
