# QMOMKit paper studies

`paper/verification/` contains standalone, explicitly run numerical studies for
the QMOMKit paper. They complement, but do not replace, the automated tests
under `tests/`. `paper/examples/` is reserved for future paper examples.

Install QMOMKit, then configure the verification project in an out-of-source
directory:

```bash
cmake -S paper/verification -B <verification-build-dir> \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH=<qmomkit-install-prefix>
cmake --build <verification-build-dir>
```

Run the studies from that build directory so that generated files remain
outside the source tree:

```bash
cd <verification-build-dir>
./qmomkit_paper_quadrature_1d
./qmomkit_paper_cqmom_affine_chain
./qmomkit_paper_adaptive_branch_quadrature
./qmomkit_paper_time_convergence
```

They write `quadrature_1d.csv`, `cqmom_affine_chain.csv`,
`adaptive_branch_quadrature.csv`, `rk_convergence.csv`, and
`splitting_convergence.csv` in the current directory.

The adaptive study loads its two-dimensional `(d,u)` hierarchy, including both
node counts and conditioning order, from `configuration.qmom`. The first node
count determines the number of size branches, and the second is the maximum
conditional velocity-node count for each branch. Configurations such as
`{1,5}`, `{2,3}`, and `{2,2}` require no source changes. Its CSV uses one row
per output time and size branch with `time`, `total_mass`, `mean_velocity`,
`velocity_variance`, `branch_index`, and `active_nodes`. The relaxation,
reinjection pulse, final time, and timestep remain fixed C++ definitions because
they specify the benchmark physics. The configured splitting, node-advection,
and output options are loaded through QMOMKit.

Create the three figures from the same build directory with NumPy and
Matplotlib:

```bash
python3 <source-dir>/paper/verification/adaptive_branch_quadrature/plot.py
python3 <source-dir>/paper/verification/time_convergence/plot.py
```

The scripts write `adaptive_branch_quadrature.pdf`, `rk_convergence.pdf`, and
`splitting_convergence.pdf` in the current directory.
