# QMOMKit paper studies

`paper/verification/` contains explicitly invoked, reproducible numerical studies for the QMOMKit paper. `paper/examples/` is reserved for later application-oriented paper examples. Neither group is part of CTest or replaces the fast automated coverage under `tests/`.

The verification project is a standalone downstream consumer. Install QMOMKit, or otherwise make its package configuration discoverable, then configure it in any out-of-source directory:

```bash
cmake -S paper/verification -B <verification-build-dir> \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH=<qmomkit-install-prefix>
cmake --build <verification-build-dir>
```

Run the executables from that build directory so generated files remain outside the source tree:

```bash
cd <verification-build-dir>
./qmomkit_paper_quadrature_1d
./qmomkit_paper_cqmom_affine_chain
./qmomkit_paper_adaptive_branch_quadrature
./qmomkit_paper_time_convergence
```

The studies write these CSV files in the current working directory:

- `quadrature_1d.csv`
- `cqmom_affine_chain.csv`
- `adaptive_branch_quadrature.csv`
- `rk_convergence.csv`
- `splitting_convergence.csv`

The adaptive study uses the established two-dimensional `(d,u)` variance-collapse and finite-pulse reinjection physics with adaptive CQMOM node counts `{2,3}` and conditioning order `{0,1}`. Branch 0 and branch 1 in its CSV are the stable first-level size paths `[0]` and `[1]`.

The one-dimensional study reconstructs N=1 through 10 reference rules. Its conventions are unit-weight Gauss--Legendre on `[-1,1]`, standard `alpha=0` Gauss--Laguerre with weight `exp(-x)`, and the normalized standard-normal (probabilists' Hermite) weight `(2*pi)^(-1/2) exp(-x^2/2)`. The latter is not the physicists' `exp(-x^2)` convention.

The affine-chain study uses three nodes at every conditional level for dimensions 2, 3, 5, 7, and 10. Its 10D, 59,049-leaf run is an explicit offline paper study, not a unit test.

The temporal study reports fixed-step orders three for SSP-RK3 and RK23, five for RK45, one for `advsrc` and `srcadv`, and two for `strang`.

NumPy and Matplotlib are the only plotting dependencies:

```bash
python3 <source-dir>/paper/verification/adaptive_branch_quadrature/plot.py
python3 <source-dir>/paper/verification/time_convergence/plot.py
```

These create `adaptive_branch_quadrature.pdf`, `rk_convergence.pdf`, and `splitting_convergence.pdf` in the current working directory by default. Both scripts accept optional input and output paths. Generated data and figures are not source files and should be regenerated from a clean checkout.
