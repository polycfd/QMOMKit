# QMOMKitMonteCarlo

This directory contains the optional QMOMKitMonteCarlo companion package. It
adds finite ensemble sampling, sample-derived moment estimation, continuous
sample evolution, and moment comparison utilities for QMOMKit validation and
analysis.

The companion package is built separately from the core QMOMKit library and
links against an installed QMOMKit package.

## Build

Install QMOMKit from the repository root first:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build
cmake --install build
```

Then build the companion package:

```bash
cmake -S tools/monte_carlo -B tools/monte_carlo/build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install" \
    -DCMAKE_INSTALL_PREFIX="$PWD/install"

cmake --build tools/monte_carlo/build
cmake --install tools/monte_carlo/build
```

This builds the library and the user-facing sample-based initialization example.
It installs the companion package into the same repository-local `install/`
prefix as QMOMKit. It does not compile the full test suite unless testing is
enabled.

## Run The Example Driver

```bash
./tools/monte_carlo/build/monte_carlo_sample_based_initialization_example \
    tools/monte_carlo/configuration/qmom_distribution_lognormal.qmom \
    tools/monte_carlo/build/sample_based_qmom
```

The driver writes QMOMKit-style `.txt` tables with commented headers and
whitespace-delimited columns:

```text
sample_based_qmom_moments.txt
sample_based_mc_moments.txt
sample_based_mc_statistics.txt
```

For the CQMOM sample case, use
`configuration/cqmom_independent_product_2d.qmom`.

## Compare Histories

```bash
python3 tools/monte_carlo/scripts/compare_moments.py \
    tools/monte_carlo/build/sample_based_qmom/sample_based_mc_moments.txt \
    tools/monte_carlo/build/sample_based_qmom/sample_based_qmom_moments.txt \
    --mc-statistics tools/monte_carlo/build/sample_based_qmom/sample_based_mc_statistics.txt \
    --output-dir tools/monte_carlo/build/sample_based_qmom/comparison
```

The script writes `moment_comparison.txt` with these error columns:

```text
signed_error = computed - reference
absolute_error = abs(signed_error)
```

## Run Tests

Use a separate build when you want the Monte Carlo tests and validation-only
drivers:

```bash
cmake -S tools/monte_carlo -B tools/monte_carlo/build-tests \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$PWD/install" \
    -DBUILD_TESTING=ON

cmake --build tools/monte_carlo/build-tests
ctest --test-dir tools/monte_carlo/build-tests --output-on-failure
```

The installed exported target is:

```cmake
qmomkit_monte_carlo::qmomkit_monte_carlo
```

## Documentation

- User guide: `docs/user/monte_carlo.md`
- Developer guide: `docs/developer/monte_carlo.md`
