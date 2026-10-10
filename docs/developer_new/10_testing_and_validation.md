# Testing and validation

Tests are plain C++ executables registered with CTest; there is no test
framework. The registered tests are the inventory: list them with `ctest`
rather than maintaining a list by hand. The rules on what every change must be
tested with are in [AGENTS.md](../../AGENTS.md#testing-and-validation).

## The test tree

| Directory | CTest label | Put a test here when it checks |
| --- | --- | --- |
| [`tests/unit/`](../../tests/unit) | `unit` | one function, type, option, status or failure path in isolation |
| [`tests/component_regression/`](../../tests/component_regression) | `component_regression` | an interaction between components, status propagation, or fragile behavior a maintainer could break while verification passes |
| [`tests/verification/`](../../tests/verification) | `verification` | a mathematical property or public workflow against an analytical, manufactured, convergence or otherwise independent reference |
| [`tests/public_headers/`](../../tests/public_headers) | `public_headers` | that each installed header compiles when included first; add an `include_<header>.cpp` with every public header |
| [`tests/package_consumer/`](../../tests/package_consumer) | not in CTest | that an external project finds, includes and links the installed package |
| [`tests/benchmarks/`](../../tests/benchmarks) | not in CTest | performance only; see [Benchmarks and performance](11_benchmarks_and_performance.md) |

Examples are never tests and are never registered with CTest.

Unit tests check with `assert()`. Verification and component-regression tests
check with the `QMOMKIT_VERIFY_*` macros of
[`tests/verification/verification_helpers.hpp`](../../tests/verification/verification_helpers.hpp)
or with explicit return codes; on failure they print the measured values and
bounds to standard error and exit non-zero, in any build type. A verification
case prints nothing when it passes; its source header states the reference,
its construction and the tolerance, and each area under `tests/verification/`
has a README.

## Debug rule

`NDEBUG` removes `assert()`, so a Release build compiles the unit checks out.
The suite is authoritative only in a Debug build; a Release CTest result is not
evidence that the unit suite passed.

## Registering a test case

Every test executable ends in a call to `qmomkit::test::dispatch_test_case`
from [`tests/test_case_dispatch.hpp`](../../tests/test_case_dispatch.hpp),
which maps selectors to functions. Without an argument the executable runs
every case; with one selector it runs that case; an unknown selector prints the
available ones and returns 2, which CTest reports as a failure.

[`tests/CMakeLists.txt`](../../tests/CMakeLists.txt) registers cases in three
ways:

| Helper | Used by | Registers |
| --- | --- | --- |
| `qmomkit_register_source_test_cases` | `tests/unit/`, most of `tests/component_regression/` | every definition `void test_<case>(` found in the source |
| `qmomkit_register_source_logical_cases` | several verification sources | every `void test_<case>()` or `void verify_<case>()` |
| `qmomkit_register_test_case` or `add_test` | the other executables, typically driven by a `.qmom` file or case name | the cases listed in the directory's `CMakeLists.txt` |

The scanned forms read the source with a regular expression and name each test
`<label>_<subsystem>_<case>`, where `<subsystem>` is the target name without a
leading `test_` or `component_regression_` (for a unit test, the file stem).
The selector is `<case>`. So:

- In a scanned file, adding a case takes a `void test_<case>()` function and a
  `{"<case>", test_<case>}` entry in the dispatch table; CMake needs no change.
  A scanned function missing from the table fails with exit code 2.
- Do not name helper functions `test_*` in a scanned file; they would be
  registered.
- A unit source file is listed in `unit_test_sources` in
  [`tests/unit/CMakeLists.txt`](../../tests/unit/CMakeLists.txt).

Every registered test carries the label of its directory and a timeout of
`QMOMKIT_TEST_TIMEOUT_SECONDS` (300 s), so a case that fails by not terminating
reports a failure. In Debug, tests compile with `-Wall -Wextra`.

## Listing and selecting tests

```bash
ctest --test-dir build-tests -N                       # list every test
ctest --test-dir build-tests -N -L verification       # list one label
ctest --test-dir build-tests -L unit --output-on-failure
ctest --test-dir build-tests -R '^unit_wheeler_' --output-on-failure
ctest --test-dir build-tests -R '^verification_qmom_analytical_gaussian_n3$' --output-on-failure
./build-tests/tests/unit/test_wheeler nonfinite_moment_reports_invalid_input
```

The labels are `unit`, `component_regression`, `verification` and
`public_headers`; each test has exactly one.

## Building and running

These commands are those of the CI workflows. The suite links the build-tree
`qmomkit` target and needs no installation:

```bash
cmake -S . -B build-tests -DCMAKE_BUILD_TYPE=Debug -DQMOMKIT_BUILD_TESTS=ON
cmake --build build-tests -j8
ctest --test-dir build-tests --output-on-failure -j8
```

The instrumented build compiles a different counter implementation and the
counter assertions in the tests; neither exists in the build above (see
[Benchmarks and performance](11_benchmarks_and_performance.md#performance-counters)):

```bash
cmake -S . -B build-tests-instrumented -DCMAKE_BUILD_TYPE=Debug -DQMOMKIT_BUILD_TESTS=ON \
  -DQMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION=ON
cmake --build build-tests-instrumented -j8
ctest --test-dir build-tests-instrumented --output-on-failure -j8
```

The package consumer and the Monte Carlo companion are separate CMake projects
built against an installed QMOMKit:

```bash
cmake -S . -B build-install -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX="$PWD/install-consumers"
cmake --build build-install -j8 --target install

cmake -S tests/package_consumer -B build-package-consumer -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_PREFIX_PATH="$PWD/install-consumers"
cmake --build build-package-consumer -j8
./build-package-consumer/consumer

cmake -S tools/monte_carlo -B build-monte-carlo -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON \
  -DCMAKE_PREFIX_PATH="$PWD/install-consumers"
cmake --build build-monte-carlo -j8
ctest --test-dir build-monte-carlo --output-on-failure -j8
```

The consumer uses only `find_package(QMOMKit CONFIG REQUIRED)` and
`qmomkit::qmomkit`, includes `<qmomkit/qmomkit.hpp>`, and runs a short QMOM and
CQMOM workflow; it returns non-zero if any step fails. Run it when changing
installation, exports, package metadata or installed headers.

## CI gates

All workflows run on a self-hosted runner, on every push to `main` and on pull
requests.

| Workflow | File | Runs | Reproduce locally |
| --- | --- | --- | --- |
| CI-Debug | [`debug.yml`](../../.github/workflows/debug.yml) | the suite in Debug, then again with instrumentation | the first two blocks of [Building and running](#building-and-running) |
| CI-Release | [`release.yml`](../../.github/workflows/release.yml) | the suite in Release | as CI-Debug with `-DCMAKE_BUILD_TYPE=Release` |
| CI-Warnings-GCC | [`warnings_gcc.yml`](../../.github/workflows/warnings_gcc.yml) | a Debug build of library and tests with the runner's `g++` and `-Werror`; nothing is run | `cmake -S . -B build-warnings-gcc -DCMAKE_BUILD_TYPE=Debug -DQMOMKIT_BUILD_TESTS=ON -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_FLAGS="-Werror"`, then build |
| CI-Indentation | [`clang-format.yml`](../../.github/workflows/clang-format.yml) | [`tools/check_indentation.sh`](../../tools/check_indentation.sh) with clang-format 22.1 | see below |
| CI-Consumers | [`consumers.yml`](../../.github/workflows/consumers.yml) | install, package consumer, Monte Carlo tests, all Debug | the last block of [Building and running](#building-and-running) |

On pull requests, every workflow skips changes confined to `docs/**`; all but
CI-Indentation also skip `examples/**`, and CI-Consumers skips `paper/**`.

The library compiles with `-Wall -Wextra -Wpedantic`. Warnings depend on the
GCC version, so a build that is clean with one `g++` can fail on the runner.

Formatting follows [`.clang-format`](../../.clang-format) and requires
clang-format 22.1 exactly; the scripts check the version.
`./tools/download_clang_format` installs it into `tools/programs/`.
[`tools/indent-all`](../../tools/indent-all) formats every `.cpp` and `.hpp`
under `include`, `src`, `tests`, `examples` and `tools`, and also fixes
trailing whitespace, line endings and final newlines; `check_indentation.sh`
runs it and fails if any file changed. Run from the repository root:

```bash
./tools/download_clang_format          # once
REPORT_ONLY=true ./tools/indent-all    # list files that would change
./tools/indent-all                     # rewrite them
```

[`.clang-tidy`](../../.clang-tidy) defines the static-analysis checks and
naming rules. No workflow runs clang-tidy; code must nevertheless be
compatible with it.

## Validation evidence

A numerical change needs evidence beyond agreement with the
implementation it replaces. What counts:

- analytical results: exact discrete quadratures, classical Gaussian rules,
  closed-form ODE, source and node trajectories;
- manufactured cases: moments generated from prescribed nodes, weights or
  conditional quadratures, which the method must recover;
- convergence: observed order under step refinement; formal order only with
  fixed steps, while adaptive cases measure controller behavior;
- independent references: separately constructed data or implementations;
- Monte Carlo references, below.

The expected values must come from the mathematical case, not from another call
into the production code under test. Benchmark timings are never evidence of
correctness. A validation change states how the reference is constructed, the
assumptions, the error metric, why the tolerance is appropriate, and how to
reproduce the result
([AGENTS.md](../../AGENTS.md#testing-and-validation)); changing a tolerance,
fallback, node ordering or status mapping needs explicit approval.

Evidence that the suite enforces lives in `tests/verification/`.
[`paper/verification/`](../../paper/verification) holds paper-specific studies;
it is not part of CTest or of this suite.

## Monte Carlo references

The companion package under [`tools/monte_carlo/`](../../tools/monte_carlo)
provides an independent finite-sample reference. It compares raw moments, not
nodes or branch structure. An ensemble of $N$ samples with number density $n_0$
estimates

```math
M_\alpha = \frac{n_0}{N} \sum_{p=1}^{N} \prod_{d} \xi_{p,d}^{\alpha_d},
```

with standard error $n_0 \sqrt{s^2/N}$, where $s^2$ is the unbiased sample
variance of the monomial (zero for $N < 2$).

`sample_initialization(...)` returns an ensemble together with the moments
estimated from that same ensemble. Starting QMOMKit from those moments while
the ensemble is evolved by `evolve_ensemble(...)` removes the initial sampling
difference, so what remains comes from the quadrature representation and its
evolution. Starting QMOMKit from analytical moments instead leaves the sampling
error of the initial ensemble in the comparison. The reference itself carries
sampling noise that shrinks as $N^{-1/2}$.

`compare_moments(reference, computed, standard_errors)` in
[`src/comparison.cpp`](../../tools/monte_carlo/src/comparison.cpp) reports, per
moment:

| Metric | Definition |
| --- | --- |
| `signed_error` | computed − reference |
| `absolute_error` | $\lvert$`signed_error`$\rvert$ |
| `relative_error` | `absolute_error` / max($\lvert$reference$\rvert$, standard error, $10^{-300}$) |
| `z_score` | `signed_error` / standard error; 0 for exact agreement and ±∞ otherwise when the standard error is zero |

The script `scripts/compare_moments.py` uses $\lvert$reference$\rvert + \epsilon$
as the relative-error denominator ($\epsilon = 10^{-30}$, `--epsilon`). A large
z-score is a diagnostic, not a failure verdict.

Conventions that matter when adding a case:

- Each sampling call seeds one `std::mt19937_64` with the explicit `seed`.
  Samples repeat for the same distribution, sample count, seed, compiler and
  standard library, but standard-library distributions differ between
  implementations; tests compare with tolerances or standard errors, never with
  exact sampled values.
- Output tables are `.txt` with a `#` header and 17 significant digits, and
  record seed, sample count and number density.
- Tests register as `monte_carlo_<subsystem>_<case>` (label `monte_carlo`) from
  `void test_<case>()` definitions, as in the core suite, and check with
  [`tests/test_checks.hpp`](../../tools/monte_carlo/tests/test_checks.hpp), not
  `assert()`.

Build and usage are in [`tools/monte_carlo/README.md`](../../tools/monte_carlo/README.md).
