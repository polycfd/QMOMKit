# Testing and validation

This page covers where a test goes, how it is registered, built and run, the
CI checks, and the evidence a numerical change needs. Tests are plain C++
executables registered with CTest; there is no test framework, and `ctest -N`
lists them all. The rules on what every change must be tested with are in
[AGENTS.md](../../AGENTS.md#testing-and-validation).

## The test tree

| Directory | CTest label | Put a test here when it checks |
| --- | --- | --- |
| [`tests/unit/`](../../tests/unit) | `unit` | one function, type, option, status or failure path in isolation |
| [`tests/component_regression/`](../../tests/component_regression) | `component_regression` | an interaction between components, status propagation, or fragile behavior a maintainer could break while verification passes |
| [`tests/verification/`](../../tests/verification) | `verification` | a mathematical property or public workflow against an analytical, manufactured, convergence or otherwise independent reference |
| [`tests/public_headers/`](../../tests/public_headers) | `public_headers` | that each installed header compiles when included first; add an `include_<header>.cpp` with every public header |
| [`tests/package_consumer/`](../../tests/package_consumer) | not in CTest | that an external project finds, includes and links the installed package |

Examples under [`examples/`](../../examples) are user documentation built
against the installed package, never tests
([AGENTS.md](../../AGENTS.md#tests-and-examples)). When a change affects an
example, run it and check its numbers.

Unit tests check with `assert()`. Verification and component-regression tests
check with the `QMOMKIT_VERIFY_*` macros of
[`tests/verification/verification_helpers.hpp`](../../tests/verification/verification_helpers.hpp)
or with explicit return codes; on failure they print the measured values and
bounds to standard error and exit non-zero, in any build type. A verification
case prints nothing when it passes; its source header states the reference,
its construction and the tolerance, and each area under `tests/verification/`
has a README.

`NDEBUG` removes `assert()`, so only a Debug build of the suite is evidence
that the unit tests passed.

## Registering a test case

Every executable registered by the scanning helpers below ends in a call to
`qmomkit::test::dispatch_test_case` from
[`tests/test_case_dispatch.hpp`](../../tests/test_case_dispatch.hpp), which
maps selectors to functions; the other executables parse their own arguments,
a `.qmom` path and a case name, or a case name. Without an argument the executable runs
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
- Do not name helper functions `test_*` in a scanned file, nor zero-argument
  helpers `verify_*` in a file scanned for logical cases; they would be
  registered.
- A unit source file is listed in `unit_test_sources` in
  [`tests/unit/CMakeLists.txt`](../../tests/unit/CMakeLists.txt), and a scanned
  component-regression source in `component_case_sources` in
  [`tests/component_regression/CMakeLists.txt`](../../tests/component_regression/CMakeLists.txt).

Every registered test carries the label of its directory and a timeout of
`QMOMKIT_TEST_TIMEOUT_SECONDS` (300 s), so a case that fails by not terminating
reports a failure. In Debug, tests compile with `-Wall -Wextra`.

## Public headers

A new header in `include/qmomkit/` is included by `qmomkit.hpp`, listed in the
`install(FILES …)` list of the root [`CMakeLists.txt`](../../CMakeLists.txt),
and given an `include_<header>.cpp` in `tests/public_headers/`, listed in that
directory's `CMakeLists.txt`. That file includes the header under test first,
before any standard header, so a missing include fails to compile; where
clang-format (`IncludeBlocks: Regroup`) would sort a standard header ahead of
it, a comment line separates the two. It may pin declarations with
`static_assert`s in namespace `public_header_pins`. An input aggregate also
enters `brace_initialization.cpp`
([Architecture](architecture.md#public-and-internal-boundary)).

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
`qmomkit` target and needs no installation. Without `CMAKE_BUILD_TYPE`, a
single-configuration generator builds `Release`, so always give `Debug` for
tests; without an install prefix, the prefix is `install/` in the source
directory.

```bash
cmake -S . -B build-tests -DCMAKE_BUILD_TYPE=Debug -DQMOMKIT_BUILD_TESTS=ON
cmake --build build-tests -j8
ctest --test-dir build-tests --output-on-failure -j8
```

The instrumented build compiles a different counter implementation and the
counter assertions in the tests; neither exists in the build above (see
[Performance counters](performance_counters.md)):

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

A numerical change needs evidence beyond agreement with the implementation
it replaces, and a validation change states its reference, assumptions, error
metric, tolerance and reproduction
([AGENTS.md](../../AGENTS.md#testing-and-validation)). The references used
here:

- analytical results: exact discrete quadratures, classical Gaussian rules,
  closed-form ODE, source and node trajectories;
- manufactured cases: moments generated from prescribed nodes, weights or
  conditional quadratures, which the method must recover;
- convergence: the observed order under step refinement, compared with the
  formal order only for fixed steps; adaptive steps test the controller;
- independent references: separately constructed data or implementations;
- Monte Carlo references, below.

Expected values come from the mathematical case, never from another call into
the code under test. Evidence that the suite enforces lives in
`tests/verification/`. [`paper/verification/`](../../paper/verification)
holds paper-specific studies outside CTest.

## Monte Carlo references

The companion package under [`tools/monte_carlo/`](../../tools/monte_carlo)
provides an independent finite-sample reference for raw moments, not for nodes
or branch structure. An ensemble of $P$ samples with number density $n_0$
estimates

```math
M_\alpha = \frac{n_0}{P} \sum_{p=1}^{P} \prod_{d} x_{p,d}^{\alpha_d},
```

with standard error $n_0 \sqrt{s^2/P}$, where $s^2$ is the unbiased sample
variance of the monomial (zero for $P < 2$). Why a comparison starts QMOMKit
from the moments of the sample is in the
[companion README](../../tools/monte_carlo/README.md).

A statistical test accepts an estimate within a few standard errors of a
reference that comes from the mathematical case, not from the companion:
[`tools/monte_carlo/tests/weighted_point_evolution.cpp`](../../tools/monte_carlo/tests/weighted_point_evolution.cpp)
samples a few weighted points, moves the samples with `evolve_ensemble(...)`,
and compares the estimated moments within six standard errors with the exact
moments of the points moved along their analytical trajectories.

Conventions that matter when adding a case:

- Each sampling call seeds its `std::mt19937_64` engines deterministically from
  the explicit `seed`. Samples repeat for the same distribution, sample count, seed, compiler and
  standard library, but standard-library distributions differ between
  implementations; tests compare with tolerances or standard errors, never with
  exact sampled values.
- Moment histories are written with QMOMKit's `MomentWriter`
  (`make_moment_history_writer(...)`), so they have the format of QMOMKit's own
  output.
- Code that only tests use, such as a reference construction, lives under
  `tools/monte_carlo/tests/`, not in the package.
- Tests register as `monte_carlo_<subsystem>_<case>` (label `monte_carlo`) from
  `void test_<case>()` definitions, as in the core suite, and check with
  [`tools/monte_carlo/tests/test_checks.hpp`](../../tools/monte_carlo/tests/test_checks.hpp), not
  `assert()`.

Build and usage are in [`tools/monte_carlo/README.md`](../../tools/monte_carlo/README.md).
