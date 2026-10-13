# QMOMKit tests

The test suite: plain C++ executables registered with CTest, grouped into
`unit/`, `component_regression/`, `verification/` and `public_headers/`, plus
the standalone installed-package check in `package_consumer/`. Each
verification area has a README with its references and tolerances.

Build and run it from the repository root in Debug, because the unit tests check
with `assert()`:

```bash
cmake -S . -B build-tests -DCMAKE_BUILD_TYPE=Debug -DQMOMKIT_BUILD_TESTS=ON
cmake --build build-tests --parallel
ctest --test-dir build-tests --output-on-failure
```

`ctest --test-dir build-tests -N` lists every test. Where a test goes, how it is
registered, the CI checks and the evidence a numerical change needs are in
[Testing](../docs/developer/testing.md).
