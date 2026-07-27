# Installed Package Consumer Smoke Test

This directory is a minimal external CMake project used to validate the
installed QMOMKit package interface. It is intentionally separate from the
normal CTest unit, verification, and public-header suites because it requires
QMOMKit to be built and installed before the consumer project is configured.

The consumer uses only:

```cmake
find_package(QMOMKit CONFIG REQUIRED)
target_link_libraries(consumer PRIVATE qmomkit::qmomkit)
```

It must not use repository source paths, `add_subdirectory`, internal headers,
or test utilities.

## Run

From the top-level QMOMKit checkout:

```bash
cmake -S . -B build
cmake --build build
cmake --install build

cmake -S tests/package_consumer \
  -B tests/package_consumer/build \
  -DCMAKE_PREFIX_PATH="install"

cmake --build tests/package_consumer/build
```

This verifies that `find_package(QMOMKit CONFIG REQUIRED)` succeeds, the
installed `qmomkit::qmomkit` target resolves, and installed public headers
compile through `<qmomkit/qmomkit.hpp>`.
