# Installed-package consumer

An external CMake project that checks that an installed QMOMKit is found, its
headers compile and its library links, and that a short QMOM and CQMOM workflow
runs through it; `consumer` returns non-zero if any step fails. It uses only
`find_package(QMOMKit CONFIG REQUIRED)`, `qmomkit::qmomkit` and
`<qmomkit/qmomkit.hpp>`. It is not registered with CTest; build and run it as
in [Building and running](../../docs/developer/testing_and_validation.md#building-and-running).
