# Realizability Check Example

This example checks a raw moment sequence against the QMOMKit realizability
utilities and reports the resulting status.

This is a diagnostic/reference example, not part of the first QMOM or CQMOM
workflow.

Requires: QMOMKit.

Build and run from this example directory:

```bash
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./realizability_check
```
