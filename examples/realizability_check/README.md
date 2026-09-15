# Realizability Check Example

This example checks a raw moment sequence against the QMOMKit realizability
utilities and reports the resulting status. It is a diagnostic of moment
validity only; it does not reconstruct a quadrature or advance a phase-space
evolution state.

The sequence used by [main.cpp](./main.cpp) is
$m_0=1$, $m_1=1$, $m_2=2$, and $m_3=6$. These are the raw, unnormalized
moments through order three of a unit-density exponential distribution on the
nonnegative real line. The check therefore uses
`qmomkit::MomentSupport::nonnegative` and is expected to report
`strictly_realizable`.

The executable prints the selected support, status, and diagnostic message. A
zero exit code means that `result.realizable()` is true; a nonzero exit code
means that the supplied sequence was not accepted as realizable.

Build and run from this example directory:

```bash
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./realizability_check
```
