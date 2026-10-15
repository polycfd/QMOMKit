# Realizability check

This example checks whether a raw moment sequence is realizable on a given
support with `check_realizability(...)`, as described in
[Moments and quadrature](../../docs/user/4_moments_and_quadrature.md#check-realizability).
The moments $m_0 = 1$, $m_1 = 1$, $m_2 = 2$ and $m_3 = 6$ are those of the
exponential distribution $n(\xi) = e^{-\xi}$ of a dimensionless coordinate
$\xi \ge 0$, so the check uses the support `MomentSupport::nonnegative`.

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./realizability_check
```

## Output

The example prints the support, the status and the message of the check. It
writes no files. It exits with code 0 when `result.realizable()` is true and 1
otherwise.

## What to expect

The status is `strictly_realizable`, and the exit code is 0.
