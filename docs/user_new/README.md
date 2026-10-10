# QMOMKit user guide

QMOMKit represents a population by a few raw moments of its number density and
reconstructs from them a quadrature: a small set of weighted nodes. One
coordinate uses QMOM, several use CQMOM.

```text
raw moments → quadrature → state → time step → state → …
```

**QMOMKit does** the numerics: it reconstructs quadratures, generates moments
from distributions, advances a state by moving its nodes and changing its
moments, checks realizability, and writes output.

**You do** everything else: the meaning and units of the coordinates, the
physics (as callbacks), the time loop and step size, what happens when a step
fails, and anything in physical space, such as meshes, transport, and
parallelism.

## Read in this order

1. [Installation](1_installation.md): build, install, and link QMOMKit.
2. [Getting started](2_getting_started.md): three short programs.
3. [Concepts](3_concepts.md): the ideas the rest of the guide uses.
4. [Moments and quadrature](4_moments_and_quadrature.md): set up a case,
   reconstruct its quadrature, check realizability.
5. [Evolving in time](5_evolving_in_time.md): the time step, callbacks, and
   output.
6. [Coupling QMOMKit](6_coupling_qmomkit.md): use QMOMKit inside a simulation
   code.
7. [Reference](7_reference.md): configuration keys, statuses, and
   troubleshooting.
8. [Going further](8_going_further.md): examples, the Monte Carlo companion,
   and the developer documentation.
