# QMOM states from changed moments

This example shows how an application hands moments it changed itself to a
`QmomState`. A state cannot be edited, because its moments and its quadrature
always agree. Instead, `QmomCase::state_candidate(...)` creates a new state from
the changed moments, with the node count and options of the case, and reports
in a status whether that worked; see
[Moments changed by your code](../../docs/user/6_coupling_qmomkit.md#moments-changed-by-your-code).

## Model

The case holds the raw moments $m_0, \dots, m_3$ of weights 0.5 at abscissas 1
and 3 and asks for two nodes; its `initial_state()` is the state that the
changed moments replace. Each changed moment set is checked against a
realizability policy with nonnegative support, as for particle sizes, before
the quadrature is reconstructed. Where there is no new state, the example
keeps the previous one; what to do instead is the application's decision.

| Changed moments | Outcome |
| --- | --- |
| $\{2, 5, 14, 41\}$ | `success`: weights 0.5 and 1.5 at abscissas 1 and 3 replace the state. |
| $\{0, 0, 0, 0\}$ | `empty`: no number density is left, so there is nothing to evolve. |
| $\{1, -1, 2, -4\}$ | `realizability_rejected`: a negative mean is impossible on $[0, \infty)$. |

## Build and run

Install QMOMKit first, as described in
[Installation](../../docs/user/1_installation.md). From this directory:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=../../install
cmake --build build
cd build
./qmom_state_candidate_example
```

## Output

The example prints the outcome for each changed moment set: the new nodes, or
why the state is kept. It writes no files.

## What to expect

One outcome per set, as in the table; the rejected set reports that its
moments are not realizable on nonnegative support.
