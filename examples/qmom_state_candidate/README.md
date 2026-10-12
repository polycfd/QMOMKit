# QMOM State Candidate Example

This example shows how an application hands moments it changed itself to a
`QmomState`. A state cannot be edited, because its moments and its quadrature
always agree. Instead, `QmomCase::state_candidate(...)` creates a new state from
the changed moments, with the node count and options of the case, and reports
in a status whether that worked.

The case holds the raw moments $m_0, \dots, m_3$ of weights 0.5 at abscissas 1
and 3 and asks for two nodes. Its `initial_state()` is the state that the
changed moments then replace. Each changed moment set is checked against a
realizability policy with nonnegative support, as for particle sizes, before
the quadrature is reconstructed:

| Changed moments | Outcome |
| --- | --- |
| $\{2, 5, 14, 41\}$ | `success`: weights 0.5 and 1.5 at abscissas 1 and 3 replace the state. |
| $\{0, 0, 0, 0\}$ | `empty`: no number density is left, so there is nothing to evolve. |
| $\{1, -1, 2, -4\}$ | `realizability_rejected`: a negative mean is impossible on $[0, \infty)$. |

The example keeps the previous state when there is no new one. What to do
instead is the application's decision.
[Moments changed by your code](../../docs/user/6_coupling_qmomkit.md#moments-changed-by-your-code)
explains the statuses and how a candidate fits into an evolution step.

Build and run from this example directory:

```bash
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./qmom_state_candidate_example
```

The executable writes no files. It prints the outcome for each changed moment
set to the terminal.
