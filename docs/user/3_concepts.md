# Concepts

This page explains the ideas and terms the other pages use.

## Raw moments

QMOMKit works with raw moments of the number density $n$ of a population. In
one coordinate and in several, they are

```math
m_k = \int \xi^k\, n(\xi)\,\mathrm{d}\xi ,
\qquad
M_{k_0, \ldots, k_{D-1}} = \int \xi_0^{k_0} \cdots \xi_{D-1}^{k_{D-1}}\, n(\boldsymbol{\xi})\,\mathrm{d}\boldsymbol{\xi} .
```

where $\xi_d$ is coordinate $d$ of $D$. The zeroth moment, $m_0$ or
$M_{0, \ldots, 0}$, is the total number density. QMOMKit never normalizes it
away.

## Quadrature

A quadrature replaces $n$ by $N$ nodes; node $i$ has a weight $w_i$ and an
abscissa $\xi_i$, chosen so that $m_k = \sum_i w_i \xi_i^k$. $N$ nodes need the
moments $m_0, \ldots, m_{2N-1}$. The weights add up to $m_0$.

## Phase space

The phase space is the set of coordinates $\boldsymbol{\xi}$ that describe
each element of the population, such as size, velocity, or temperature.
Coordinate names are labels: QMOMKit gives them no meaning and no units. The
order in which you declare the coordinates is used everywhere, for moment
indices, node coordinates, callbacks, and output.

## QMOM and CQMOM

With one coordinate, QMOM [[McGraw, 1997](https://www.doi.org/10.1080/02786829708965471)]
reconstructs the quadrature with the Wheeler algorithm
[[Wheeler, 1974](https://www.doi.org/10.1216/RMJ-1974-4-2-287)]. With several,
CQMOM [[Yuan & Fox, 2011](https://www.doi.org/10.1016/j.jcp.2011.07.020)]
reconstructs one coordinate after another, each conditioned on the previous
ones; QMOMKit does this recursively
[[Bouguezzoul, 2026](https://publications.polymtl.ca/76505/)]. The
*conditioning order* chooses which coordinate comes first. It changes the
reconstruction, never the order of the coordinates. Each resulting *terminal
node* has one weight and one value per coordinate. CQMOM needs only some of the
moments; QMOMKit derives which ones from the node counts and the conditioning
order.

## Case and state

A *case* is the fixed setup: the initial moments, node counts, and options.
`QmomCase` is the QMOM case. For CQMOM, `CqmomPreparedCase` prepares the
branch structure once and shares it with everything made from it.
`compute_qmom(...)` and `compute_cqmom(...)` reconstruct the quadrature of a
case.

A *state* is what evolves: moments together with the quadrature that matches
them. `initial_state()` creates one from a case, and each time step returns a
new one. A state cannot be edited, so moments that your application changes
become a state through `state_candidate(...)`, as described in
[Coupling QMOMKit](6_coupling_qmomkit.md).

## Adaptive and fixed node counts

Reconstruction is adaptive by default: node counts are maxima, and a
quadrature uses fewer nodes where the moments do not support more. The branches
of a CQMOM quadrature can then differ in size, so iterate over `weights()` and
`abscissas()` rather than indexing nodes on a grid. A fixed reconstruction,
`[qmom] adaptive = false`, returns exactly the requested nodes or fails.

## Results and statuses

Numerical operations report failure in their result rather than throwing. Every
result has a `status` and a `message`, and all but a realizability check have
`ok()`; check it before using the result. A failed time step returns the state
it started from, and what to do next, such as retrying with a smaller step, is
your decision. Loading a configuration is checked separately:
`if (!configured)`, then `configured.error()`. A few functions throw instead,
notably `initial_state()`, when the moments of a case cannot be reconstructed,
and output writers, when a file cannot be written. The
[Reference](7_reference.md) lists all statuses.

## Terms

| Term | Meaning |
| --- | --- |
| abscissa | The position of a node in phase space. |
| boundary contribution | A moment rate from exchange across a phase-space boundary, such as nucleation; added to the source. |
| bounds | An interval the nodes of a coordinate must lie in, `[phase_space] support_lower` and `support_upper`. |
| branch | In CQMOM, the one-dimensional reconstruction of one coordinate for given values of the coordinates earlier in the conditioning order. |
| branch structure | The branches a CQMOM hierarchy can have and the moments each reads, prepared once by `CqmomPreparedCase`. |
| callback | A function you supply with the physics: a velocity, a source, or a boundary contribution. |
| case | The fixed setup of a reconstruction: moments, node counts, options. |
| conditioning order | The order in which CQMOM reconstructs the coordinates. |
| hierarchy | The node counts of a CQMOM case and its conditioning order, `PhaseSpaceHierarchy`. |
| layout | The moments a CQMOM hierarchy requires, in a fixed order, `CqmomMomentLayout`. |
| node | A weighted point of a quadrature. |
| node advection | Moving the nodes with the velocity over a step. |
| number density | Elements per unit volume of physical space and of phase space; its zeroth moment is the total number density. |
| realizability | Whether moments can come from a non-negative number density on a given support. |
| post-source reconstruction | The reconstruction of the quadrature from the moments a source step produced. |
| source | A rate of change of the moments, such as breakage or aggregation. |
| source step | Integrating the moments with the source over a step. |
| state | Moments and their matching quadrature, advanced in time. |
| support | The range a coordinate can take, such as the non-negative numbers, against which a realizability check tests the moments. |
| terminal node | A node of a CQMOM quadrature, with one weight and one value per coordinate. |
| velocity | The rate $\dot{\boldsymbol{\xi}}$ at which nodes move through phase space, such as growth or drag. |
| weight | The share of the number density carried by a node. |
