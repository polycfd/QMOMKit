# QMOMKit user guide

Many flows carry a population of dispersed elements, such as particles,
droplets, or bubbles, whose properties differ from element to element. A number
density $n$ describes how many elements there are at each position
$\mathbf{x}$ in physical space with each set of properties
$\boldsymbol{\xi}$, such as size, velocity, or temperature. These properties
span the *phase space*. A population balance equation (PBE) governs how $n$
changes [[Marchisio & Fox, 2013](https://www.doi.org/10.1017/CBO9781139016599)]:

```math
\frac{\partial n}{\partial t}
+ \underbrace{\nabla_{\mathbf{x}} \cdot (\mathbf{u}\, n)}_{\text{physical space}}
+ \underbrace{\nabla_{\boldsymbol{\xi}} \cdot (\dot{\boldsymbol{\xi}}\, n)}_{\text{phase space}}
= S .
```

Elements move through physical space with a velocity $\mathbf{u}$ and through
phase space with a velocity $\dot{\boldsymbol{\xi}}$, for example as they grow
or are dragged by the surrounding flow. The source $S$ creates and removes
elements, for example by breakage or aggregation. Elements that enter or leave
across a boundary of the phase space, as in nucleation, form a boundary
condition.

Solving for $n$ itself is expensive, because the phase space adds dimensions
to the problem. Moment methods track a few raw moments of $n$ instead, which
carry its total, mean, and spread. To evaluate the phase-space terms from
moments, QMOMKit represents $n$ by a quadrature: a small set of weighted nodes
that reproduces the moments. With one phase-space coordinate this is QMOM; with
several, CQMOM.

**QMOMKit solves the phase-space part** of the PBE, the terms with
$\dot{\boldsymbol{\xi}}$ and $S$ and the phase-space boundaries, at one
location: given the moments there, it reconstructs the quadrature, moves its
nodes with $\dot{\boldsymbol{\xi}}$, and integrates the moments with the
source:

```text
raw moments → quadrature → state → time step → state → …
```

It also generates moments from distributions, checks whether moments are
realizable, and writes output.

**Your code does everything else:** the physics, as callbacks that return
$\dot{\boldsymbol{\xi}}$, the moment rates of $S$, and those of the exchange
across phase-space boundaries; the meaning and units of the coordinates; the
time loop and step size; what happens when a step fails; and everything in
physical space, such as meshes, transport with $\mathbf{u}$, and parallelism.
[Coupling QMOMKit](6_coupling_qmomkit.md) shows how the two parts fit
together.

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
7. [Reference](7_reference.md): configuration keys, C++ entry points and
   callbacks, statuses and exceptions, troubleshooting, and output formats.
8. [Going further](8_going_further.md): examples, the Monte Carlo companion,
   citing QMOMKit, and the developer documentation.
