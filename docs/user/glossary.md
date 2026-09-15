# QMOMKit Glossary

This glossary explains common scientific terms used in the QMOMKit user guide.

## Mathematical Concepts

| Term | Explanation | More information |
|---|---|---|
| distribution | A mathematical description of how population mass is spread over one or more coordinates. In QMOMKit, distributions are mainly used to generate initial raw moments. | [Distributions](distributions.md) |
| generalized moment | A moment of integer or fractional order evaluated from a valid quadrature using finite real exponents. Generalized moments are derived post-processing values, not reconstruction inputs or phase-space-evolved state. | [One-Dimensional QMOM](qmom.md) |
| moment | A weighted integral or sum of powers of coordinates that summarizes a distribution. Moment-method workflows store and evolve these summaries instead of the full distribution. | [One-Dimensional QMOM](qmom.md) |
| phase space | The ordered set of coordinates used to describe the represented population. Coordinate names are labels; the chosen order controls moment indices, quadrature coordinates, callbacks, and output. | [Configuration](configuration.md) |
| raw moment | An uncentered, unnormalized moment of the number density function. QMOMKit preserves the total number density instead of automatically scaling moments to probability moments. | [Distributions](distributions.md) |
| realizability | The question of whether a supplied moment sequence can come from a nonnegative distribution on a chosen support. It matters because impossible moments may fail quadrature reconstruction or need to be rejected during evolution. | [Realizability](realizability.md) |
| zeroth moment | The moment with all orders equal to zero, representing total number density or total quadrature weight. In one dimension this is $m_0$; in multiple dimensions it is $M_{[0,...,0]}$. | [Configuration](configuration.md) |

## Quadrature Concepts

| Term | Explanation | More information |
|---|---|---|
| abscissa | The coordinate value of a quadrature node in one phase-space dimension. One-dimensional QMOM returns abscissas sorted in ascending order. | [One-Dimensional QMOM](qmom.md) |
| quadrature | A compact weighted-node representation of moments. When reconstruction succeeds, the quadrature reproduces the requested raw moments. | [One-Dimensional QMOM](qmom.md) |
| quadrature node | One representative support point in the quadrature. In CQMOM, a terminal node has one coordinate value in each phase-space dimension. | [Multidimensional CQMOM](cqmom.md) |
| quadrature reconstruction | The numerical step that converts raw moments into quadrature nodes and weights. Always check whether reconstruction was accepted before using the quadrature. | [One-Dimensional QMOM](qmom.md) |
| quadrature weight | The amount of number density carried by a quadrature node. The sum of weights reproduces the zeroth moment for a successful reconstruction. | [One-Dimensional QMOM](qmom.md) |
| terminal node | In CQMOM, a quadrature node carrying one coordinate in every phase-space dimension, reached by conditioning through each dimension in turn. The terminal node count is the product of the per-dimension node counts, so it grows quickly with phase-space dimension. | [Multidimensional CQMOM](cqmom.md) |
| moment reconstruction | The post-reconstruction step that evaluates native integer-order raw moments from an accepted quadrature, for output, validation, or handoff to another workflow. | [One-Dimensional QMOM](qmom.md) |

## QMOM Concepts

| Term | Explanation | More information |
|---|---|---|
| adaptive branch counts | A CQMOM option that allows individual branches to use fewer active nodes when the local moments only support a lower-order rule. Adaptive results should be read as terminal nodes rather than assumed to be rectangular. | [Multidimensional CQMOM](cqmom.md) |
| conditional moment | A moment used by CQMOM to reconstruct later dimensions after earlier conditioning nodes have been selected. It lets a multidimensional reconstruction proceed one conditioned dimension at a time. | [Multidimensional CQMOM](cqmom.md) |
| conditioning order | The numerical order in which CQMOM processes phase-space dimensions. It does not rename coordinates or change tensor, callback, or output order. | [Multidimensional CQMOM](cqmom.md) |
| CQMOM | The Conditional Quadrature Method of Moments for multidimensional phase spaces. Use it when the represented population depends on several ordered coordinates. | [Multidimensional CQMOM](cqmom.md) |
| QMOM | The Quadrature Method of Moments for one-dimensional moment sequences. Use it when a population is represented by one phase-space coordinate. | [One-Dimensional QMOM](qmom.md) |

## Simulation Concepts

| Term | Explanation | More information |
|---|---|---|
| boundary contribution | An additive raw-moment rate from phase-space boundary exchange such as injection, removal, or similar production and loss terms. It is not a mesh boundary-condition framework. | [Phase-space evolution](phase_space_evolution.md) |
| moment source | A rate of change applied directly to raw moments. Use it for source-like effects that are better represented as `dm/dt` than as smooth node motion. | [Phase-space evolution](phase_space_evolution.md) |
| node advection | The phase-space evolution operation that moves quadrature abscissas through a velocity field. It changes where nodes sit in phase space and leaves their weights unchanged. | [Phase-space evolution](phase_space_evolution.md) |
| operator splitting | The scheme by which one phase-space evolution step composes node advection with the moment-source update. The default, Strang splitting, applies an advection half-step, a full source stage, and a second advection half-step. | [Phase-space evolution](phase_space_evolution.md) |
| Monte Carlo reference | A stochastic moment history computed from a sample ensemble and used for comparison with deterministic QMOMKit results. It is a reference with statistical uncertainty, not exact ground truth. | [Monte Carlo Companion Package](monte_carlo.md) |
| physics callback | Application-supplied code that gives QMOMKit a velocity, source, or boundary-rate model. QMOMKit advances moments and quadrature, while the application defines the physical meaning. | [Callbacks](callbacks.md) |
| sample ensemble | A finite collection of samples drawn from a distribution or supplied externally. Sample-derived moments contain sampling fluctuation, even with a fixed random seed. | [Monte Carlo Companion Package](monte_carlo.md) |
| SSP-RK3 | The strong-stability-preserving third-order Runge-Kutta scheme. It advances moment-source updates and is one of the selectable node-advection schemes. Each stage can be checked against an acceptance predicate, so a step that would leave a state unusable is retried with a smaller timestep. | [Phase-space evolution](phase_space_evolution.md) |
| phase-space evolution state | Runtime data advanced by a phase-space evolution workflow. QMOM phase-space evolution advances a `QmomState` created from a `QmomCase`; normal CQMOM phase-space evolution advances a `CqmomState` created from a `CqmomCase`. | [Phase-space evolution](phase_space_evolution.md) |
