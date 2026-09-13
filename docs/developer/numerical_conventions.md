# Numerical Conventions

This page is the stable developer reference for QMOMKit numerical conventions.
It describes the conventions that algorithms, tests, examples, configuration
parsing, and validation workflows must share.

QMOMKit is a numerical C++ library. These conventions do not imply a CFD mesh,
solver loop, boundary-condition model, particle model, physical unit system, or
application-specific coordinate meaning.

Algorithm-specific details belong in:

- [qmom_wheeler.md](qmom_wheeler.md) for Wheeler inversion and adaptive one-dimensional QMOM.
- [cqmom_recursive.md](cqmom_recursive.md) for recursive CQMOM traversal and branch-local data.
- [phase_space_evolution_design.md](phase_space_evolution_design.md) for time integration, node advection, moment sources,
  and phase-space evolution composition.
- [configuration_design.md](configuration_design.md) for parser behavior and configuration-to-request
  mapping.

## Moment Conventions

QMOMKit uses raw moments of the number density function (NDF). A
one-dimensional raw moment is
$$m_k = \int \xi^k n(\xi) \ \mathrm{d}\xi$$
where $n(\xi)$ is the number density function and $m_0$ is the total number
density.

For a discrete quadrature representation, the same moment is
$$m_k = \sum_i w_i x_i^k$$
where $x_i$ is an abscissa, $w_i$ is its quadrature weight, and $k$ is the
moment order. The weights are quadrature contributions to the number density
distribution.

Some initialization helpers start from a normalized distribution $p(\xi)$ and a
specified total number density $n_0$. In that case:

$$n(\xi) = n_0 p(\xi)$$
$$m_k = n_0 E[\xi^k]$$
$$m_0 = n_0$$

This relation is an initialization convenience. The QMOM moments themselves
are moments of the NDF, not probability-density moments with an implicit
density factor.

The value $m_0$ is the total number density. It is not normalized away. A
moment sequence with $m_0 = 1$ is a unit-number-density sequence only because
the caller or generator chose that total number density.

QMOMKit algorithms assume unnormalized raw moments unless a specific API says
otherwise. Do not divide moments by $m_0$ before passing them to QMOMKit unless
the intended mathematical problem is explicitly a unit-number-density problem.

`qmomkit::MomentSet` stores one-dimensional moments in ascending order:

```text
moments[k] == m_k
```

Native QMOM and CQMOM moments are integer-order polynomial moments. Generalized
moments with real exponents may be evaluated from an existing quadrature rule
as derived quantities, but they are not `MomentSet`, `RawMomentTensor`, or
CQMOM inversion inputs.

For an $N$-node one-dimensional QMOM inversion, callers provide at least
$2N$ moments: $m_0, m_1, ..., m_{2N-1}$.

All numerical storage currently uses `double` values in runtime-sized
containers. Algorithms report invalid input or numerical failure through their
documented result/status types; non-empty output data alone is not a success
criterion.

## Tensor Moment Conventions

A multidimensional raw moment uses a phase-space-order multi-index:
$$\alpha = [\alpha_0, \alpha_1, ..., \alpha_{d-1}]$$

The corresponding raw tensor moment is
$$m_\alpha = \sum_j W_j \prod_r \xi_{r,j}^{\alpha_r}$$
where $r$ is a phase-space dimension, $W_j$ is the global quadrature weight of
terminal node $j$, and $\xi_{r,j}$ is that node's coordinate in dimension $r$.

`qmomkit::RawMomentTensor` stores multidimensional raw moments as:
```text
extents[r] = number of stored moment orders in phase-space dimension r
values[flatten_multi_index(alpha, extents)] = m_alpha
```

Tensor storage is row-major in phase-space dimension order, with the last
dimension varying fastest. For

```text
alpha   = [i0, i1, i2]
extents = [n0, n1, n2]
```

the flat index is

```text
flat = ((i0 * n1) + i1) * n2 + i2
```

For `extents = [2, 3]`, enumeration order is:

```text
[0, 0]
[0, 1]
[0, 2]
[1, 0]
[1, 1]
[1, 2]
```

Tensor indices always refer to moment orders in phase-space dimension order.
They are not ordered by conditioning level, callback order, file column order,
or any recursive traversal detail unless a specific internal data structure
explicitly documents another layout.

## Phase-Space Conventions

A phase-space dimension is the stable coordinate axis identifier used by
QMOMKit. Dimension indices are

```text
0, 1, ..., d - 1
```

The same phase-space dimension order is used for:

- raw moment tensor multi-indices;
- tensor extents;
- node counts;
- quadrature abscissa arrays;
- phase-space coordinate vectors passed to velocity callbacks;
- output labels and selected moment labels when those labels refer to tensor
  order.

Coordinate names are metadata. A name such as `size`, `velocity`, `mass`, or
`temperature` can help humans read configuration and output, but it does not
change numerical ordering.

QMOMKit does not assign physical meaning to any coordinate. Velocity
components are ordinary phase-space dimensions. Scalar or internal coordinates
such as particle size, composition, mass, or a reaction coordinate are also
ordinary phase-space dimensions. A workflow that uses both velocity and
internal coordinates must define one explicit phase-space order and use it
consistently in moments, quadrature, callbacks, output, and validation.

For phase-space evolution callbacks, a node coordinate vector is ordered by phase-space
dimension:

```text
coordinates[r] corresponds to phase-space dimension r
```

The returned phase-space velocity vector uses the same order and size. See
`phase_space_evolution_design.md` for phase-space evolution-specific callback behavior.

## Conditioning Order Conventions

Conditioning order is separate from phase-space dimension order. It defines
the order in which CQMOM conditions dimensions:

```text
phase_space_dimension = conditioning_order[level]
```

The default conditioning order is identity order:

```text
[0, 1, ..., d - 1]
```

A non-identity order is a permutation of phase-space dimensions. For example:

```text
phase-space dimensions = [0, 1, 2]
conditioning_order    = [2, 0, 1]
```

means CQMOM processes dimension `2` first, then dimension `0`, then dimension
`1`.

Conditioning order affects the numerical conditioning sequence used by CQMOM.
It does not change:

- tensor storage order;
- the meaning of phase-space dimensions;
- node-count ownership;
- coordinate vector order;
- final quadrature abscissa dimension order.

Node counts are attached to phase-space dimensions:

```text
node_counts[r] == requested node count for phase-space dimension r
```

`PhaseSpaceHierarchy` validates that conditioning order has length `d` and is a
permutation with no duplicate, missing, or out-of-range dimensions.

Recursive path storage, branch-local moment slices, and traversal details are
implementation topics for `cqmom_recursive.md`.

## Quadrature Conventions

QMOMKit uses the terms node and abscissa as follows:

- A quadrature node is one weighted support point in the reconstructed rule.
- An abscissa is the coordinate value of a node in one phase-space dimension.
- A weight is the quadrature contribution associated with a node in the number
  density distribution.

In one dimension, a quadrature rule represents the supplied moments by
$$m_k = \sum_i w_i x_i^k$$
In multiple dimensions, terminal quadrature nodes represent moments by
$$m_\alpha = \sum_j W_j \prod_r \xi_{r,j}^{\alpha_r}$$
where $W_j$ is the global terminal-node weight.

One-dimensional Wheeler output is sorted by ascending abscissa, with weights
permuted consistently:

```text
weights[i] corresponds to abscissas[i]
```

The order of source nodes used to generate a moment sequence is irrelevant;
raw moments are sums over nodes.

For multidimensional CQMOM quadrature, final abscissas are stored by
phase-space dimension. For terminal node `j`:

```text
weights[j]                    = global terminal-node weight
abscissas[dimension][j]       = coordinate in phase-space dimension
```

Fixed rectangular CQMOM quadratures may also support phase-space-order
multi-index access. Adaptive or ragged quadrature representations should be
interpreted by iterating over exposed terminal nodes, weights, abscissas, and
paths rather than assuming rectangular storage.

Local conditional weights that arise inside CQMOM reconstruction are not the
same as final global weights. The distinction between local conditional
weights and accumulated terminal weights belongs in [cqmom_recursive.md](cqmom_recursive.md).

## Distribution And Initialization Conventions

Distribution helpers are initialization tools. They generate raw moments of an
NDF and raw NDF moment tensors; the rest of the QMOMKit workflow consumes those
moments. QMOM and CQMOM reconstruction do not retain a probability-density
object.

One-dimensional analytical distribution helpers interpret the supplied
analytical distribution as a normalized shape $p(xi)$ and the supplied $n_0$ as
the total number density:
$$n(\xi) = n_0 p(\xi)$$
$$m_k = n_0 E[\xi^k]$$
$$m_0 = n_0$$
The value $n_0$ is an absolute total number density. It is not a fraction or a
probability unless the caller intentionally chooses a unit-number-density
setup.

One-dimensional mixtures use absolute component number densities:
$$m_k = \sum_s n_s E_s[\xi^k]$$
$$m_0 = \sum_s n_s$$

Component number densities are not normalized and are not required to sum to
one. Zero-density components contribute zero.

Independent multidimensional product initialization uses normalized
one-dimensional marginals and one global total number density:
$$n(\xi_0, ..., \xi_{d-1}) = n_0 \prod_r p_r(\xi_r)$$
$$\int p_r(\xi_r) \ \mathrm{d}\xi_r = 1$$
$$m_\alpha = n_0 \prod_r E_r[\xi_r^{\alpha_r}]$$
$$m_{[0, ..., 0]} = n_0$$

Marginals are ordered by phase-space dimension. The global number density is
applied once. Per-dimension number densities, conditioning-order-dependent
marginals, and correlated multidimensional distributions are not part of the
current initialization convention.

Requested maximum orders are inclusive. In one dimension, `maximum_order = 3`
generates `m_0` through `m_3`. For product tensors,
`maximum_orders[r] = 3` gives `extents[r] = 4`.

Configuration files preserve the same conventions. Explicit moments are
already raw moments. Distribution-based initialization is converted to raw
moments before QMOM case objects or CQMOM case objects are constructed.
Parser-specific
precedence and field mapping belong in [configuration_design.md](configuration_design.md).

## Rules For Maintainers

Preserve these conventions unless a deliberate convention change is reviewed,
documented, and covered by tests:

- raw moments are moments of the number density function;
- `m_0` is total number density;
- one-dimensional moments are stored in ascending order;
- tensor multi-indices are ordered by phase-space dimension;
- tensor flat storage is row-major with the last phase-space dimension fastest;
- phase-space dimension order is independent of conditioning order;
- node counts are attached to phase-space dimensions;
- coordinate vectors use phase-space dimension order, and velocity callbacks
  receive all nodes packed node-major with that order within each node;
- one-dimensional quadrature abscissas are sorted ascending with corresponding
  weights permuted consistently;
- final multidimensional quadrature weights are global terminal weights;
- analytical mixtures use absolute component number densities;
- independent products use normalized marginals and one global total number
  density;
- no algorithm should infer physical units, supports, velocity meanings, or
  particle-specific semantics from dimension names.