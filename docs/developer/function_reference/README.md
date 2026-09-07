# Implementation and function reference

This reference connects mathematical operations to their implementation, control
flow, and numerical contracts. It explains the current source rather than
repeating API signatures. User documentation explains supported workflows;
public headers define interfaces; the developer algorithm pages explain the
mathematical design. These pages connect those layers to the functions and
callers that implement them.

When this prose and the source diverge, the source is authoritative. The pages
describe implemented behavior and its important limitations; they do not make
universal numerical guarantees from bounded test cases.

## Functional index

| Reference | Responsibility |
|---|---|
| [QMOM reconstruction and scaling](qmom_reconstruction.md) | Raw sequence → recurrence → eigenproblem → accepted node count; bounds and coordinate scaling |
| [Realizability](realizability.md) | Support-specific matrix classification, range checks, acceptance predicates |
| [Conditional systems and linear algebra](conditional_systems.md) | Weighted Vandermonde equations, BP, dense elimination, QR, SVD, eigenkernel |
| [Recursive CQMOM](cqmom_reconstruction.md) | Required layouts, conditional slices, branch traversal, active-count constraints |
| [State and representations](state_and_representations.md) | Phase-space indices, conditioning paths, dense/lean storage, accepted state and metadata ownership |
| [Moment evaluation and rates](moment_evaluation.md) | Integer/generalized quadrature moments and chain-rule advection rates |
| [Sources and source contexts](moment_sources.md) | Stage-local moment access, boundary addition, source interval completion, marginal predicates |
| [Time integration and node motion](time_integration.md) | SSP-RK3, embedded RK23/RK45, controllers, node packing |
| [Phase-space evolution](phase_space_evolution.md) | Operator clocks, transactional state transitions, topology changes, batch calls |
| [Initialization and configuration](initialization_and_configuration.md) | Distribution formulas, mixtures/products, parser decisions and policy ownership |
| [Output and support utilities](utilities.md) | Writer contracts, finite-value validation, instrumentation |
| [Coverage inventory](coverage_inventory.md) | Definition-to-subsystem routing and overload locations |

## How to use this reference

**Level A** sections explain numerical transformations, important loops, branch
conditions, failures, and sensitivities. Related small functions sometimes share
one section when they implement a single mathematical step. **Level B** entries
cover orchestration and policy handoffs. **Level C** tables group elementary
access, formatting, validation, and forwarding functions by actual responsibility.
The inventory is navigation and coverage bookkeeping, not a substitute for these
explanations. Anonymous namespaces are identified as local scope; unqualified
public symbols are in `qmomkit`, and shared private symbols in `qmomkit::internal`.

Source and test links provide traceability for the stated contracts. A named test
shows the behavior exercised by that test; it is not a numerical certification.

“The implementation establishes” denotes observable control flow or arithmetic.
An inferred purpose is explicitly marked. “No direct test identified” means the
reference does not identify a targeted assertion; it does not exclude indirect
coverage.

Accepted reconstruction, mathematical realizability, configured capacity, active
counts, and integration acceptance are distinct concepts. Consult [numerical
conventions](../numerical_conventions.md), [testing](../testing.md), and
[validation](../validation.md) alongside it.
