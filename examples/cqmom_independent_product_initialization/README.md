# CQMOM Independent Product Initialization

This example loads an independent-product CQMOM setup from configuration and
shows the complete normal initialization lifecycle:

```text
configuration
    -> load_cqmom(...)
    -> CqmomCase
    -> compute_cqmom(...)
    -> quadrature nodes and weights

CqmomCase
    -> CqmomCase::initial_state()
    -> CqmomState
```

`CqmomCase` is the immutable configuration and problem definition.
`CqmomState` is the runtime state that a transport workflow evolves later.
This example creates the state but does not add transport physics.

Requires: QMOMKit.

Required: `configuration.qmom`, `main.cpp`, and an installed QMOMKit package.
Optional: adapt the distribution, `node_counts`, or `conditioning_order` in
`configuration.qmom`.

Build and run from this example directory:

```bash
cmake -S . -B build \
    -DCMAKE_PREFIX_PATH="../../install"
cmake --build build
cd build
./cqmom_independent_product_initialization_example
```

The executable prints the CQMOM hierarchy and reconstructed node weights and
abscissas.
