# Output verification

## Purpose

`MomentWriter` and `QuadratureWriter` write text tables whose column order and
labels are a contract that readers of the files rely on
([Configuration and output](../../../docs/developer/configuration_and_output.md#output-writers)).
These cases write each table and read it back as such a reader would.

## Cases

`verification_output_writer_qmom` and `verification_output_writer_cqmom` load
`qmom/two_node_discrete.qmom` and `cqmom/independent_product.qmom`, write each
table at one time, $t = 0.25$, and check that:

- the header is `# time` followed by the expected labels, and every row has
  one value per label;
- the QMOM moment and quadrature rows, the CQMOM quadrature rows (node index,
  conditioning path, weight, coordinates) and the `CqmomState` moment row,
  labeled in `moment_indices()` order, carry the values handed to the writer;
- requested CQMOM moments keep their requested, non-monotonic order and equal
  the moments of the analytical product rule of `independent_product.qmom`,
  evaluated in the test.

## Tolerances

Values must agree to $10^{-13}$ absolute, indices exactly. The round trip
through the written text is observed to be exact, and the requested CQMOM
moments differ from the analytical ones by about $10^{-15}$, the round-off of
reconstruction; the bound is a safety bound.

## Reproduce

```bash
ctest --test-dir build-tests -R '^verification_output_' --output-on-failure
```

See [Building and running](../../../docs/developer/testing_and_validation.md#building-and-running).
