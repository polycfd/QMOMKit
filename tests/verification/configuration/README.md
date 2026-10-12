# Configuration failure verification

## Purpose

A faulty `.qmom` file must be refused by `load_qmom(...)` or
`load_cqmom_prepared(...)` with an error that names the file and the cause, and
no case may be returned that a caller could go on to use
([Configuration and output](../../../docs/developer/configuration_and_output.md#parser-structure)).

## Cases

`verification_configuration_<fixture>` loads one fixture holding one
deliberate fault:

- `malformed_syntax`\*: the section header `[phase_space` is not closed;
- `qmom_invalid_adaptive`\*: `adaptive = yes`;
- `qmom_missing_moments`: neither `[moments]` nor an `[initial]` distribution
  or mixture;
- `qmom_insufficient_moments`: three moments for two nodes;
- `cqmom_invalid_vandermonde_solver`\*: `solver = structured`;
- `cqmom_missing_phase_space`: no `[phase_space]` section;
- `cqmom_invalid_conditioning_order`\*: `conditioning_order = 0, 0`;
- `cqmom_incomplete_tensor`: $M_{1,1}$ is missing.

Each case checks that loading fails, that no value is returned, and that the
error contains the file path and a fragment of the expected message; checking
failure alone would let an unrelated error pass. For the fixtures marked \*,
whose fault sits on one line, the error must read `<file>:<line>: ...` with
the line of the fault.

## Tolerances

None: every check is exact.

## Reproduce

```bash
ctest --test-dir build-tests -R '^verification_configuration_' --output-on-failure
```

See [Building and running](../../../docs/developer/testing_and_validation.md#building-and-running).
