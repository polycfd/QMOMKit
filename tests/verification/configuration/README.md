# Configuration Failure Verification

Configuration files are user input. Invalid `.qmom` files should fail early
through the public loading APIs before a partial `QmomCase` or `CqmomPreparedCase`
can be used by numerical algorithms.

`configuration_failures.cpp` verifies representative public failure paths for:

- malformed syntax;
- invalid user option values;
- missing required sections or entries;
- mathematically inconsistent case setup.

Each fixture is a small real `.qmom` file. The verification executable loads
the files with `load_qmom(...)` or `load_cqmom_prepared(...)`, checks that loading
fails, checks that no case value is returned, and verifies that structured
errors include useful codes, filename/line information, and human-readable
message fragments. The goal is the user-facing configuration contract, not the
parser's private implementation details.
