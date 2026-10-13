# Development philosophy and verification

This document records how QMOMKit is developed and why its correctness should be
judged the way it is. It is written for anyone deciding whether to rely on the
library, and for contributors who want to understand the standard their work is
held to.

## Use of generative AI

The generative AI tools Claude (Anthropic) and Codex (OpenAI) have been used
extensively during the development of QMOMKit, including for code drafting and
refactoring, numerical investigations, test development, documentation, and
repository review.

The use of these tools has reinforced a central principle of the project:
confidence in scientific software should come from evidence rather than from the
origin of an implementation. Neither human-written nor AI-assisted code is
inherently correct.

This is stated plainly because users of scientific software are entitled to know
how it was produced.

## Evidence over provenance

Because AI-assisted code generation can accelerate development without providing
evidence of correctness, QMOMKit has adopted a deliberately rigorous approach to
testing and numerical verification. New functionality is supported, as
appropriate, by unit and regression tests, manufactured or analytical
verification cases, independent numerical comparisons, and physically meaningful
examples.

This philosophy has resulted in a comprehensive test suite and a broad collection
of examples that serve both to demonstrate the public functionality of QMOMKit
and to provide additional evidence that its numerical methods behave as intended.

The evidence itself is what matters, and it is meant to be inspected rather than
taken on trust:

- the [test suite](tests/) and the principles behind it, in the
  [testing guide](docs/developer/testing.md);
- the [verification evidence](docs/developer/testing.md#verification-evidence)
  a numerical change needs: analytical, manufactured, and convergence-based
  verification;
- the [examples](examples/README.md), which apply the library to concrete
  physical problems through its public API;
- the [Monte Carlo companion](tools/monte_carlo/README.md), which provides
  independent finite-sample reference solutions for comparison.

## Documentation as part of the evidence

Documentation is treated as another essential part of this process. Numerical
assumptions, algorithms, interfaces, and intended usage are documented at several
levels, from focused comments in the implementation to dedicated
[developer](docs/developer/README.md) and [user](docs/user/README.md)
documentation. An algorithm whose assumptions are not written down cannot be
checked against them.

## Scientific claims

The same evidentiary standard applies to scientific content. Scientific claims,
references, numerical methods, and verification results must be checked against
primary sources or independently established evidence. AI-generated
explanations, references, or numerical conclusions are not authoritative sources.

## Contributions

QMOMKit will continue to follow this development philosophy irrespective of
whether future contributions are produced manually or with the assistance of
generative AI. The relevant standards are the correctness, traceability,
verification, documentation, and maintainability of the resulting software.

What this means in practice for a contribution, including the responsibilities
that remain with a contributor using AI assistance, is set out in
[CONTRIBUTING.md](CONTRIBUTING.md). AI coding agents working in this repository
should additionally read and follow [AGENTS.md](AGENTS.md).
