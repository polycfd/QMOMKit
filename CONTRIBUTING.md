# Contributing to QMOMKit

Contributions to QMOMKit are welcome. 

QMOMKit aims to remain a focused scientific computing library with clear interfaces, 
minimal dependencies, and carefully tested numerical behaviour. Any changes should 
therefore preserve the project's focus on correctness, clarity, numerical robustness, 
and maintainable implementation.

## Before making changes

Please:

- read the relevant documentation and existing implementation before modifying numerical behaviour;
- keep changes focused on the problem being addressed;
- prefer existing QMOMKit functionality and configuration options over application-specific reimplementations;
- avoid unrelated cleanup in the same contribution;
- add or update tests and verification cases where they provide meaningful evidence for the change; and
- update the relevant user and developer documentation whenever behaviour, interfaces, configuration, or numerical methods change.

## Numerical and scientific changes

Changes to numerical algorithms should be supported by appropriate evidence.
Depending on the change, this may include:

- unit tests for well-defined local behaviour;
- regression tests protecting existing behaviour;
- manufactured or analytical verification cases;
- independent numerical comparisons;
- convergence studies;
- comparison with established literature results.

Passing the existing test suite is necessary but is not, by itself, evidence
that a new numerical method or model is correct.

Scientific claims, references, and validation data should be checked against
primary sources wherever practical.

## New features

New QMOMKit features must be accompanied by evidence that they work as
intended and by an example showing how they are used in a physically meaningful
setting.

A new feature should therefore include:

- suitable unit, regression, and/or verification tests, chosen according to
  the numerical and software behaviour being introduced;
- a physical example demonstrating the feature through the public QMOMKit
  interface and configuration system; and
- appropriate documentation for both users and developers.

Tests should establish the relevant numerical and software contracts rather
than merely exercise code paths. In particular, new numerical functionality
should normally include verification or comparison against independently known
behaviour where such evidence can reasonably be constructed.

The physical demonstration may be a new example or a clearly motivated
variant or extension of an existing example. It should exercise the feature
through the same public API and configuration options that users are
expected to employ. Example code should not duplicate, bypass, override, or
reimplement behaviour that is already provided by QMOMKit.

A feature should not be considered complete until its implementation, tests,
example, and documentation are mutually consistent.

## Documentation

All contributions should be documented at the level appropriate to their
scope.

User-facing changes should update the relevant [user documentation](docs/user/README.md), 
including configuration options, public API behaviour, assumptions, units, limitations,
and examples where applicable.

Changes to numerical algorithms, internal architecture, extension points, or
non-obvious implementation behaviour should also update the relevant 
[developer documentation](docs/developer/README.md).

Documentation should describe the behaviour that QMOMKit actually provides,
rather than only the implementation used in a particular example. Examples
should complement the documentation by showing the public functionality in a
concrete physical setting.

## Use of generative AI

The use of generative AI tools is permitted when contributing to QMOMKit. The
use of AI does not lower or raise the technical standard applied to a
contribution.

Contributors remain fully responsible for AI-assisted contributions and should:

- understand and be able to explain the submitted changes;
- review generated code and documentation rather than accepting them automatically;
- run the relevant tests and numerical verification;
- verify generated scientific claims, references, and numerical conclusions independently;
- ensure that generated material does not introduce incompatible, unlicensed, or improperly attributed third-party content;
- never provide confidential, proprietary, or otherwise restricted material to an AI service without appropriate authorization; and
- instruct any AI coding agent working in the repository to read and follow [AGENTS.md](AGENTS.md) before making changes.

Material use of generative AI should be mentioned in a _pull request_ when that
information is useful for reviewing the provenance, implementation process, or
validation of the contribution.

## Testing

Before submitting a contribution, run the tests relevant to the change and
report what was run.

For changes that affect common numerical infrastructure, broader regression
and verification testing may be appropriate.

New functionality must include suitable tests as part of the contribution;
relying solely on pre-existing tests is not sufficient.

Do not weaken tests, tolerances, or validation criteria merely to make a new
implementation pass without establishing why the change is justified.

## Pull requests

A pull request should explain:

- what changed;
- why the change is needed;
- any important numerical or API decisions;
- how the change was validated;
- which example demonstrates new functionality, where applicable;
- which user and developer documentation was added or updated; and
- any known limitations or follow-up work.

Keep pull requests as narrow as reasonably possible.