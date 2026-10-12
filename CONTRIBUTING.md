# Contributing to QMOMKit

Contributions to QMOMKit are welcome.

QMOMKit aims to remain a focused scientific computing library with clear
interfaces, minimal dependencies, and carefully tested numerical behavior. Any
changes should therefore preserve the project's focus on correctness, clarity,
numerical robustness, and maintainable implementation.

## Before making changes

Please:

- read the relevant documentation and existing implementation before modifying numerical behavior;
- keep changes focused on the problem being addressed;
- prefer existing QMOMKit functionality and configuration options over application-specific reimplementations;
- avoid unrelated cleanup in the same contribution;
- add or update tests and verification cases where they provide meaningful evidence for the change; and
- update the relevant user and developer documentation whenever behavior, interfaces, configuration, or numerical methods change.

## Numerical and scientific changes

Changes to numerical algorithms should be supported by appropriate evidence.
Depending on the change, this may include:

- unit tests for well-defined local behavior;
- regression tests protecting existing behavior;
- manufactured or analytical verification cases;
- independent numerical comparisons;
- convergence studies;
- comparison with established literature results.

Passing the existing test suite is necessary but is not, by itself, evidence
that a new numerical method or model is correct.

Scientific claims, references, and validation data should be checked against
primary sources wherever practical.
[Validation evidence](docs/developer/testing_and_validation.md#validation-evidence)
lists what a validation case must state.

## New features

New QMOMKit features must be accompanied by evidence that they work as
intended and by an example showing how they are used in a physically meaningful
setting.

A new feature should therefore include:

- suitable unit, regression, and/or verification tests, chosen according to
  the numerical and software behavior being introduced;
- a physical example demonstrating the feature through the public QMOMKit
  interface and configuration system; and
- appropriate documentation for both users and developers.

Tests should establish the relevant numerical and software contracts rather
than merely exercise code paths. In particular, new numerical functionality
should normally include verification or comparison against independently known
behavior where such evidence can reasonably be constructed.

The physical demonstration may be a new example or a clearly motivated
variant or extension of an existing example. It should exercise the feature
through the same public API and configuration options that users are
expected to employ. Example code should not duplicate, bypass, override, or
reimplement behavior that is already provided by QMOMKit.

A feature should not be considered complete until its implementation, tests,
example, and documentation are mutually consistent.

## Documentation

All contributions should be documented at the level appropriate to their
scope.

User-facing changes should update the relevant
[user documentation](docs/user/README.md), including configuration options,
public API behavior, assumptions, units, limitations, and examples where
applicable.

Changes to numerical algorithms, internal architecture, extension points, or
non-obvious implementation behavior should also update the relevant
[developer documentation](docs/developer/README.md).

Documentation should describe the behavior that QMOMKit actually provides,
rather than only the implementation used in a particular example. Examples
should complement the documentation by showing the public functionality in a
concrete physical setting.

## Use of generative AI

The use of generative AI tools is permitted when contributing to QMOMKit. The
use of AI does not lower or raise the technical standard applied to a
contribution. QMOMKit itself has been developed with substantial AI assistance;
[DEVELOPMENT_PHILOSOPHY.md](DEVELOPMENT_PHILOSOPHY.md) explains why that makes
evidence of correctness, rather than the origin of an implementation, the thing
that counts.

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

## Building, testing and formatting

Before submitting a contribution, run the tests relevant to the change and
report what was run. Run the core suite in a Debug build, because the unit tests
use `assert()`;
[Building and running](docs/developer/testing_and_validation.md#building-and-running)
has the commands, also for the installed-package and Monte Carlo checks. For
changes that affect common numerical infrastructure, broader regression and
verification testing may be appropriate.

Run `./tools/indent-all` before committing; CI fails if it would change a file.
It needs clang-format 22.1, which `./tools/download_clang_format` installs. CI
also builds with GCC and `-Werror`, so a compiler warning fails the build.

New functionality must include suitable tests as part of the contribution;
relying solely on pre-existing tests is not sufficient.

Do not weaken tests, tolerances, or validation criteria merely to make a new
implementation pass without establishing why the change is justified.

## Pull requests

A pull request uses the [PR template](.github/pull_request_template.md), which
lists what its description should explain and the license terms under which a
contribution is made. Keep pull requests as narrow as reasonably possible.