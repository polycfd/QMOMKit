#pragma once

#include <optional>

namespace qmomkit
{

/**
 * Outcome of one conditional-moment solve inside recursive CQMOM.
 *
 * Reconstruction reports the status of each conditioning branch through
 * `CqmomDiagnostics::conditional_moment_status`. `invalid_input` covers a
 * malformed branch slice or parent quadrature; `singular_system` means the
 * weighted Vandermonde solve for that branch failed.
 */
enum class ConditionalMomentStatus
{
    success,
    invalid_input,
    singular_system
};

/**
 * Numerical controls for the conditional-moment solves performed by CQMOM.
 *
 * `singular_tolerance` is the relative pivot floor for the conditional solve
 * and for its condition estimate. `conditioning_warning_threshold` is compared
 * against that estimate, so it has an effect only where the estimate is
 * computed. `moore_penrose_rcond` is the truncation cutoff used by
 * `VandermondeSolver::moore_penrose` and is ignored by the other solvers.
 *
 * `singular_tolerance` and `conditioning_warning_threshold` must be finite and
 * positive, and `moore_penrose_rcond` must be finite and in `[0, 1)`. A case
 * carrying values outside those ranges is rejected as invalid input before
 * reconstruction begins.
 */
struct ConditionalMomentOptions
{
    double singular_tolerance = 1.0e-14;
    double conditioning_warning_threshold = 1.0e12;
    double moore_penrose_rcond = 1.0e-14;
};

/**
 * Optional per-branch diagnostics from a conditional-moment solve.
 *
 * Every member is populated only at `CqmomDiagnosticLevel::detailed`, because
 * the condition estimate costs a full inverse of the conditional Vandermonde
 * matrix. Below that level the optionals stay empty and `conditioning_warning`
 * stays false, so a false warning is not by itself evidence of good
 * conditioning.
 */
struct ConditionalMomentDiagnostics
{
    std::optional<double> residual_norm;
    std::optional<double> condition_estimate;
    bool conditioning_warning = false;
};

} // namespace qmomkit
