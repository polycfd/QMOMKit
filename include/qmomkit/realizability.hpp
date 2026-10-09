#pragma once

#include "qmomkit/moment_set.hpp"

#include <cstddef>
#include <string>

namespace qmomkit
{

/**
 * One-dimensional support used for truncated raw-moment realizability checks.
 */
enum class MomentSupport
{
    real_line,
    nonnegative,
    bounded_interval
};

/**
 * Explicit closed interval for the one-dimensional Hausdorff moment problem.
 */
struct BoundedInterval
{
    double lower = 0.0;
    double upper = 1.0;
};

/**
 * Classification for a one-dimensional truncated raw-moment sequence.
 */
enum class RealizabilityStatus
{
    strictly_realizable,
    boundary,
    not_realizable,
    invalid_input,
    numerically_indeterminate
};

/**
 * Numerical controls for support-aware realizability analysis.
 *
 * With relative_tolerance, the moment matrices are tested in the standardized
 * coordinate (x - m1/m0) / s, with s about the spread of the distribution, and
 * tolerance is relative to those matrices; the classification then depends
 * neither on the location of the distribution nor on the unit of its
 * coordinate. Eigenvalues are also compared against the rounding of the input
 * moments propagated into that coordinate. Without relative_tolerance, the raw
 * matrices are tested against tolerance as an absolute threshold. In either
 * mode, an eigenvalue within 1000 epsilon times the norm of the tested matrix
 * counts as zero, so a tolerance below that resolves no finer.
 */
struct RealizabilityOptions
{
    double tolerance = 1.0e-12;
    bool relative_tolerance = true;
};

/**
 * Diagnostic result for a support-aware one-dimensional realizability check.
 *
 * minimum_pivot is the limiting eigenvalue, or kernel projection, and
 * tolerance_used the threshold it was compared against; with
 * relative_tolerance both refer to the standardized matrices.
 */
struct RealizabilityResult
{
    RealizabilityStatus status = RealizabilityStatus::invalid_input;
    std::size_t checked_order = 0;
    std::size_t failed_matrix_order = 0;
    std::string failed_matrix;
    double minimum_pivot = 0.0;
    double tolerance_used = 0.0;
    std::string message;

    [[nodiscard]] bool realizable() const noexcept;
    [[nodiscard]] bool on_boundary() const noexcept;
};

/**
 * Check a one-dimensional raw moment sequence against an explicit support.
 *
 * The input is interpreted as the inclusive sequence m_0, ..., m_N. This
 * function is opt-in and is independent of Wheeler/QMOM inversion.
 */
[[nodiscard]] RealizabilityResult check_realizability(const MomentSet& moments, MomentSupport support, const RealizabilityOptions& options = {});

/**
 * Check a one-dimensional raw moment sequence on an explicit bounded interval.
 *
 * This overload implies MomentSupport::bounded_interval.
 */
[[nodiscard]] RealizabilityResult check_realizability(const MomentSet& moments, BoundedInterval interval, const RealizabilityOptions& options = {});

} // namespace qmomkit
