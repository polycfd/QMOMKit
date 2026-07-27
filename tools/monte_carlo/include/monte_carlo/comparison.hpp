#pragma once

#include "monte_carlo/status.hpp"
#include "qmomkit/conditional_moments.hpp"
#include "qmomkit/moment_set.hpp"
#include "qmomkit/multi_index.hpp"
#include "qmomkit/runge_kutta.hpp"

#include <cstddef>
#include <iosfwd>
#include <optional>
#include <string>
#include <vector>

namespace qmomkit::monte_carlo
{

/**
 * One raw-moment comparison between a deterministic reference value and a Monte
 * Carlo sample estimate.
 *
 * `reference_value` stores the first `compare_moments` argument, named
 * `reference`; `computed_value` stores the sample estimate. The signed error is:
 *
 *     computed_value - reference_value
 *
 * `absolute_error` stores the magnitude of that signed error. A z-score measures
 * the signed error in units of the supplied Monte Carlo standard error. Large
 * z-scores are diagnostic signals, not automatic proof that either model is
 * wrong.
 */
struct MomentComparisonRow
{
    qmomkit::MultiIndex multi_index;
    std::size_t order = 0;
    double reference_value = 0.0;
    double computed_value = 0.0;
    double signed_error = 0.0;
    double absolute_error = 0.0;
    double relative_error = 0.0;
    double standard_error = 0.0;
    double z_score = 0.0;
};

/**
 * Collection of moment-comparison rows plus validation status.
 *
 * On failure, `rows` is empty and `message` describes the invalid input or
 * non-finite compared quantity.
 */
struct MomentComparisonResult
{
    std::vector<MomentComparisonRow> rows;
    ComparisonStatus status = ComparisonStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == ComparisonStatus::success;
    }
};

/**
 * Metadata written with comparison rows to make a sample-based result
 * reproducible and interpretable.
 */
struct ComparisonOutputMetadata
{
    double time = 0.0;
    std::size_t sample_count = 0;
    unsigned long long seed = 0;
    double number_density = 0.0;
    std::optional<qmomkit::AdaptiveRk23Options> rk23_options;
};

/**
 * Compare one-dimensional deterministic reference moments with Monte Carlo
 * sample estimates.
 *
 * `reference` is the deterministic or otherwise chosen baseline moment set.
 * `monte_carlo` is the sample-estimated moment set. `standard_errors` must
 * contain the Monte Carlo standard error for each compared moment. The returned
 * rows use the sign convention `computed - reference` for `signed_error` and
 * z-scores, and `absolute_error` stores the magnitude of that difference.
 */
[[nodiscard]] MomentComparisonResult compare_moments(const qmomkit::MomentSet& reference, const qmomkit::MomentSet& monte_carlo,
                                                     const std::vector<double>& standard_errors);

/**
 * Compare multidimensional deterministic reference moments with Monte Carlo
 * sample estimates.
 *
 * Tensor extents and flat storage order must match exactly. Multi-indices in
 * the result follow QMOMKit phase-space-order enumeration. The sign convention
 * is the same as the one-dimensional overload: `computed - reference`.
 */
[[nodiscard]] MomentComparisonResult compare_moments(const qmomkit::RawMomentTensor& reference, const qmomkit::RawMomentTensor& monte_carlo,
                                                     const std::vector<double>& standard_errors);

/**
 * Write comparison rows and reproducibility metadata without changing the
 * statistical interpretation: the reference value is the first compared input,
 * and the Monte Carlo value is the finite-sample estimate.
 */
void write_comparison_text(std::ostream& output, const MomentComparisonResult& result, const ComparisonOutputMetadata& metadata);

} // namespace qmomkit::monte_carlo
