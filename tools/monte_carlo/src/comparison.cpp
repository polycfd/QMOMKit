#include "monte_carlo/comparison.hpp"

#include "qmomkit/multi_index.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <limits>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

namespace qmomkit::monte_carlo
{
namespace
{

[[nodiscard]] double relative_error(const double absolute_error, const double reference, const double standard_error)
{
    // Use the larger of reference scale and sampling uncertainty so near-zero
    // moments with visible Monte Carlo noise do not produce misleading ratios.
    constexpr double floor = 1.0e-300;
    const double scale = std::max({std::abs(reference), standard_error, floor});
    return absolute_error / scale;
}

[[nodiscard]] double z_score(const double error, const double standard_error)
{
    // Z-scores scale the signed sample-reference difference by Monte Carlo
    // sampling uncertainty; they are diagnostics, not pass/fail decisions.
    if (standard_error == 0.0)
    {
        if (error == 0.0)
        {
            return 0.0;
        }
        return std::copysign(std::numeric_limits<double>::infinity(), error);
    }
    return error / standard_error;
}

[[nodiscard]] bool all_finite_moments(const std::vector<double>& values)
{
    return std::all_of(values.begin(), values.end(), [](const double value) { return std::isfinite(value); });
}

[[nodiscard]] bool has_zero_extent(const std::vector<std::size_t>& extents)
{
    return std::any_of(extents.begin(), extents.end(), [](const std::size_t extent) { return extent == 0; });
}

[[nodiscard]] ComparisonStatus validate_standard_errors(const std::vector<double>& standard_errors, const std::size_t expected_size,
                                                        std::string& message)
{
    if (standard_errors.size() != expected_size)
    {
        message = "Standard-error count must match moment count.";
        return ComparisonStatus::invalid_input;
    }
    for (const double standard_error : standard_errors)
    {
        if (!std::isfinite(standard_error) || standard_error < 0.0)
        {
            message = "Standard errors must be finite and non-negative.";
            return ComparisonStatus::invalid_input;
        }
    }
    return ComparisonStatus::success;
}

void write_rk23_metadata(std::ostream& output, const std::optional<qmomkit::AdaptiveRk23Options>& options)
{
    if (!options.has_value())
    {
        return;
    }
    output << "# rk23_relative_tolerance=" << options->relative_tolerance << '\n';
    output << "# rk23_absolute_tolerance=" << options->absolute_tolerance << '\n';
    output << "# rk23_max_step=" << options->max_step << '\n';
    output << "# rk23_min_step=" << options->min_step << '\n';
    output << "# rk23_max_steps=" << options->max_steps << '\n';
    if (options->initial_step.has_value())
    {
        output << "# rk23_initial_step=" << *options->initial_step << '\n';
    }
}

} // namespace

MomentComparisonResult compare_moments(const qmomkit::MomentSet& reference, const qmomkit::MomentSet& monte_carlo,
                                       const std::vector<double>& standard_errors)
{
    if (reference.empty() || monte_carlo.empty())
    {
        return MomentComparisonResult{{}, ComparisonStatus::invalid_input, "Moment sets must not be empty."};
    }
    if (reference.size() != monte_carlo.size())
    {
        return MomentComparisonResult{{}, ComparisonStatus::invalid_input, "Reference and Monte Carlo moment counts must match."};
    }

    std::string message;
    const ComparisonStatus error_status = validate_standard_errors(standard_errors, reference.size(), message);
    if (error_status != ComparisonStatus::success)
    {
        return MomentComparisonResult{{}, error_status, std::move(message)};
    }
    if (!all_finite_moments(reference.values()) || !all_finite_moments(monte_carlo.values()))
    {
        return MomentComparisonResult{{}, ComparisonStatus::moment_failed, "Compared moment values must be finite."};
    }

    std::vector<MomentComparisonRow> rows;
    rows.reserve(reference.size());
    for (std::size_t order = 0; order < reference.size(); ++order)
    {
        const double reference_value = reference[order];
        const double computed_value = monte_carlo[order];
        const double signed_error = computed_value - reference_value;
        const double absolute = std::abs(signed_error);
        const double standard_error = standard_errors[order];
        rows.push_back(MomentComparisonRow{qmomkit::MultiIndex{order}, order, reference_value, computed_value, signed_error, absolute,
                                           relative_error(absolute, reference_value, standard_error), standard_error,
                                           z_score(signed_error, standard_error)});
    }

    return MomentComparisonResult{std::move(rows), ComparisonStatus::success, ""};
}

MomentComparisonResult compare_moments(const qmomkit::RawMomentTensor& reference, const qmomkit::RawMomentTensor& monte_carlo,
                                       const std::vector<double>& standard_errors)
{
    if (reference.extents.empty() || monte_carlo.extents.empty())
    {
        return MomentComparisonResult{{}, ComparisonStatus::invalid_input, "Moment tensors must not be empty."};
    }
    if (has_zero_extent(reference.extents) || has_zero_extent(monte_carlo.extents))
    {
        return MomentComparisonResult{{}, ComparisonStatus::invalid_input, "Moment tensor extents must be positive."};
    }
    if (reference.extents != monte_carlo.extents)
    {
        return MomentComparisonResult{{}, ComparisonStatus::invalid_input, "Reference and Monte Carlo tensor extents must match."};
    }
    std::size_t moment_count = 0;
    try
    {
        moment_count = qmomkit::product_extents(reference.extents);
    }
    catch (const std::exception& error)
    {
        return MomentComparisonResult{{}, ComparisonStatus::invalid_input, error.what()};
    }
    if (reference.values.size() != moment_count || monte_carlo.values.size() != moment_count)
    {
        return MomentComparisonResult{{}, ComparisonStatus::invalid_input, "Moment tensor value count must match extents."};
    }

    std::string message;
    const ComparisonStatus error_status = validate_standard_errors(standard_errors, moment_count, message);
    if (error_status != ComparisonStatus::success)
    {
        return MomentComparisonResult{{}, error_status, std::move(message)};
    }
    if (!all_finite_moments(reference.values) || !all_finite_moments(monte_carlo.values))
    {
        return MomentComparisonResult{{}, ComparisonStatus::moment_failed, "Compared tensor moment values must be finite."};
    }

    const std::vector<qmomkit::MultiIndex> orders = qmomkit::enumerate_multi_indices(reference.extents);
    std::vector<MomentComparisonRow> rows;
    rows.reserve(moment_count);
    for (std::size_t flat = 0; flat < moment_count; ++flat)
    {
        const double reference_value = reference.values[flat];
        const double computed_value = monte_carlo.values[flat];
        const double signed_error = computed_value - reference_value;
        const double absolute = std::abs(signed_error);
        const double standard_error = standard_errors[flat];
        rows.push_back(MomentComparisonRow{orders[flat], flat, reference_value, computed_value, signed_error, absolute,
                                           relative_error(absolute, reference_value, standard_error), standard_error,
                                           z_score(signed_error, standard_error)});
    }

    return MomentComparisonResult{std::move(rows), ComparisonStatus::success, ""};
}

void write_comparison_text(std::ostream& output, const MomentComparisonResult& result, const ComparisonOutputMetadata& metadata)
{
    output << std::setprecision(17);
    output << "# time=" << metadata.time << '\n';
    output << "# sample_count=" << metadata.sample_count << '\n';
    output << "# seed=" << metadata.seed << '\n';
    output << "# number_density=" << metadata.number_density << '\n';
    write_rk23_metadata(output, metadata.rk23_options);
    // The repeated metadata columns make each comparison row self-contained
    // when rows are filtered by moment or time.
    output << "# time sample_count seed number_density";

    std::size_t dimensions = 0;
    if (!result.rows.empty())
    {
        dimensions = result.rows.front().multi_index.size();
    }
    for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
    {
        output << " alpha_" << dimension;
    }
    output << " reference_value computed_value signed_error absolute_error relative_error standard_error z_score\n";

    if (!result.ok())
    {
        return;
    }
    for (const MomentComparisonRow& row : result.rows)
    {
        output << metadata.time << ' ' << metadata.sample_count << ' ' << metadata.seed << ' ' << metadata.number_density;
        for (const std::size_t component : row.multi_index)
        {
            output << ' ' << component;
        }
        output << ' ' << row.reference_value << ' ' << row.computed_value << ' ' << row.signed_error << ' ' << row.absolute_error << ' '
               << row.relative_error << ' ' << row.standard_error << ' ' << row.z_score << '\n';
    }
}

} // namespace qmomkit::monte_carlo
