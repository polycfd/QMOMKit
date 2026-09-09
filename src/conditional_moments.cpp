#include "qmomkit/conditional_moments.hpp"

#include "internal/linear_algebra.hpp"
#include "internal/numeric_utils.hpp"
#include "internal/vandermonde.hpp"
#include "qmomkit/multi_index.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace qmomkit
{
namespace
{

using internal::DenseMatrix;
using internal::LinearSolveResult;

[[nodiscard]] ConditionalMomentResult invalid_result(const ConditionalMomentRequest& request, ConditionalMomentStatus status, std::string message)
{
    ConditionalMomentResult result;
    result.status = status;
    result.message = std::move(message);
    const std::size_t dimensions = request.hierarchy.number_of_dimensions();
    if (request.conditioning_level < dimensions && request.conditioning_level + 1 < dimensions)
    {
        result.parent_dimension = request.hierarchy.conditioning_order()[request.conditioning_level];
        result.target_dimension = request.hierarchy.conditioning_order()[request.conditioning_level + 1];
    }
    result.parent_node_count = request.parent_abscissas.size();
    return result;
}

[[nodiscard]] bool valid_options(const ConditionalMomentOptions& options)
{
    return std::isfinite(options.singular_tolerance) && options.singular_tolerance > 0.0 && std::isfinite(options.conditioning_warning_threshold) &&
           options.conditioning_warning_threshold > 0.0 && std::isfinite(options.moore_penrose_rcond) && options.moore_penrose_rcond >= 0.0 &&
           options.moore_penrose_rcond < 1.0;
}

[[nodiscard]] DenseMatrix extract_rhs(const ConditionalMomentRequest& request, const std::size_t parent_dimension, const std::size_t target_dimension)
{
    // RHS rows are parent moment orders and columns are target moment orders;
    // all other dimensions stay at order zero for this one-level solve.
    const std::size_t n_parent = request.parent_abscissas.size();
    const std::size_t n_target_moments = request.raw_moments.extents[target_dimension];
    DenseMatrix rhs(n_parent, n_target_moments);

    MultiIndex moment_index(request.raw_moments.extents.size(), 0);
    for (std::size_t parent_order = 0; parent_order < n_parent; ++parent_order)
    {
        for (std::size_t target_order = 0; target_order < n_target_moments; ++target_order)
        {
            std::fill(moment_index.begin(), moment_index.end(), 0);
            moment_index[parent_dimension] = parent_order;
            moment_index[target_dimension] = target_order;
            rhs(parent_order, target_order) = request.raw_moments.values[flatten_multi_index(moment_index, request.raw_moments.extents)];
        }
    }

    return rhs;
}

} // namespace

double ConditionalMomentResult::moment(const std::size_t parent_node, const std::size_t target_order) const
{
    if (parent_node >= parent_node_count || target_order >= target_moment_count)
    {
        throw std::out_of_range("Conditional moment index is out of range.");
    }
    return conditional_moments.at(parent_node * target_moment_count + target_order);
}

ConditionalMomentResult reconstruct_conditional_moments(const ConditionalMomentRequest& request)
{
    if (!valid_options(request.options))
    {
        return invalid_result(request, ConditionalMomentStatus::invalid_input, "Conditional moment options are invalid.");
    }
    const std::size_t dimensions = request.hierarchy.number_of_dimensions();
    if (request.conditioning_level >= dimensions || request.conditioning_level + 1 >= dimensions)
    {
        return invalid_result(request, ConditionalMomentStatus::invalid_input, "Conditioning level must have a next target dimension.");
    }
    if (request.parent_abscissas.empty() || request.parent_abscissas.size() != request.parent_weights.size())
    {
        return invalid_result(request, ConditionalMomentStatus::invalid_input, "Parent abscissas and weights must be non-empty and same-sized.");
    }
    if (!internal::all_finite(request.parent_abscissas) || !internal::all_finite(request.parent_weights) ||
        !internal::all_finite(request.raw_moments.values))
    {
        return invalid_result(request, ConditionalMomentStatus::invalid_input, "Conditional moment inputs must be finite.");
    }
    if (request.raw_moments.extents.size() != request.hierarchy.number_of_dimensions())
    {
        return invalid_result(request, ConditionalMomentStatus::invalid_input, "Mixed moment extents must match the hierarchy dimension.");
    }

    std::size_t expected_moments = 0;
    try
    {
        expected_moments = product_extents(request.raw_moments.extents);
    }
    catch (const std::invalid_argument& ex)
    {
        return invalid_result(request, ConditionalMomentStatus::invalid_input, ex.what());
    }
    if (request.raw_moments.values.size() != expected_moments)
    {
        return invalid_result(request, ConditionalMomentStatus::invalid_input, "Mixed moment data size does not match its extents.");
    }

    const std::size_t parent_dimension = request.hierarchy.conditioning_order()[request.conditioning_level];
    const std::size_t target_dimension = request.hierarchy.conditioning_order()[request.conditioning_level + 1];
    const std::size_t n_parent = request.parent_abscissas.size();
    if (request.raw_moments.extents[parent_dimension] < n_parent)
    {
        return invalid_result(request, ConditionalMomentStatus::invalid_input, "Mixed moments do not contain enough parent moment orders.");
    }

    // Solve weighted parent-node interpolation for the conditional moments of
    // the next target dimension at each parent abscissa.
    DenseMatrix system = internal::make_conditional_vandermonde_matrix(request.parent_abscissas, request.parent_weights);
    DenseMatrix rhs = extract_rhs(request, parent_dimension, target_dimension);
    LinearSolveResult solve = internal::solve_dense(system, rhs, request.options.singular_tolerance, "Weighted Vandermonde system is singular.");
    if (!solve.success)
    {
        return invalid_result(request, ConditionalMomentStatus::singular_system, solve.message);
    }

    const double residual = internal::residual_norm(system, solve.solution, rhs);
    const double condition = internal::condition_estimate(system, request.options.singular_tolerance);
    const std::size_t target_moment_count = request.raw_moments.extents[target_dimension];

    ConditionalMomentResult result;
    result.parent_dimension = parent_dimension;
    result.target_dimension = target_dimension;
    result.parent_node_count = n_parent;
    result.target_moment_count = target_moment_count;
    result.status = ConditionalMomentStatus::success;
    result.message = "Conditional moments reconstructed.";
    result.diagnostics.residual_norm = residual;
    result.diagnostics.condition_estimate = condition;
    result.diagnostics.conditioning_warning = condition > request.options.conditioning_warning_threshold;
    result.conditional_moments.assign(n_parent * target_moment_count, 0.0);
    for (std::size_t node = 0; node < n_parent; ++node)
    {
        for (std::size_t order = 0; order < target_moment_count; ++order)
        {
            result.conditional_moments[node * target_moment_count + order] = solve.solution(node, order);
        }
    }
    return result;
}

} // namespace qmomkit
