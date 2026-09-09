#pragma once

#include "qmomkit/phase_space_hierarchy.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace qmomkit
{

enum class ConditionalMomentStatus
{
    success,
    invalid_input,
    singular_system
};

struct ConditionalMomentOptions
{
    double singular_tolerance = 1.0e-14;
    double conditioning_warning_threshold = 1.0e12;
    double moore_penrose_rcond = 1.0e-14;
};

/**
 * Global raw moment tensor flattened in phase-space dimension order.
 *
 * `values[flatten_multi_index(alpha, extents)]` stores the raw mixed moment
 * with phase-space-order multi-index `alpha`.
 */
struct RawMomentTensor
{
    std::vector<std::size_t> extents;
    std::vector<double> values;
};

/**
 * Inputs for one-level conditional-moment reconstruction.
 *
 * `raw_moments` is a global raw moment tensor flattened in phase-space dimension order.
 * The selected `conditioning_level` identifies the parent dimension as
 * `conditioning_order[level]` and the target next dimension as
 * `conditioning_order[level + 1]`. All non-parent/non-target moment indices
 * are read at order zero.
 */
struct ConditionalMomentRequest
{
    PhaseSpaceHierarchy hierarchy;
    std::size_t conditioning_level = 0;
    std::vector<double> parent_abscissas;
    std::vector<double> parent_weights;
    RawMomentTensor raw_moments;
    ConditionalMomentOptions options;
};

struct ConditionalMomentDiagnostics
{
    /**
     * Expensive diagnostic quantities are reported when the caller requested
     * them. They are empty when diagnostic computation is disabled, so callers
     * should check availability before reading the values.
     */
    std::optional<double> residual_norm;
    std::optional<double> condition_estimate;
    bool conditioning_warning = false;
};

struct ConditionalMomentResult
{
    std::vector<double> conditional_moments;
    std::size_t parent_dimension = 0;
    std::size_t target_dimension = 0;
    std::size_t parent_node_count = 0;
    std::size_t target_moment_count = 0;
    ConditionalMomentStatus status = ConditionalMomentStatus::success;
    std::string message;
    ConditionalMomentDiagnostics diagnostics;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == ConditionalMomentStatus::success;
    }

    [[nodiscard]] double moment(std::size_t parent_node, std::size_t target_order) const;
};

/**
 * Reconstruct conditional moments for the next conditioning variable.
 *
 * Solves the direct weighted Vandermonde system
 *
 *     V(parent_abscissas)^T * diag(parent_weights) * C = M
 *
 * where rows of `M` are mixed raw moments for parent orders
 * `0..parent_node_count-1`, and columns correspond to target-dimension moment
 * orders. The returned conditional moments are flattened as
 * `{parent_node, target_order}` in row-major order.
 */
[[nodiscard]] ConditionalMomentResult reconstruct_conditional_moments(const ConditionalMomentRequest& request);

} // namespace qmomkit
