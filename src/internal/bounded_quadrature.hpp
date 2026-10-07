#pragma once

#include "qmomkit/abscissa_bounds.hpp"
#include "qmomkit/quadrature.hpp"
#include "qmomkit/wheeler.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace qmomkit::internal
{

struct BoundedQuadratureCandidate
{
    std::size_t node_count = 0;
    std::string attempt_type;
    Quadrature1D quadrature;
    bool inversion_attempted = false;
    bool accepted = false;
    std::string message;
};

struct BoundedQuadratureSelectionRequest
{
    std::size_t initial_node_count = 0;
    std::size_t configured_node_count = 0;
    std::optional<std::size_t> previous_active_node_count;
    bool adaptive = false;
    bool allow_growth = false;
};

/**
 * Outcome of a bounded-support quadrature selection.
 *
 * `quadrature.status` is the single record of whether a rule was selected: an
 * unsuccessful selection carries a non-success status and a message describing
 * why. The remaining fields are diagnostic and describe how the selection was
 * reached; they do not restate the outcome.
 */
struct BoundedQuadratureSelection
{
    Quadrature1D quadrature;
    std::size_t selected_node_count_ceiling = 0;
    bool reduced = false;
    BoundedQuadratureCandidate diagnostic_candidate;
    std::size_t accepted_fallback_node_count = 0;
};

/**
 * Wheeler inversion of conditional moments, which a solve produced from other
 * moments. Their rounding is relative to `coordinate_scale`, the scale of
 * their coordinate, rather than to their own magnitudes, and adaptive node
 * reduction allows for it. A zero scale treats the moments as data, as the
 * public `compute_quadrature(...)` does.
 */
[[nodiscard]] Quadrature1D compute_conditional_quadrature(const MomentSet& moments, std::size_t max_nodes, const WheelerOptions& options,
                                                          double coordinate_scale);

/**
 * `coordinate_scale` is passed to every Wheeler inversion of the selection;
 * see `compute_conditional_quadrature(...)`.
 */
[[nodiscard]] BoundedQuadratureSelection select_bounded_quadrature(const std::vector<double>& moment_sequence,
                                                                   const BoundedQuadratureSelectionRequest& request, AbscissaBounds interval,
                                                                   const WheelerOptions& wheeler_options, double coordinate_scale = 0.0);

} // namespace qmomkit::internal
