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
 * `success` and `quadrature.status` always agree: an unsuccessful selection
 * carries a non-success quadrature status and a message describing why. Callers
 * that only forward the quadrature therefore still report the failure.
 */
struct BoundedQuadratureSelection
{
    Quadrature1D quadrature;
    std::size_t selected_node_count_ceiling = 0;
    bool reduced = false;
    bool success = false;
    BoundedQuadratureCandidate diagnostic_candidate;
    std::size_t accepted_fallback_node_count = 0;
};

[[nodiscard]] BoundedQuadratureSelection select_bounded_quadrature(const std::vector<double>& moment_sequence,
                                                                   const BoundedQuadratureSelectionRequest& request, AbscissaBounds interval,
                                                                   const WheelerOptions& wheeler_options);

} // namespace qmomkit::internal
