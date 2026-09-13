#pragma once

#include "qmomkit/multi_index.hpp"

#include <cstddef>
#include <vector>

namespace qmomkit
{

/**
 * Describes the phase-space dimensions and recursive conditioning order.
 *
 * Node counts are stored by phase-space dimension. Flat node-index
 * conversion also uses phase-space dimension order. The conditioning order is a
 * separate permutation of phase-space dimension indices and controls deterministic
 * recursive traversal order.
 */
class PhaseSpaceHierarchy
{
public:
    explicit PhaseSpaceHierarchy(std::vector<std::size_t> node_counts);
    PhaseSpaceHierarchy(std::vector<std::size_t> node_counts, std::vector<std::size_t> conditioning_order);

    // Defined here rather than in the library so that callers can inline them.
    // These are read in the innermost reconstruction loops, where an
    // out-of-line call for a stored value is a measurable share of the work.
    [[nodiscard]] std::size_t number_of_dimensions() const noexcept
    {
        return node_counts_.size();
    }
    [[nodiscard]] const std::vector<std::size_t>& node_counts() const noexcept
    {
        return node_counts_;
    }
    [[nodiscard]] const std::vector<std::size_t>& conditioning_order() const noexcept
    {
        return conditioning_order_;
    }
    [[nodiscard]] std::size_t total_nodes() const noexcept
    {
        return total_nodes_;
    }

    [[nodiscard]] std::size_t node_count(std::size_t dimension) const;

    /**
     * Extents for the raw moment tensor needed by recursive CQMOM.
     *
     * Recursive CQMOM passes 1D slices of length 2*n to Wheeler in each
     * dimension, so the required tensor convention is 2 * node_count per
     * phase-space dimension.
     */
    [[nodiscard]] std::vector<std::size_t> moment_extents() const;

    /**
     * Convert between flat and multi-index node notation in phase-space dimension
     * order, independent of the conditioning order.
     */
    [[nodiscard]] std::size_t flatten_node_index(MultiIndexView node_index) const;
    [[nodiscard]] MultiIndex unflatten_node_index(std::size_t flat_index) const;

private:
    std::vector<std::size_t> node_counts_;
    std::vector<std::size_t> conditioning_order_;
    std::size_t total_nodes_ = 0;
};

/**
 * Enumerate quadrature-node multi-indices in recursive conditioning order.
 */
[[nodiscard]] std::vector<MultiIndex> enumerate_quadrature_node_multi_indices(const PhaseSpaceHierarchy& hierarchy);

} // namespace qmomkit
