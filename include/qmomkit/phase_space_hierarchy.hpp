#pragma once

#include "qmomkit/multi_index.hpp"
#include "qmomkit/qmom.hpp"

#include <cstddef>
#include <vector>

namespace qmomkit
{

/**
 * The QMOM family of one coordinate of a CQMOM hierarchy, and the settings of
 * the generalized family, as `QmomCase` holds them for QMOM. The settings of
 * the generalized family are ignored for any other family. CQMOM accepts
 * `classical` and `hyperbolic` coordinates.
 */
struct CoordinateFamily
{
    QmomFamily family = QmomFamily::classical;
    GeneralizedOptions generalized_options{};
};

/**
 * The families of the coordinates of a CQMOM hierarchy. `coordinates` is
 * either empty, which makes every coordinate `classical`, or holds one entry
 * per phase-space dimension, in phase-space dimension order.
 */
struct HierarchyFamilies
{
    std::vector<CoordinateFamily> coordinates{};
};

/**
 * Describes the phase-space dimensions, their families and the recursive
 * conditioning order.
 *
 * Node counts are stored by phase-space dimension. The node count of a
 * coordinate, `node_counts()` and `node_count(dimension)`, is the node count
 * of its reconstruction, which fixes the moments the coordinate requires,
 * `moment_extents()`. The maximum node count, `maximum_node_counts()`, is the
 * largest number of nodes the coordinate's quadrature holds, which fixes the
 * node structure: `total_nodes()`, the flat node indices and
 * `enumerate_quadrature_node_multi_indices`. For a `classical` coordinate, the
 * two counts are equal. A `hyperbolic` coordinate with node count n holds up
 * to n + 1 nodes. Flat node-index conversion uses phase-space dimension
 * order. The conditioning order is a separate permutation of phase-space
 * dimension indices and controls deterministic recursive traversal order.
 */
class PhaseSpaceHierarchy
{
public:
    explicit PhaseSpaceHierarchy(std::vector<std::size_t> node_counts);
    PhaseSpaceHierarchy(std::vector<std::size_t> node_counts, std::vector<std::size_t> conditioning_order);

    /**
     * A hierarchy whose coordinates have the given families. Throws
     * `std::invalid_argument` when `families.coordinates` is neither empty nor
     * of one entry per phase-space dimension, and when a coordinate's family
     * is `generalized`.
     */
    PhaseSpaceHierarchy(std::vector<std::size_t> node_counts, std::vector<std::size_t> conditioning_order, HierarchyFamilies families);

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
    [[nodiscard]] const std::vector<std::size_t>& maximum_node_counts() const noexcept
    {
        return maximum_node_counts_.empty() ? node_counts_ : maximum_node_counts_;
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
     * The families of the coordinates, which are empty when every coordinate is
     * `classical`, also when the hierarchy was built with `classical` entries,
     * and otherwise hold one entry per phase-space dimension.
     */
    [[nodiscard]] const HierarchyFamilies& families() const noexcept;

    /**
     * Extents for the raw moment tensor needed by recursive CQMOM: per
     * phase-space dimension, the number of moments the family of the
     * coordinate reads, 2 * node_count for `classical` and
     * 2 * node_count + 1 for `hyperbolic`.
     */
    [[nodiscard]] std::vector<std::size_t> moment_extents() const;

    /**
     * Convert between flat and multi-index node notation in phase-space dimension
     * order, independent of the conditioning order, over the maximum node
     * counts.
     */
    [[nodiscard]] std::size_t flatten_node_index(MultiIndexView node_index) const;
    [[nodiscard]] MultiIndex unflatten_node_index(std::size_t flat_index) const;

private:
    // A hierarchy is copied into every layout, quadrature and state, several
    // times per step of an evolution. When every coordinate is classical,
    // families_ and maximum_node_counts_ stay empty, so that a copy allocates
    // no more than the node counts and the conditioning order.
    std::vector<std::size_t> node_counts_;
    std::vector<std::size_t> conditioning_order_;
    HierarchyFamilies families_;
    std::vector<std::size_t> maximum_node_counts_;
    std::size_t total_nodes_ = 0;
};

/**
 * Enumerate quadrature-node multi-indices in recursive conditioning order.
 */
[[nodiscard]] std::vector<MultiIndex> enumerate_quadrature_node_multi_indices(const PhaseSpaceHierarchy& hierarchy);

} // namespace qmomkit
