#pragma once

#include "performance.hpp"
#include "qmomkit/cqmom.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace qmomkit::internal
{

[[nodiscard]] inline ConditioningPath parent_path_at_level(const ConditioningPathView terminal_path, const std::size_t conditioning_level)
{
    ConditioningPath parent;
    for (std::size_t level = 0; level < conditioning_level; ++level)
    {
        parent = parent.with_child(terminal_path.node_index_at_level(level));
    }
    return parent;
}

struct ActiveBranchAccumulator
{
    std::size_t conditioning_level = 0;
    ConditioningPath parent_path;
    std::size_t phase_space_dimension = 0;
    std::size_t configured_node_count = 0;
    std::vector<std::size_t> child_node_indices;
};

[[nodiscard]] inline ActiveBranchAccumulator* find_active_branch_accumulator(std::vector<ActiveBranchAccumulator>& counts,
                                                                             const std::size_t conditioning_level,
                                                                             const ConditioningPath& parent_path,
                                                                             const std::size_t phase_space_dimension)
{
    for (ActiveBranchAccumulator& count : counts)
    {
        if (count.conditioning_level == conditioning_level && count.phase_space_dimension == phase_space_dimension &&
            count.parent_path == parent_path)
        {
            return &count;
        }
    }
    return nullptr;
}

[[nodiscard]] inline bool contains_node_index(const std::vector<std::size_t>& nodes, const std::size_t node_index)
{
    return std::find(nodes.begin(), nodes.end(), node_index) != nodes.end();
}

inline void validate_contiguous_child_indices(const ActiveBranchAccumulator& accumulator)
{
    std::vector<std::size_t> sorted = accumulator.child_node_indices;
    std::sort(sorted.begin(), sorted.end());
    for (std::size_t expected = 0; expected < sorted.size(); ++expected)
    {
        if (sorted[expected] != expected)
        {
            throw std::invalid_argument("CQMOM adaptive quadrature terminal paths must use contiguous branch-local child node indices.");
        }
    }
}

[[nodiscard]] inline std::vector<CqmomActiveBranchCount> active_branch_counts_from_quadrature(const CqmomQuadrature& quadrature)
{
    if (!quadrature.adaptive_branch_counts_enabled())
    {
        return {};
    }
    const PhaseSpaceHierarchy& hierarchy = quadrature.hierarchy();
    const std::size_t dimensions = hierarchy.number_of_dimensions();
    if (dimensions == 0)
    {
        throw std::invalid_argument("CQMOM adaptive quadrature must have at least one phase-space dimension.");
    }
    if (quadrature.terminal_paths().size() != quadrature.weights().size())
    {
        throw std::invalid_argument("CQMOM adaptive quadrature terminal-path count must match its weight count.");
    }

    std::vector<ActiveBranchAccumulator> accumulators;
    for (const ConditioningPathView terminal_path : quadrature.terminal_paths())
    {
        if (terminal_path.size() != dimensions)
        {
            throw std::invalid_argument("CQMOM adaptive quadrature terminal paths must contain one node index per conditioning level.");
        }
        for (std::size_t level = 0; level < dimensions; ++level)
        {
            const std::size_t phase_space_dimension = hierarchy.conditioning_order()[level];
            const std::size_t child_node_index = terminal_path.node_index_at_level(level);
            const std::size_t configured_node_count = hierarchy.node_count(phase_space_dimension);
            if (child_node_index >= configured_node_count)
            {
                throw std::invalid_argument("CQMOM adaptive quadrature terminal path contains a child node index outside the configured capacity.");
            }

            const ConditioningPath parent_path = parent_path_at_level(terminal_path, level);
            ActiveBranchAccumulator* accumulator = find_active_branch_accumulator(accumulators, level, parent_path, phase_space_dimension);
            if (accumulator == nullptr)
            {
                accumulators.push_back({level, parent_path, phase_space_dimension, configured_node_count, {}});
                accumulator = &accumulators.back();
            }
            if (!contains_node_index(accumulator->child_node_indices, child_node_index))
            {
                accumulator->child_node_indices.push_back(child_node_index);
            }
        }
    }

    std::vector<CqmomActiveBranchCount> counts;
    counts.reserve(accumulators.size());
    for (const ActiveBranchAccumulator& accumulator : accumulators)
    {
        validate_contiguous_child_indices(accumulator);
        counts.push_back({accumulator.conditioning_level, accumulator.parent_path, accumulator.phase_space_dimension,
                          accumulator.configured_node_count, accumulator.child_node_indices.size()});
    }
    performance::record_active_count_derivation();
    return counts;
}

} // namespace qmomkit::internal
