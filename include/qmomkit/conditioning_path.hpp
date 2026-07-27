#pragma once

#include "qmomkit/multi_index.hpp"
#include "qmomkit/phase_space_hierarchy.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace qmomkit
{

/**
 * Node path through recursive conditioning levels.
 *
 * Entries are ordered by conditioning level, not phase-space dimension. For a
 * hierarchy with conditioning order `{2, 0, 1}`, `node_indices_by_level()[1]`
 * is the selected node index in phase-space dimension `0`.
 */
class ConditioningPath
{
public:
    ConditioningPath() = default;
    explicit ConditioningPath(std::vector<std::size_t> node_indices_by_level);

    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] const std::vector<std::size_t>& node_indices_by_level() const noexcept;
    [[nodiscard]] std::size_t node_index_at_level(std::size_t conditioning_level) const;

    [[nodiscard]] ConditioningPath with_child(std::size_t child_node_index) const;
    [[nodiscard]] MultiIndex phase_space_node_index(const PhaseSpaceHierarchy& hierarchy) const;
    [[nodiscard]] std::string to_string() const;
    [[nodiscard]] std::string describe(const PhaseSpaceHierarchy& hierarchy) const;

private:
    std::vector<std::size_t> node_indices_by_level_;
};

[[nodiscard]] bool operator==(const ConditioningPath& lhs, const ConditioningPath& rhs) noexcept;
[[nodiscard]] bool operator!=(const ConditioningPath& lhs, const ConditioningPath& rhs) noexcept;
[[nodiscard]] bool operator<(const ConditioningPath& lhs, const ConditioningPath& rhs) noexcept;

/**
 * Map a complete recursive path to the phase-space-order quadrature node index.
 *
 * The path must contain one node index per conditioning level. The returned
 * multi-index is ordered by phase-space dimension and is suitable for
 * rectangular `CqmomQuadrature` accessors.
 */
[[nodiscard]] MultiIndex recursive_path_to_node_index(const PhaseSpaceHierarchy& hierarchy, const ConditioningPath& path);

[[nodiscard]] std::size_t recursive_path_to_flat_node_index(const PhaseSpaceHierarchy& hierarchy, const ConditioningPath& path);

} // namespace qmomkit
