#include "qmomkit/conditioning_path.hpp"

#include <algorithm>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace qmomkit
{

ConditioningPath::ConditioningPath(std::vector<std::size_t> node_indices_by_level) : node_indices_by_level_(std::move(node_indices_by_level)) {}

bool ConditioningPath::empty() const noexcept
{
    return node_indices_by_level_.empty();
}

std::size_t ConditioningPath::size() const noexcept
{
    return node_indices_by_level_.size();
}

const std::vector<std::size_t>& ConditioningPath::node_indices_by_level() const noexcept
{
    return node_indices_by_level_;
}

std::size_t ConditioningPath::node_index_at_level(const std::size_t conditioning_level) const
{
    if (conditioning_level >= node_indices_by_level_.size())
    {
        throw std::out_of_range("Conditioning path level is out of range.");
    }
    return node_indices_by_level_[conditioning_level];
}

ConditioningPath ConditioningPath::with_child(const std::size_t child_node_index) const
{
    std::vector<std::size_t> child_path = node_indices_by_level_;
    child_path.push_back(child_node_index);
    return ConditioningPath(std::move(child_path));
}

MultiIndex ConditioningPath::phase_space_node_index(const PhaseSpaceHierarchy& hierarchy) const
{
    if (node_indices_by_level_.size() > hierarchy.number_of_dimensions())
    {
        throw std::invalid_argument("Conditioning path is deeper than the hierarchy.");
    }

    MultiIndex node_index(hierarchy.number_of_dimensions(), 0);
    for (std::size_t level = 0; level < node_indices_by_level_.size(); ++level)
    {
        const std::size_t phase_space_dimension = hierarchy.conditioning_order()[level];
        const std::size_t node = node_indices_by_level_[level];
        if (node >= hierarchy.node_count(phase_space_dimension))
        {
            throw std::out_of_range("Conditioning path node index is outside the hierarchy.");
        }
        node_index[phase_space_dimension] = node;
    }
    return node_index;
}

std::string ConditioningPath::to_string() const
{
    std::ostringstream stream;
    stream << "[";
    for (std::size_t level = 0; level < node_indices_by_level_.size(); ++level)
    {
        if (level != 0)
        {
            stream << ", ";
        }
        stream << "level " << level << ": node " << node_indices_by_level_[level];
    }
    stream << "]";
    return stream.str();
}

std::string ConditioningPath::describe(const PhaseSpaceHierarchy& hierarchy) const
{
    (void)phase_space_node_index(hierarchy);

    std::ostringstream stream;
    stream << "[";
    for (std::size_t level = 0; level < node_indices_by_level_.size(); ++level)
    {
        if (level != 0)
        {
            stream << ", ";
        }
        stream << "level " << level << " / dim " << hierarchy.conditioning_order()[level] << " / node " << node_indices_by_level_[level];
    }
    stream << "]";
    return stream.str();
}

bool operator==(const ConditioningPath& lhs, const ConditioningPath& rhs) noexcept
{
    return lhs.node_indices_by_level() == rhs.node_indices_by_level();
}

bool operator!=(const ConditioningPath& lhs, const ConditioningPath& rhs) noexcept
{
    return !(lhs == rhs);
}

bool operator<(const ConditioningPath& lhs, const ConditioningPath& rhs) noexcept
{
    return std::lexicographical_compare(lhs.node_indices_by_level().begin(), lhs.node_indices_by_level().end(), rhs.node_indices_by_level().begin(),
                                        rhs.node_indices_by_level().end());
}

MultiIndex recursive_path_to_node_index(const PhaseSpaceHierarchy& hierarchy, const ConditioningPath& path)
{
    if (path.size() != hierarchy.number_of_dimensions())
    {
        throw std::invalid_argument("A recursive quadrature storage path must include every conditioning level.");
    }
    return path.phase_space_node_index(hierarchy);
}

std::size_t recursive_path_to_flat_node_index(const PhaseSpaceHierarchy& hierarchy, const ConditioningPath& path)
{
    return hierarchy.flatten_node_index(recursive_path_to_node_index(hierarchy, path));
}

} // namespace qmomkit
