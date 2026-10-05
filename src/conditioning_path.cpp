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

MultiIndex ConditioningPathView::phase_space_node_index(const PhaseSpaceHierarchy& hierarchy) const
{
    if (size() > hierarchy.number_of_dimensions())
    {
        throw std::invalid_argument("Conditioning path is deeper than the hierarchy.");
    }

    MultiIndex node_index(hierarchy.number_of_dimensions(), 0);
    for (std::size_t level = 0; level < size(); ++level)
    {
        const std::size_t phase_space_dimension = hierarchy.conditioning_order()[level];
        const std::size_t node = node_indices_[level];
        if (node >= hierarchy.node_count(phase_space_dimension))
        {
            throw std::out_of_range("Conditioning path node index is outside the hierarchy.");
        }
        node_index[phase_space_dimension] = node;
    }
    return node_index;
}

MultiIndex ConditioningPath::phase_space_node_index(const PhaseSpaceHierarchy& hierarchy) const
{
    return ConditioningPathView(*this).phase_space_node_index(hierarchy);
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
    // Discards the computed index; called only for its checks that this path is
    // no deeper than the hierarchy and that its node indices are within it.
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

std::size_t ConditioningPathView::node_index_at_level(const std::size_t conditioning_level) const
{
    if (conditioning_level >= levels_)
    {
        throw std::out_of_range("Conditioning path level is out of range.");
    }
    return node_indices_[conditioning_level];
}

bool operator==(const ConditioningPathView& lhs, const ConditioningPathView& rhs) noexcept
{
    return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

bool operator!=(const ConditioningPathView& lhs, const ConditioningPathView& rhs) noexcept
{
    return !(lhs == rhs);
}

ConditioningPathTable::ConditioningPathTable(const std::size_t levels, std::vector<std::size_t> node_indices)
    : node_indices_(std::move(node_indices)), levels_(levels)
{
    if (levels_ == 0)
    {
        if (!node_indices_.empty())
        {
            throw std::invalid_argument("A conditioning-path table without levels must be empty.");
        }
        return;
    }
    if (node_indices_.size() % levels_ != 0)
    {
        throw std::invalid_argument("Conditioning-path table storage must contain a whole number of equal-depth paths.");
    }
}

ConditioningPathTable::ConditioningPathTable(const std::vector<ConditioningPath>& paths)
{
    if (paths.empty())
    {
        return;
    }
    levels_ = paths.front().size();
    node_indices_.reserve(paths.size() * levels_);
    for (const ConditioningPath& path : paths)
    {
        if (path.size() != levels_)
        {
            throw std::invalid_argument("Conditioning-path table entries must all reach the same conditioning level.");
        }
        const std::vector<std::size_t>& entries = path.node_indices_by_level();
        node_indices_.insert(node_indices_.end(), entries.begin(), entries.end());
    }
}

ConditioningPathView ConditioningPathTable::at(const std::size_t path) const
{
    if (path >= size())
    {
        throw std::out_of_range("Conditioning-path table row is outside the table.");
    }
    return (*this)[path];
}

void ConditioningPathTable::push_back(const ConditioningPath& path)
{
    push_back(ConditioningPathView(path));
}

void ConditioningPathTable::push_back(const ConditioningPathView path)
{
    if (node_indices_.empty() && levels_ == 0)
    {
        levels_ = path.size();
    }
    if (path.size() != levels_)
    {
        throw std::invalid_argument("Conditioning-path table entries must all reach the same conditioning level.");
    }
    node_indices_.insert(node_indices_.end(), path.begin(), path.end());
}

void ConditioningPathTable::clear() noexcept
{
    node_indices_.clear();
}

void ConditioningPathTable::reserve(const std::size_t paths)
{
    node_indices_.reserve(paths * levels_);
}

std::vector<ConditioningPath> ConditioningPathTable::to_conditioning_paths() const
{
    std::vector<ConditioningPath> paths;
    paths.reserve(size());
    for (std::size_t path = 0; path < size(); ++path)
    {
        paths.push_back((*this)[path].to_conditioning_path());
    }
    return paths;
}

bool operator==(const ConditioningPathTable& lhs, const ConditioningPathTable& rhs) noexcept
{
    return lhs.levels() == rhs.levels() && lhs.node_indices() == rhs.node_indices();
}

bool operator!=(const ConditioningPathTable& lhs, const ConditioningPathTable& rhs) noexcept
{
    return !(lhs == rhs);
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

MultiIndex recursive_path_to_node_index(const PhaseSpaceHierarchy& hierarchy, const ConditioningPathView path)
{
    if (path.size() != hierarchy.number_of_dimensions())
    {
        throw std::invalid_argument("A recursive quadrature storage path must include every conditioning level.");
    }
    return path.phase_space_node_index(hierarchy);
}

} // namespace qmomkit
