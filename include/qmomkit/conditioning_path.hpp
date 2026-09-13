#pragma once

#include "qmomkit/multi_index.hpp"
#include "qmomkit/phase_space_hierarchy.hpp"

#include <cstddef>
#include <initializer_list>
#include <iterator>
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

/**
 * Non-owning view of one path stored inside a `ConditioningPathTable`.
 *
 * Entries are ordered by conditioning level, exactly as in `ConditioningPath`.
 * The view is valid only while the owning table is alive and unmodified.
 */
class ConditioningPathView
{
public:
    ConditioningPathView() = default;
    ConditioningPathView(const std::size_t* node_indices, std::size_t levels) noexcept : node_indices_(node_indices), levels_(levels) {}

    /**
     * View an owning path.
     *
     * Implicit so that a function taking a view accepts both a stored path and
     * a table row. As with any view, the referenced path must outlive it.
     */
    ConditioningPathView(const ConditioningPath& path) noexcept
        : node_indices_(path.node_indices_by_level().data()), levels_(path.node_indices_by_level().size())
    {
    }

    [[nodiscard]] bool empty() const noexcept
    {
        return levels_ == 0;
    }
    [[nodiscard]] std::size_t size() const noexcept
    {
        return levels_;
    }
    [[nodiscard]] std::size_t node_index_at_level(std::size_t conditioning_level) const;
    [[nodiscard]] MultiIndex phase_space_node_index(const PhaseSpaceHierarchy& hierarchy) const;
    [[nodiscard]] const std::size_t* begin() const noexcept
    {
        return node_indices_;
    }
    [[nodiscard]] const std::size_t* end() const noexcept
    {
        return node_indices_ + levels_;
    }
    [[nodiscard]] const std::size_t* data() const noexcept
    {
        return node_indices_;
    }

    /** Copy the viewed entries into an owning path. */
    [[nodiscard]] ConditioningPath to_conditioning_path() const
    {
        return ConditioningPath(std::vector<std::size_t>(node_indices_, node_indices_ + levels_));
    }

private:
    const std::size_t* node_indices_ = nullptr;
    std::size_t levels_ = 0;
};

[[nodiscard]] bool operator==(const ConditioningPathView& lhs, const ConditioningPathView& rhs) noexcept;
[[nodiscard]] bool operator!=(const ConditioningPathView& lhs, const ConditioningPathView& rhs) noexcept;

/**
 * Contiguous table of terminal conditioning paths.
 *
 * Every terminal path reaches the deepest conditioning level, so all rows have
 * the same length even when branch node counts differ. The entries are stored
 * in one flat array of `size() * levels()` values, so a table of `n` paths
 * costs one allocation instead of `n + 1` and is copied in one allocation.
 * Quadratures carry these tables through node advection unchanged, where a
 * `std::vector<ConditioningPath>` allocated once per path on every copy.
 *
 * A default-constructed table has depth zero and adopts the depth of the first
 * path appended to it. `clear()` removes the paths but keeps the depth, so two
 * empty tables of different depth are not equal.
 */
class ConditioningPathTable
{
public:
    /** Row iterator yielding `ConditioningPathView` by value. */
    class const_iterator
    {
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = ConditioningPathView;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using reference = ConditioningPathView;

        const_iterator() = default;
        const_iterator(const std::size_t* node_indices, std::size_t levels) noexcept : node_indices_(node_indices), levels_(levels) {}

        [[nodiscard]] reference operator*() const noexcept
        {
            return ConditioningPathView(node_indices_, levels_);
        }
        const_iterator& operator++() noexcept
        {
            node_indices_ += levels_;
            return *this;
        }
        const_iterator operator++(int) noexcept
        {
            const_iterator previous = *this;
            ++(*this);
            return previous;
        }
        [[nodiscard]] bool operator==(const const_iterator& other) const noexcept
        {
            return node_indices_ == other.node_indices_;
        }
        [[nodiscard]] bool operator!=(const const_iterator& other) const noexcept
        {
            return node_indices_ != other.node_indices_;
        }

    private:
        const std::size_t* node_indices_ = nullptr;
        std::size_t levels_ = 0;
    };

    ConditioningPathTable() = default;

    /** Build from flat storage; `node_indices.size()` must be a multiple of `levels`. */
    ConditioningPathTable(std::size_t levels, std::vector<std::size_t> node_indices);

    /** Build from individual paths, which must all reach the same level. */
    explicit ConditioningPathTable(const std::vector<ConditioningPath>& paths);

    /** Build from a braced list of paths; the entries are copied into flat storage. */
    ConditioningPathTable(std::initializer_list<ConditioningPath> paths) : ConditioningPathTable(std::vector<ConditioningPath>(paths)) {}

    [[nodiscard]] std::size_t size() const noexcept
    {
        return levels_ == 0 ? 0 : node_indices_.size() / levels_;
    }
    /** Number of conditioning levels in every path, zero for a default-constructed table. */
    [[nodiscard]] std::size_t levels() const noexcept
    {
        return levels_;
    }
    [[nodiscard]] bool empty() const noexcept
    {
        return node_indices_.empty();
    }
    [[nodiscard]] ConditioningPathView operator[](const std::size_t path) const noexcept
    {
        return ConditioningPathView(node_indices_.data() + path * levels_, levels_);
    }
    [[nodiscard]] ConditioningPathView at(std::size_t path) const;

    [[nodiscard]] const_iterator begin() const noexcept
    {
        return const_iterator(node_indices_.data(), levels_);
    }
    [[nodiscard]] const_iterator end() const noexcept
    {
        return const_iterator(node_indices_.data() + node_indices_.size(), levels_);
    }

    /** Flat storage, `size() * levels()` values ordered by conditioning level within each path. */
    [[nodiscard]] const std::vector<std::size_t>& node_indices() const noexcept
    {
        return node_indices_;
    }

    /** Append one path, which must match the table depth unless the table is empty. */
    void push_back(const ConditioningPath& path);
    void push_back(ConditioningPathView path);

    void reserve(std::size_t paths);

    /** Remove every path; the conditioning depth is retained so the table can be refilled. */
    void clear() noexcept;

    /** Copy the table into individual owning paths. */
    [[nodiscard]] std::vector<ConditioningPath> to_conditioning_paths() const;

private:
    std::vector<std::size_t> node_indices_;
    std::size_t levels_ = 0;
};

[[nodiscard]] bool operator==(const ConditioningPathTable& lhs, const ConditioningPathTable& rhs) noexcept;
[[nodiscard]] bool operator!=(const ConditioningPathTable& lhs, const ConditioningPathTable& rhs) noexcept;

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
[[nodiscard]] MultiIndex recursive_path_to_node_index(const PhaseSpaceHierarchy& hierarchy, ConditioningPathView path);

} // namespace qmomkit
