#pragma once

#include <cstddef>
#include <vector>

namespace qmomkit
{

using MultiIndex = std::vector<std::size_t>;

/**
 * Non-owning view of one multi-index stored inside a `MultiIndexTable`.
 *
 * The view is valid only while the owning table is alive and unmodified. It
 * offers the subset of `MultiIndex` operations that reading code needs, so a
 * `const MultiIndex&` can usually be replaced by a `MultiIndexView` without
 * touching the expressions that consume it.
 */
class MultiIndexView
{
public:
    MultiIndexView() = default;
    MultiIndexView(const std::size_t* entries, std::size_t rank) noexcept : entries_(entries), rank_(rank) {}

    [[nodiscard]] std::size_t size() const noexcept
    {
        return rank_;
    }
    [[nodiscard]] bool empty() const noexcept
    {
        return rank_ == 0;
    }
    [[nodiscard]] std::size_t operator[](const std::size_t axis) const noexcept
    {
        return entries_[axis];
    }
    [[nodiscard]] const std::size_t* begin() const noexcept
    {
        return entries_;
    }
    [[nodiscard]] const std::size_t* end() const noexcept
    {
        return entries_ + rank_;
    }
    [[nodiscard]] const std::size_t* data() const noexcept
    {
        return entries_;
    }

    /** Copy the viewed entries into an owning multi-index. */
    [[nodiscard]] MultiIndex to_multi_index() const
    {
        return MultiIndex(entries_, entries_ + rank_);
    }

private:
    const std::size_t* entries_ = nullptr;
    std::size_t rank_ = 0;
};

[[nodiscard]] bool operator==(const MultiIndexView& lhs, const MultiIndexView& rhs) noexcept;
[[nodiscard]] bool operator!=(const MultiIndexView& lhs, const MultiIndexView& rhs) noexcept;
[[nodiscard]] bool operator==(const MultiIndexView& lhs, const MultiIndex& rhs) noexcept;
[[nodiscard]] bool operator==(const MultiIndex& lhs, const MultiIndexView& rhs) noexcept;
[[nodiscard]] bool operator!=(const MultiIndexView& lhs, const MultiIndex& rhs) noexcept;
[[nodiscard]] bool operator!=(const MultiIndex& lhs, const MultiIndexView& rhs) noexcept;

/**
 * Contiguous table of equal-rank multi-indices.
 *
 * Entries are stored in one flat array of `size() * rank()` values, row-major,
 * so a table of `n` indices costs one allocation instead of `n + 1`, and a
 * traversal walks memory in order. This matters because moment layouts and
 * terminal-path tables are copied on every reconstruction, where a
 * `std::vector<MultiIndex>` allocates once per row.
 *
 * All rows have the same rank, which is the phase-space dimension count for
 * moment layouts and the conditioning-level count for terminal paths. An empty
 * table has rank zero.
 */
class MultiIndexTable
{
public:
    MultiIndexTable() = default;

    /** Build from flat row-major storage; `values.size()` must be a multiple of `rank`. */
    MultiIndexTable(std::size_t rank, std::vector<std::size_t> values);

    /** Build from individual indices, which must all have the same rank. */
    explicit MultiIndexTable(const std::vector<MultiIndex>& indices);

    [[nodiscard]] std::size_t size() const noexcept
    {
        return rank_ == 0 ? 0 : values_.size() / rank_;
    }
    [[nodiscard]] std::size_t rank() const noexcept
    {
        return rank_;
    }
    [[nodiscard]] bool empty() const noexcept
    {
        return values_.empty();
    }
    [[nodiscard]] MultiIndexView operator[](const std::size_t row) const noexcept
    {
        return MultiIndexView(values_.data() + row * rank_, rank_);
    }
    [[nodiscard]] MultiIndexView at(std::size_t row) const;

    /** Flat row-major storage, `size() * rank()` values. */
    [[nodiscard]] const std::vector<std::size_t>& values() const noexcept
    {
        return values_;
    }

    /** Append one index, which must match the table rank unless the table is empty. */
    void push_back(const MultiIndex& index);
    void push_back(MultiIndexView index);

    void reserve(std::size_t rows);
    void clear() noexcept;

    /** Copy the table into individual owning indices. */
    [[nodiscard]] std::vector<MultiIndex> to_multi_indices() const;

private:
    std::vector<std::size_t> values_;
    std::size_t rank_ = 0;
};

[[nodiscard]] bool operator==(const MultiIndexTable& lhs, const MultiIndexTable& rhs) noexcept;
[[nodiscard]] bool operator!=(const MultiIndexTable& lhs, const MultiIndexTable& rhs) noexcept;

/**
 * Return the row-major product of extents.
 *
 * The last index is the fastest-varying index. All extents must be positive.
 */
[[nodiscard]] std::size_t product_extents(const std::vector<std::size_t>& extents);

/**
 * Convert a row-major multi-index to a flat index.
 *
 * The supplied multi-index and extents use the same axis ordering. The last
 * axis is fastest varying, matching C/C++ row-major storage.
 */
[[nodiscard]] std::size_t flatten_multi_index(const MultiIndex& multi_index, const std::vector<std::size_t>& extents);

/**
 * Convert a row-major flat index to a multi-index.
 */
[[nodiscard]] MultiIndex unflatten_multi_index(std::size_t flat_index, const std::vector<std::size_t>& extents);

/**
 * Enumerate all row-major multi-indices for the supplied extents.
 */
[[nodiscard]] std::vector<MultiIndex> enumerate_multi_indices(const std::vector<std::size_t>& extents);

} // namespace qmomkit
