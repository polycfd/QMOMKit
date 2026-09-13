#pragma once

#include <cstddef>
#include <initializer_list>
#include <iterator>
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

    /**
     * View an owning multi-index.
     *
     * Implicit so that a function taking a `MultiIndexView` accepts both a
     * stored index and a table row. As with any view, the referenced index must
     * outlive the view; binding one to a temporary dangles.
     */
    MultiIndexView(const MultiIndex& index) noexcept : entries_(index.data()), rank_(index.size()) {}

    /**
     * View a braced list of components, as in `moment({1, 0, 2})`.
     *
     * Intended for the argument position, where the list outlives the call. A
     * view stored from a braced list dangles immediately, exactly as one stored
     * from a temporary index does.
     */
    MultiIndexView(std::initializer_list<std::size_t> entries) noexcept : entries_(entries.begin()), rank_(entries.size()) {}

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
 * moment layouts and the conditioning-level count for terminal paths. A
 * default-constructed table has rank zero and adopts the rank of the first
 * index appended to it. `clear()` removes the rows but keeps the rank, so a
 * cleared table can be refilled with indices of the same shape; two empty
 * tables of different rank are therefore not equal.
 */
class MultiIndexTable
{
public:
    /**
     * Row iterator yielding `MultiIndexView` by value.
     *
     * Because it hands back a view rather than a reference into a container of
     * indices, it satisfies the input-iterator requirements rather than the
     * forward-iterator ones. That is enough for range-based `for` and for
     * `std::find`, which is how the moment layouts are traversed.
     */
    class const_iterator
    {
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = MultiIndexView;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using reference = MultiIndexView;

        const_iterator() = default;
        const_iterator(const std::size_t* entries, std::size_t rank) noexcept : entries_(entries), rank_(rank) {}

        [[nodiscard]] reference operator*() const noexcept
        {
            return MultiIndexView(entries_, rank_);
        }
        const_iterator& operator++() noexcept
        {
            entries_ += rank_;
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
            return entries_ == other.entries_;
        }
        [[nodiscard]] bool operator!=(const const_iterator& other) const noexcept
        {
            return entries_ != other.entries_;
        }
        [[nodiscard]] difference_type operator-(const const_iterator& other) const noexcept
        {
            return rank_ == 0 ? 0 : (entries_ - other.entries_) / static_cast<difference_type>(rank_);
        }

    private:
        const std::size_t* entries_ = nullptr;
        std::size_t rank_ = 0;
    };

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

    [[nodiscard]] const_iterator begin() const noexcept
    {
        return const_iterator(values_.data(), rank_);
    }
    [[nodiscard]] const_iterator end() const noexcept
    {
        return const_iterator(values_.data() + values_.size(), rank_);
    }

    /** Flat row-major storage, `size() * rank()` values. */
    [[nodiscard]] const std::vector<std::size_t>& values() const noexcept
    {
        return values_;
    }

    /** Append one index, which must match the table rank unless the table is empty. */
    void push_back(const MultiIndex& index);
    void push_back(MultiIndexView index);

    /** Remove the last index; the table must not be empty. */
    void pop_back();

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
[[nodiscard]] std::size_t flatten_multi_index(MultiIndexView multi_index, const std::vector<std::size_t>& extents);

/**
 * Convert a row-major flat index to a multi-index.
 */
[[nodiscard]] MultiIndex unflatten_multi_index(std::size_t flat_index, const std::vector<std::size_t>& extents);

/**
 * Enumerate all row-major multi-indices for the supplied extents.
 */
[[nodiscard]] std::vector<MultiIndex> enumerate_multi_indices(const std::vector<std::size_t>& extents);

} // namespace qmomkit
