#include <qmomkit/conditioning_path.hpp>

// Standard headers follow the public header, so they cannot supply an include it is missing.
#include <cstddef>
#include <initializer_list>
#include <type_traits>
#include <vector>

// The view is a non-owning pair of pointer and length, so passing one by value
// must stay as cheap as passing a reference to an owning path.
static_assert(std::is_trivially_copyable_v<qmomkit::ConditioningPathView>);
static_assert(std::is_nothrow_default_constructible_v<qmomkit::ConditioningPathView>);
static_assert(sizeof(qmomkit::ConditioningPathView) <= 2 * sizeof(void*));
static_assert(std::is_convertible_v<const qmomkit::ConditioningPath&, qmomkit::ConditioningPathView>);

static_assert(std::is_default_constructible_v<qmomkit::ConditioningPathTable>);
static_assert(std::is_copy_constructible_v<qmomkit::ConditioningPathTable>);
static_assert(std::is_nothrow_move_constructible_v<qmomkit::ConditioningPathTable>);

// A braced list of paths is copied into flat storage, so unlike a view it is
// safe to build a table from one.
static_assert(std::is_constructible_v<qmomkit::ConditioningPathTable, std::initializer_list<qmomkit::ConditioningPath>>);
static_assert(std::is_constructible_v<qmomkit::ConditioningPathTable, std::size_t, std::vector<std::size_t>>);
static_assert(std::is_constructible_v<qmomkit::ConditioningPathTable, const std::vector<qmomkit::ConditioningPath>&>);

using PathRow = qmomkit::ConditioningPathView (qmomkit::ConditioningPathTable::*)(std::size_t) const;
using PathStorage = const std::vector<std::size_t>& (qmomkit::ConditioningPathTable::*)() const noexcept;
using PathToOwning = std::vector<qmomkit::ConditioningPath> (qmomkit::ConditioningPathTable::*)() const;

static_assert(std::is_same_v<decltype(&qmomkit::ConditioningPathTable::at), PathRow>);
static_assert(std::is_same_v<decltype(&qmomkit::ConditioningPathTable::node_indices), PathStorage>);
static_assert(std::is_same_v<decltype(&qmomkit::ConditioningPathTable::to_conditioning_paths), PathToOwning>);

// Entries stay ordered by conditioning level, so the view exposes the same
// level-indexed accessor as an owning path.
using ViewLevel = std::size_t (qmomkit::ConditioningPathView::*)(std::size_t) const;
static_assert(std::is_same_v<decltype(&qmomkit::ConditioningPathView::node_index_at_level), ViewLevel>);
static_assert(
    std::is_same_v<decltype(&qmomkit::ConditioningPath::node_index_at_level), std::size_t (qmomkit::ConditioningPath::*)(std::size_t) const>);
