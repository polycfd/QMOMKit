#include <cstddef>
#include <qmomkit/multi_index.hpp>
#include <type_traits>
#include <vector>

static_assert(std::is_same_v<qmomkit::MultiIndex, std::vector<std::size_t>>);

// The view is a non-owning pair of pointer and length, so passing one by value
// must stay as cheap as passing a reference to an owning index.
static_assert(std::is_trivially_copyable_v<qmomkit::MultiIndexView>);
static_assert(std::is_nothrow_default_constructible_v<qmomkit::MultiIndexView>);
static_assert(sizeof(qmomkit::MultiIndexView) <= 2 * sizeof(void*));

static_assert(std::is_default_constructible_v<qmomkit::MultiIndexTable>);
static_assert(std::is_copy_constructible_v<qmomkit::MultiIndexTable>);
static_assert(std::is_nothrow_move_constructible_v<qmomkit::MultiIndexTable>);

using TableFromFlat = qmomkit::MultiIndexTable (*)(std::size_t, std::vector<std::size_t>);
using TableRow = qmomkit::MultiIndexView (qmomkit::MultiIndexTable::*)(std::size_t) const;
using TableValues = const std::vector<std::size_t>& (qmomkit::MultiIndexTable::*)() const noexcept;

static_assert(std::is_same_v<decltype(&qmomkit::MultiIndexTable::at), TableRow>);
static_assert(std::is_same_v<decltype(&qmomkit::MultiIndexTable::values), TableValues>);

[[maybe_unused]] constexpr bool table_is_constructible_from_flat_storage =
    std::is_constructible_v<qmomkit::MultiIndexTable, std::size_t, std::vector<std::size_t>>;
static_assert(table_is_constructible_from_flat_storage);
static_assert(std::is_constructible_v<qmomkit::MultiIndexTable, const std::vector<qmomkit::MultiIndex>&>);

// Index-reading free functions accept both an owning index and a table row.
using FlattenIndex = std::size_t (*)(qmomkit::MultiIndexView, const std::vector<std::size_t>&);
static_assert(std::is_same_v<decltype(&qmomkit::flatten_multi_index), FlattenIndex>);

// A view is constructible from an owning index and from a braced list, which is
// what keeps read-only call sites unchanged.
static_assert(std::is_convertible_v<const qmomkit::MultiIndex&, qmomkit::MultiIndexView>);
static_assert(std::is_constructible_v<qmomkit::MultiIndexView, std::initializer_list<std::size_t>>);
