#include <qmomkit/qmom.hpp>
#include <type_traits>
#include <utility>

static_assert(std::is_same_v<decltype(std::declval<const qmomkit::QmomState&>().maximum_node_count()), std::size_t>);
static_assert(std::is_same_v<decltype(std::declval<const qmomkit::QmomState&>().qmom_options()), const qmomkit::QmomOptions&>);
