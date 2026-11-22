#include <cstddef>
#include <qmomkit/phase_space_hierarchy.hpp>
#include <type_traits>
#include <utility>
#include <vector>

// The family of each coordinate, with the settings of the generalized family,
// is classical by default, and both structs are aggregates.
static_assert(std::is_aggregate_v<qmomkit::CoordinateFamily>);
static_assert(std::is_aggregate_v<qmomkit::HierarchyFamilies>);
static_assert(std::is_same_v<decltype(qmomkit::CoordinateFamily::family), qmomkit::QmomFamily>);
static_assert(std::is_same_v<decltype(qmomkit::CoordinateFamily::generalized_options), qmomkit::GeneralizedOptions>);
static_assert(qmomkit::CoordinateFamily{}.family == qmomkit::QmomFamily::classical);
static_assert(std::is_same_v<decltype(qmomkit::HierarchyFamilies::coordinates), std::vector<qmomkit::CoordinateFamily>>);

// A hierarchy is built from node counts, a conditioning order and the families
// of its coordinates, and it returns the families and the maximum node counts.
static_assert(std::is_constructible_v<qmomkit::PhaseSpaceHierarchy, std::vector<std::size_t>, std::vector<std::size_t>, qmomkit::HierarchyFamilies>);
static_assert(std::is_same_v<decltype(std::declval<const qmomkit::PhaseSpaceHierarchy&>().families()), const qmomkit::HierarchyFamilies&>);
static_assert(std::is_same_v<decltype(std::declval<const qmomkit::PhaseSpaceHierarchy&>().maximum_node_counts()), const std::vector<std::size_t>&>);
