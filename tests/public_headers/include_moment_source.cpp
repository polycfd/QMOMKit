#include <qmomkit/moment_source.hpp>
#include <type_traits>

static_assert(std::is_default_constructible_v<qmomkit::MomentBoundaryRate>);
static_assert(std::is_default_constructible_v<qmomkit::MomentBoundaryContribution>);
static_assert(std::is_default_constructible_v<qmomkit::CqmomMomentBoundaryContribution>);
