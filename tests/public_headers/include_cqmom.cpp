#include <qmomkit/cqmom.hpp>
#include <type_traits>

static_assert(std::is_same_v<qmomkit::AbscissaBounds, qmomkit::CqmomAbscissaBounds>);
