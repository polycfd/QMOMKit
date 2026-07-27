#include <qmomkit/moment_rates.hpp>
#include <type_traits>

static_assert(std::is_default_constructible_v<qmomkit::QmomContinuousMomentRateResult>);
static_assert(std::is_default_constructible_v<qmomkit::CqmomContinuousMomentRateResult>);
static_assert(std::is_default_constructible_v<qmomkit::QmomDiscontinuousMomentRateResult>);
static_assert(std::is_default_constructible_v<qmomkit::CqmomDiscontinuousMomentRateResult>);
