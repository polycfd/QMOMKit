#include <qmomkit/moment_rates.hpp>
#include <type_traits>
#include <vector>

static_assert(std::is_default_constructible_v<qmomkit::QmomContinuousMomentRateResult>);
static_assert(std::is_default_constructible_v<qmomkit::CqmomContinuousMomentRateResult>);
static_assert(std::is_default_constructible_v<qmomkit::QmomDiscontinuousMomentRateResult>);
static_assert(std::is_default_constructible_v<qmomkit::CqmomDiscontinuousMomentRateResult>);

namespace
{

[[maybe_unused]] qmomkit::QmomDiscontinuousMomentRateResult evaluate_qmom_rates_with_boundary()
{
    const qmomkit::MomentSource source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };
    const qmomkit::MomentBoundaryContribution boundary = [](double, const std::vector<double>& state)
    { return qmomkit::MomentBoundaryRate{std::vector<double>(state.size(), 0.0)}; };
    return qmomkit::evaluate_qmom_discontinuous_moment_rates(0.0, {1.0}, source, boundary);
}

} // namespace
