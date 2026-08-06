#include <qmomkit/moment_rates.hpp>
#include <type_traits>
#include <vector>

static_assert(std::is_default_constructible_v<qmomkit::QmomAdvectionMomentRateResult>);
static_assert(std::is_default_constructible_v<qmomkit::CqmomAdvectionMomentRateResult>);
static_assert(std::is_default_constructible_v<qmomkit::QmomMomentSourceRateResult>);
static_assert(std::is_default_constructible_v<qmomkit::CqmomMomentSourceRateResult>);

namespace
{

[[maybe_unused]] qmomkit::QmomMomentSourceRateResult evaluate_qmom_rates_with_boundary()
{
    const qmomkit::MomentSource source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };
    const qmomkit::MomentBoundaryContribution boundary = [](double, const std::vector<double>& state)
    { return qmomkit::MomentBoundaryRate{std::vector<double>(state.size(), 0.0)}; };
    return qmomkit::evaluate_qmom_moment_source_rates(0.0, {1.0}, source, boundary);
}

} // namespace
