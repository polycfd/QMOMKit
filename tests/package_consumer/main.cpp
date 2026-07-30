#include <qmomkit/qmomkit.hpp>
#include <vector>

namespace
{

using CqmomContinuousStage = qmomkit::CqmomTransportResult (*)(const qmomkit::CqmomState&, double, double, const qmomkit::PhaseSpaceVelocity&,
                                                               const qmomkit::CqmomTransportOptions&);
using CqmomDiscontinuousStage = qmomkit::CqmomTransportResult (*)(const qmomkit::CqmomState&, double, double, const qmomkit::CqmomMomentSource&,
                                                                  const qmomkit::CqmomTransportOptions&);
using CqmomBoundaryDiscontinuousStage = qmomkit::CqmomTransportResult (*)(const qmomkit::CqmomState&, double, double,
                                                                          const qmomkit::CqmomMomentSource&,
                                                                          const qmomkit::CqmomMomentBoundaryContribution&,
                                                                          const qmomkit::CqmomTransportOptions&);

[[maybe_unused]] CqmomContinuousStage cqmom_continuous_stage = qmomkit::advance_cqmom_continuous;
[[maybe_unused]] CqmomDiscontinuousStage cqmom_discontinuous_stage = qmomkit::advance_cqmom_discontinuous;
[[maybe_unused]] CqmomBoundaryDiscontinuousStage cqmom_boundary_discontinuous_stage = qmomkit::advance_cqmom_discontinuous;

} // namespace

int main()
{
    const qmomkit::MomentSet moments({1.0, 0.8, 2.8, 4.4});
    const qmomkit::QmomResult result = qmomkit::compute_qmom(moments, 2);

    if (!result.ok())
    {
        return 1;
    }
    if (result.weights.size() != 2 || result.abscissas.size() != 2)
    {
        return 1;
    }

    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>& coordinate)
    { return std::vector<double>{0.0 * coordinate.front()}; };
    const qmomkit::QmomContinuousMomentRateResult continuous_rates =
        qmomkit::evaluate_qmom_continuous_moment_rates(moments.values(), 2, 0.0, velocity);
    if (!continuous_rates.ok() || continuous_rates.rates.size() != moments.values().size())
    {
        return 1;
    }

    const qmomkit::MomentSource source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };
    const qmomkit::QmomDiscontinuousMomentRateResult discontinuous_rates =
        qmomkit::evaluate_qmom_discontinuous_moment_rates(0.0, moments.values(), source);
    if (!discontinuous_rates.ok() || discontinuous_rates.rates.size() != moments.values().size())
    {
        return 1;
    }

    return 0;
}
