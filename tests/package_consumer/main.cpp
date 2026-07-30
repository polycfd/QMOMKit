#include <qmomkit/qmomkit.hpp>

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

    return 0;
}
