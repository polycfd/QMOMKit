#include <qmomkit/transport_integrator.hpp>

namespace
{

using QmomContinuousStage = qmomkit::QmomTransportResult (*)(const qmomkit::QmomState&, double, double, const qmomkit::PhaseSpaceVelocity&,
                                                             const qmomkit::QmomTransportOptions&);
using QmomDiscontinuousStage = qmomkit::QmomTransportResult (*)(const qmomkit::QmomState&, std::size_t, double, double, const qmomkit::MomentSource&,
                                                                const qmomkit::QmomTransportOptions&);
using QmomBoundaryDiscontinuousStage = qmomkit::QmomTransportResult (*)(const qmomkit::QmomState&, std::size_t, double, double,
                                                                        const qmomkit::MomentSource&, const qmomkit::MomentBoundaryContribution&,
                                                                        const qmomkit::QmomTransportOptions&);

[[maybe_unused]] QmomContinuousStage qmom_continuous_stage = qmomkit::advance_qmom_continuous;
[[maybe_unused]] QmomDiscontinuousStage qmom_discontinuous_stage = qmomkit::advance_qmom_discontinuous;
[[maybe_unused]] QmomBoundaryDiscontinuousStage qmom_boundary_discontinuous_stage = qmomkit::advance_qmom_discontinuous;

} // namespace
