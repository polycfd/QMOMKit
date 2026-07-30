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
using CqmomContinuousStage = qmomkit::CqmomTransportResult (*)(const qmomkit::CqmomState&, double, double, const qmomkit::PhaseSpaceVelocity&,
                                                               const qmomkit::CqmomTransportOptions&);
using CqmomDiscontinuousStage = qmomkit::CqmomTransportResult (*)(const qmomkit::CqmomState&, double, double, const qmomkit::CqmomMomentSource&,
                                                                  const qmomkit::CqmomTransportOptions&);
using CqmomBoundaryDiscontinuousStage = qmomkit::CqmomTransportResult (*)(const qmomkit::CqmomState&, double, double,
                                                                          const qmomkit::CqmomMomentSource&,
                                                                          const qmomkit::CqmomMomentBoundaryContribution&,
                                                                          const qmomkit::CqmomTransportOptions&);

[[maybe_unused]] QmomContinuousStage qmom_continuous_stage = qmomkit::advance_qmom_continuous;
[[maybe_unused]] QmomDiscontinuousStage qmom_discontinuous_stage = qmomkit::advance_qmom_discontinuous;
[[maybe_unused]] QmomBoundaryDiscontinuousStage qmom_boundary_discontinuous_stage = qmomkit::advance_qmom_discontinuous;
[[maybe_unused]] CqmomContinuousStage cqmom_continuous_stage = qmomkit::advance_cqmom_continuous;
[[maybe_unused]] CqmomDiscontinuousStage cqmom_discontinuous_stage = qmomkit::advance_cqmom_discontinuous;
[[maybe_unused]] CqmomBoundaryDiscontinuousStage cqmom_boundary_discontinuous_stage = qmomkit::advance_cqmom_discontinuous;

} // namespace
