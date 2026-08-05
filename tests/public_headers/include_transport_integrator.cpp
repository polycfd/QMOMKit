#include <cstddef>
#include <functional>
#include <qmomkit/transport_integrator.hpp>
#include <type_traits>
#include <vector>

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

using QmomContinuousBatch = std::vector<qmomkit::QmomTransportResult> (*)(const std::vector<qmomkit::QmomState>&, double, double,
                                                                          const qmomkit::PhaseSpaceVelocity&, const qmomkit::QmomTransportOptions&);
using CqmomContinuousBatch = std::vector<qmomkit::CqmomTransportResult> (*)(const std::vector<qmomkit::CqmomState>&, double, double,
                                                                            const qmomkit::PhaseSpaceVelocity&,
                                                                            const qmomkit::CqmomTransportOptions&);
using QmomDiscontinuousBatch = std::vector<qmomkit::QmomTransportResult> (*)(const std::vector<qmomkit::QmomState>&, std::size_t, double, double,
                                                                             const qmomkit::MomentSource&, const qmomkit::QmomTransportOptions&);
using CqmomDiscontinuousBatch = std::vector<qmomkit::CqmomTransportResult> (*)(const std::vector<qmomkit::CqmomState>&, double, double,
                                                                               const qmomkit::CqmomMomentSource&,
                                                                               const qmomkit::CqmomTransportOptions&);
using QmomBoundaryDiscontinuousBatch = std::vector<qmomkit::QmomTransportResult> (*)(const std::vector<qmomkit::QmomState>&, std::size_t, double,
                                                                                     double, const qmomkit::MomentSource&,
                                                                                     const qmomkit::MomentBoundaryContribution&,
                                                                                     const qmomkit::QmomTransportOptions&);
using CqmomBoundaryDiscontinuousBatch = std::vector<qmomkit::CqmomTransportResult> (*)(const std::vector<qmomkit::CqmomState>&, double, double,
                                                                                       const qmomkit::CqmomMomentSource&,
                                                                                       const qmomkit::CqmomMomentBoundaryContribution&,
                                                                                       const qmomkit::CqmomTransportOptions&);
using QmomTransportBatch = std::vector<qmomkit::QmomTransportResult> (*)(const std::vector<qmomkit::QmomState>&, std::size_t, double, double,
                                                                         const qmomkit::PhaseSpaceVelocity&, const qmomkit::MomentSource&,
                                                                         const qmomkit::QmomTransportOptions&);
using CqmomTransportBatch = std::vector<qmomkit::CqmomTransportResult> (*)(const std::vector<qmomkit::CqmomState>&, double, double,
                                                                           const qmomkit::PhaseSpaceVelocity&, const qmomkit::CqmomMomentSource&,
                                                                           const qmomkit::CqmomTransportOptions&);
using QmomBoundaryTransportBatch = std::vector<qmomkit::QmomTransportResult> (*)(const std::vector<qmomkit::QmomState>&, std::size_t, double, double,
                                                                                 const qmomkit::PhaseSpaceVelocity&, const qmomkit::MomentSource&,
                                                                                 const qmomkit::MomentBoundaryContribution&,
                                                                                 const qmomkit::QmomTransportOptions&);
using CqmomBoundaryTransportBatch = std::vector<qmomkit::CqmomTransportResult> (*)(const std::vector<qmomkit::CqmomState>&, double, double,
                                                                                   const qmomkit::PhaseSpaceVelocity&,
                                                                                   const qmomkit::CqmomMomentSource&,
                                                                                   const qmomkit::CqmomMomentBoundaryContribution&,
                                                                                   const qmomkit::CqmomTransportOptions&);

using QmomBatchContinuousBatch = std::vector<qmomkit::QmomTransportResult> (*)(const std::vector<qmomkit::QmomState>&, double, double,
                                                                               const qmomkit::BatchPhaseSpaceVelocity&,
                                                                               const qmomkit::QmomTransportOptions&);
using CqmomBatchContinuousBatch = std::vector<qmomkit::CqmomTransportResult> (*)(const std::vector<qmomkit::CqmomState>&, double, double,
                                                                                 const qmomkit::BatchPhaseSpaceVelocity&,
                                                                                 const qmomkit::CqmomTransportOptions&);
using QmomBatchDiscontinuousBatch = std::vector<qmomkit::QmomTransportResult> (*)(const std::vector<qmomkit::QmomState>&, std::size_t, double, double,
                                                                                  const qmomkit::BatchMomentSource&,
                                                                                  const qmomkit::QmomTransportOptions&);
using CqmomBatchDiscontinuousBatch = std::vector<qmomkit::CqmomTransportResult> (*)(const std::vector<qmomkit::CqmomState>&, double, double,
                                                                                    const qmomkit::BatchCqmomMomentSource&,
                                                                                    const qmomkit::CqmomTransportOptions&);
using QmomBatchBoundaryDiscontinuousBatch = std::vector<qmomkit::QmomTransportResult> (*)(const std::vector<qmomkit::QmomState>&, std::size_t, double,
                                                                                          double, const qmomkit::BatchMomentSource&,
                                                                                          const qmomkit::BatchMomentBoundaryContribution&,
                                                                                          const qmomkit::QmomTransportOptions&);
using CqmomBatchBoundaryDiscontinuousBatch = std::vector<qmomkit::CqmomTransportResult> (*)(const std::vector<qmomkit::CqmomState>&, double, double,
                                                                                            const qmomkit::BatchCqmomMomentSource&,
                                                                                            const qmomkit::BatchCqmomMomentBoundaryContribution&,
                                                                                            const qmomkit::CqmomTransportOptions&);
using QmomBatchTransportBatch = std::vector<qmomkit::QmomTransportResult> (*)(const std::vector<qmomkit::QmomState>&, std::size_t, double, double,
                                                                              const qmomkit::BatchPhaseSpaceVelocity&,
                                                                              const qmomkit::BatchMomentSource&,
                                                                              const qmomkit::QmomTransportOptions&);
using CqmomBatchTransportBatch = std::vector<qmomkit::CqmomTransportResult> (*)(const std::vector<qmomkit::CqmomState>&, double, double,
                                                                                const qmomkit::BatchPhaseSpaceVelocity&,
                                                                                const qmomkit::BatchCqmomMomentSource&,
                                                                                const qmomkit::CqmomTransportOptions&);
using QmomBatchBoundaryTransportBatch = std::vector<qmomkit::QmomTransportResult> (*)(const std::vector<qmomkit::QmomState>&, std::size_t, double,
                                                                                      double, const qmomkit::BatchPhaseSpaceVelocity&,
                                                                                      const qmomkit::BatchMomentSource&,
                                                                                      const qmomkit::BatchMomentBoundaryContribution&,
                                                                                      const qmomkit::QmomTransportOptions&);
using CqmomBatchBoundaryTransportBatch = std::vector<qmomkit::CqmomTransportResult> (*)(const std::vector<qmomkit::CqmomState>&, double, double,
                                                                                        const qmomkit::BatchPhaseSpaceVelocity&,
                                                                                        const qmomkit::BatchCqmomMomentSource&,
                                                                                        const qmomkit::BatchCqmomMomentBoundaryContribution&,
                                                                                        const qmomkit::CqmomTransportOptions&);

static_assert(std::is_same_v<qmomkit::BatchPhaseSpaceVelocity, std::function<std::vector<double>(std::size_t, double, const std::vector<double>&)>>);
static_assert(std::is_same_v<qmomkit::BatchMomentSource, std::function<std::vector<double>(std::size_t, double, const std::vector<double>&)>>);
static_assert(std::is_same_v<qmomkit::BatchMomentBoundaryContribution,
                             std::function<qmomkit::MomentBoundaryRate(std::size_t, double, const std::vector<double>&)>>);
static_assert(std::is_same_v<qmomkit::BatchCqmomMomentSource,
                             std::function<qmomkit::CqmomMomentRate(std::size_t, double, const qmomkit::CqmomSourceContext&)>>);
static_assert(std::is_same_v<qmomkit::BatchCqmomMomentBoundaryContribution,
                             std::function<qmomkit::CqmomMomentRate(std::size_t, double, const qmomkit::CqmomSourceContext&)>>);

struct DualPhaseSpaceVelocity
{
    [[nodiscard]] std::vector<double> operator()(double, const std::vector<double>&) const
    {
        return {};
    }

    [[nodiscard]] std::vector<double> operator()(std::size_t, double, const std::vector<double>&) const
    {
        return {};
    }
};

struct DualCqmomMomentSource
{
    [[nodiscard]] qmomkit::CqmomMomentRate operator()(double, const qmomkit::CqmomSourceContext&) const
    {
        return qmomkit::CqmomMomentRate{std::vector<double>{}};
    }

    [[nodiscard]] qmomkit::CqmomMomentRate operator()(std::size_t, double, const qmomkit::CqmomSourceContext&) const
    {
        return qmomkit::CqmomMomentRate{std::vector<double>{}};
    }
};

[[maybe_unused]] void verify_batch_callback_overload_resolution()
{
    const std::vector<qmomkit::QmomState> qmom_states;
    const std::vector<qmomkit::CqmomState> cqmom_states;
    const qmomkit::QmomTransportOptions qmom_options;
    const qmomkit::CqmomTransportOptions cqmom_options;

    const auto generic_velocity = [](auto&&...) { return std::vector<double>{}; };
    (void)qmomkit::advance_qmom_continuous_batch(qmom_states, 0.0, 0.1, generic_velocity);
    (void)qmomkit::advance_qmom_continuous_batch_with_index(qmom_states, 0.0, 0.1, generic_velocity);
    (void)qmomkit::advance_qmom_continuous_batch(qmom_states, 0.0, 0.1, DualPhaseSpaceVelocity{});
    (void)qmomkit::advance_qmom_continuous_batch_with_index(qmom_states, 0.0, 0.1, DualPhaseSpaceVelocity{});
    (void)qmomkit::advance_qmom_continuous_batch(qmom_states, 0.0, 0.1, nullptr);
    (void)qmomkit::advance_qmom_continuous_batch_with_index(qmom_states, 0.0, 0.1, nullptr);

    const auto generic_qmom_source = [](auto&&...) { return std::vector<double>{}; };
    const auto generic_qmom_boundary = [](auto&&...) { return qmomkit::MomentBoundaryRate{std::vector<double>{}}; };
    (void)qmomkit::advance_qmom_discontinuous_batch(qmom_states, 1, 0.0, 0.1, generic_qmom_source);
    (void)qmomkit::advance_qmom_discontinuous_batch_with_index(qmom_states, 1, 0.0, 0.1, generic_qmom_source);
    (void)qmomkit::advance_qmom_discontinuous_batch(qmom_states, 1, 0.0, 0.1, generic_qmom_source, generic_qmom_boundary, qmom_options);
    (void)qmomkit::advance_qmom_discontinuous_batch_with_index(qmom_states, 1, 0.0, 0.1, generic_qmom_source, generic_qmom_boundary, qmom_options);
    (void)qmomkit::advance_qmom_discontinuous_batch(qmom_states, 1, 0.0, 0.1, nullptr);
    (void)qmomkit::advance_qmom_discontinuous_batch_with_index(qmom_states, 1, 0.0, 0.1, nullptr);
    (void)qmomkit::advance_qmom_discontinuous_batch(qmom_states, 1, 0.0, 0.1, generic_qmom_source, nullptr, qmom_options);
    (void)qmomkit::advance_qmom_discontinuous_batch_with_index(qmom_states, 1, 0.0, 0.1, generic_qmom_source, nullptr, qmom_options);

    const auto generic_cqmom_source = [](auto&&...) { return qmomkit::CqmomMomentRate{std::vector<double>{}}; };
    const auto generic_cqmom_boundary = [](auto&&...) { return qmomkit::CqmomMomentRate{std::vector<double>{}}; };
    (void)qmomkit::advance_cqmom_discontinuous_batch(cqmom_states, 0.0, 0.1, generic_cqmom_source);
    (void)qmomkit::advance_cqmom_discontinuous_batch_with_index(cqmom_states, 0.0, 0.1, generic_cqmom_source);
    (void)qmomkit::advance_cqmom_discontinuous_batch(cqmom_states, 0.0, 0.1, generic_cqmom_source, generic_cqmom_boundary, cqmom_options);
    (void)qmomkit::advance_cqmom_discontinuous_batch_with_index(cqmom_states, 0.0, 0.1, generic_cqmom_source, generic_cqmom_boundary, cqmom_options);
    (void)qmomkit::advance_cqmom_discontinuous_batch(cqmom_states, 0.0, 0.1, nullptr);
    (void)qmomkit::advance_cqmom_discontinuous_batch_with_index(cqmom_states, 0.0, 0.1, nullptr);
    (void)qmomkit::advance_cqmom_discontinuous_batch(cqmom_states, 0.0, 0.1, generic_cqmom_source, nullptr, cqmom_options);
    (void)qmomkit::advance_cqmom_discontinuous_batch_with_index(cqmom_states, 0.0, 0.1, generic_cqmom_source, nullptr, cqmom_options);
    (void)qmomkit::advance_cqmom_discontinuous_batch(cqmom_states, 0.0, 0.1, DualCqmomMomentSource{});
    (void)qmomkit::advance_cqmom_discontinuous_batch_with_index(cqmom_states, 0.0, 0.1, DualCqmomMomentSource{});
}

[[maybe_unused]] QmomContinuousStage qmom_continuous_stage = qmomkit::advance_qmom_continuous;
[[maybe_unused]] QmomDiscontinuousStage qmom_discontinuous_stage = qmomkit::advance_qmom_discontinuous;
[[maybe_unused]] QmomBoundaryDiscontinuousStage qmom_boundary_discontinuous_stage = qmomkit::advance_qmom_discontinuous;
[[maybe_unused]] CqmomContinuousStage cqmom_continuous_stage = qmomkit::advance_cqmom_continuous;
[[maybe_unused]] CqmomDiscontinuousStage cqmom_discontinuous_stage = qmomkit::advance_cqmom_discontinuous;
[[maybe_unused]] CqmomBoundaryDiscontinuousStage cqmom_boundary_discontinuous_stage = qmomkit::advance_cqmom_discontinuous;
[[maybe_unused]] QmomContinuousBatch qmom_continuous_batch = qmomkit::advance_qmom_continuous_batch;
[[maybe_unused]] CqmomContinuousBatch cqmom_continuous_batch = qmomkit::advance_cqmom_continuous_batch;
[[maybe_unused]] QmomDiscontinuousBatch qmom_discontinuous_batch = qmomkit::advance_qmom_discontinuous_batch;
[[maybe_unused]] CqmomDiscontinuousBatch cqmom_discontinuous_batch = qmomkit::advance_cqmom_discontinuous_batch;
[[maybe_unused]] QmomBoundaryDiscontinuousBatch qmom_boundary_discontinuous_batch = qmomkit::advance_qmom_discontinuous_batch;
[[maybe_unused]] CqmomBoundaryDiscontinuousBatch cqmom_boundary_discontinuous_batch = qmomkit::advance_cqmom_discontinuous_batch;
[[maybe_unused]] QmomTransportBatch qmom_transport_batch = qmomkit::advance_transport_qmom_batch;
[[maybe_unused]] CqmomTransportBatch cqmom_transport_batch = qmomkit::advance_transport_cqmom_batch;
[[maybe_unused]] QmomBoundaryTransportBatch qmom_boundary_transport_batch = qmomkit::advance_transport_qmom_batch;
[[maybe_unused]] CqmomBoundaryTransportBatch cqmom_boundary_transport_batch = qmomkit::advance_transport_cqmom_batch;
[[maybe_unused]] QmomBatchContinuousBatch qmom_batch_continuous_batch = qmomkit::advance_qmom_continuous_batch_with_index;
[[maybe_unused]] CqmomBatchContinuousBatch cqmom_batch_continuous_batch = qmomkit::advance_cqmom_continuous_batch_with_index;
[[maybe_unused]] QmomBatchDiscontinuousBatch qmom_batch_discontinuous_batch = qmomkit::advance_qmom_discontinuous_batch_with_index;
[[maybe_unused]] CqmomBatchDiscontinuousBatch cqmom_batch_discontinuous_batch = qmomkit::advance_cqmom_discontinuous_batch_with_index;
[[maybe_unused]] QmomBatchBoundaryDiscontinuousBatch qmom_batch_boundary_discontinuous_batch = qmomkit::advance_qmom_discontinuous_batch_with_index;
[[maybe_unused]] CqmomBatchBoundaryDiscontinuousBatch cqmom_batch_boundary_discontinuous_batch =
    qmomkit::advance_cqmom_discontinuous_batch_with_index;
[[maybe_unused]] QmomBatchTransportBatch qmom_batch_transport_batch = qmomkit::advance_transport_qmom_batch_with_index;
[[maybe_unused]] CqmomBatchTransportBatch cqmom_batch_transport_batch = qmomkit::advance_transport_cqmom_batch_with_index;
[[maybe_unused]] QmomBatchBoundaryTransportBatch qmom_batch_boundary_transport_batch = qmomkit::advance_transport_qmom_batch_with_index;
[[maybe_unused]] CqmomBatchBoundaryTransportBatch cqmom_batch_boundary_transport_batch = qmomkit::advance_transport_cqmom_batch_with_index;

} // namespace
