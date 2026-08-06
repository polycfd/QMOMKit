#include <cstddef>
#include <functional>
#include <qmomkit/phase_space_evolution.hpp>
#include <type_traits>
#include <vector>

namespace
{

using QmomAdvectionStage = qmomkit::QmomPhaseSpaceEvolutionResult (*)(const qmomkit::QmomState&, double, double, const qmomkit::PhaseSpaceVelocity&,
                                                                      const qmomkit::QmomPhaseSpaceEvolutionOptions&);
using QmomMomentSourceStage = qmomkit::QmomPhaseSpaceEvolutionResult (*)(const qmomkit::QmomState&, std::size_t, double, double,
                                                                         const qmomkit::MomentSource&,
                                                                         const qmomkit::QmomPhaseSpaceEvolutionOptions&);
using QmomBoundaryMomentSourceStage = qmomkit::QmomPhaseSpaceEvolutionResult (*)(const qmomkit::QmomState&, std::size_t, double, double,
                                                                                 const qmomkit::MomentSource&,
                                                                                 const qmomkit::MomentBoundaryContribution&,
                                                                                 const qmomkit::QmomPhaseSpaceEvolutionOptions&);
using CqmomAdvectionStage = qmomkit::CqmomPhaseSpaceEvolutionResult (*)(const qmomkit::CqmomState&, double, double,
                                                                        const qmomkit::PhaseSpaceVelocity&,
                                                                        const qmomkit::CqmomPhaseSpaceEvolutionOptions&);
using CqmomMomentSourceStage = qmomkit::CqmomPhaseSpaceEvolutionResult (*)(const qmomkit::CqmomState&, double, double,
                                                                           const qmomkit::CqmomMomentSource&,
                                                                           const qmomkit::CqmomPhaseSpaceEvolutionOptions&);
using CqmomBoundaryMomentSourceStage = qmomkit::CqmomPhaseSpaceEvolutionResult (*)(const qmomkit::CqmomState&, double, double,
                                                                                   const qmomkit::CqmomMomentSource&,
                                                                                   const qmomkit::CqmomMomentBoundaryContribution&,
                                                                                   const qmomkit::CqmomPhaseSpaceEvolutionOptions&);

using QmomAdvectionBatch = std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> (*)(const std::vector<qmomkit::QmomState>&, double, double,
                                                                                   const qmomkit::PhaseSpaceVelocity&,
                                                                                   const qmomkit::QmomPhaseSpaceEvolutionOptions&);
using CqmomAdvectionBatch = std::vector<qmomkit::CqmomPhaseSpaceEvolutionResult> (*)(const std::vector<qmomkit::CqmomState>&, double, double,
                                                                                     const qmomkit::PhaseSpaceVelocity&,
                                                                                     const qmomkit::CqmomPhaseSpaceEvolutionOptions&);
using QmomMomentSourceBatch = std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> (*)(const std::vector<qmomkit::QmomState>&, std::size_t, double,
                                                                                      double, const qmomkit::MomentSource&,
                                                                                      const qmomkit::QmomPhaseSpaceEvolutionOptions&);
using CqmomMomentSourceBatch = std::vector<qmomkit::CqmomPhaseSpaceEvolutionResult> (*)(const std::vector<qmomkit::CqmomState>&, double, double,
                                                                                        const qmomkit::CqmomMomentSource&,
                                                                                        const qmomkit::CqmomPhaseSpaceEvolutionOptions&);
using QmomBoundaryMomentSourceBatch = std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> (*)(const std::vector<qmomkit::QmomState>&, std::size_t,
                                                                                              double, double, const qmomkit::MomentSource&,
                                                                                              const qmomkit::MomentBoundaryContribution&,
                                                                                              const qmomkit::QmomPhaseSpaceEvolutionOptions&);
using CqmomBoundaryMomentSourceBatch = std::vector<qmomkit::CqmomPhaseSpaceEvolutionResult> (*)(const std::vector<qmomkit::CqmomState>&, double,
                                                                                                double, const qmomkit::CqmomMomentSource&,
                                                                                                const qmomkit::CqmomMomentBoundaryContribution&,
                                                                                                const qmomkit::CqmomPhaseSpaceEvolutionOptions&);
using QmomPhaseSpaceEvolutionBatch = std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> (*)(const std::vector<qmomkit::QmomState>&, std::size_t,
                                                                                             double, double, const qmomkit::PhaseSpaceVelocity&,
                                                                                             const qmomkit::MomentSource&,
                                                                                             const qmomkit::QmomPhaseSpaceEvolutionOptions&);
using CqmomPhaseSpaceEvolutionBatch = std::vector<qmomkit::CqmomPhaseSpaceEvolutionResult> (*)(const std::vector<qmomkit::CqmomState>&, double,
                                                                                               double, const qmomkit::PhaseSpaceVelocity&,
                                                                                               const qmomkit::CqmomMomentSource&,
                                                                                               const qmomkit::CqmomPhaseSpaceEvolutionOptions&);
using QmomBoundaryPhaseSpaceEvolutionBatch = std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> (*)(
    const std::vector<qmomkit::QmomState>&, std::size_t, double, double, const qmomkit::PhaseSpaceVelocity&, const qmomkit::MomentSource&,
    const qmomkit::MomentBoundaryContribution&, const qmomkit::QmomPhaseSpaceEvolutionOptions&);
using CqmomBoundaryPhaseSpaceEvolutionBatch = std::vector<qmomkit::CqmomPhaseSpaceEvolutionResult> (*)(
    const std::vector<qmomkit::CqmomState>&, double, double, const qmomkit::PhaseSpaceVelocity&, const qmomkit::CqmomMomentSource&,
    const qmomkit::CqmomMomentBoundaryContribution&, const qmomkit::CqmomPhaseSpaceEvolutionOptions&);

using QmomBatchAdvectionBatch = std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> (*)(const std::vector<qmomkit::QmomState>&, double, double,
                                                                                        const qmomkit::BatchPhaseSpaceVelocity&,
                                                                                        const qmomkit::QmomPhaseSpaceEvolutionOptions&);
using CqmomBatchAdvectionBatch = std::vector<qmomkit::CqmomPhaseSpaceEvolutionResult> (*)(const std::vector<qmomkit::CqmomState>&, double, double,
                                                                                          const qmomkit::BatchPhaseSpaceVelocity&,
                                                                                          const qmomkit::CqmomPhaseSpaceEvolutionOptions&);
using QmomBatchMomentSourceBatch = std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> (*)(const std::vector<qmomkit::QmomState>&, std::size_t,
                                                                                           double, double, const qmomkit::BatchMomentSource&,
                                                                                           const qmomkit::QmomPhaseSpaceEvolutionOptions&);
using CqmomBatchMomentSourceBatch = std::vector<qmomkit::CqmomPhaseSpaceEvolutionResult> (*)(const std::vector<qmomkit::CqmomState>&, double, double,
                                                                                             const qmomkit::BatchCqmomMomentSource&,
                                                                                             const qmomkit::CqmomPhaseSpaceEvolutionOptions&);
using QmomBatchBoundaryMomentSourceBatch = std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> (*)(const std::vector<qmomkit::QmomState>&,
                                                                                                   std::size_t, double, double,
                                                                                                   const qmomkit::BatchMomentSource&,
                                                                                                   const qmomkit::BatchMomentBoundaryContribution&,
                                                                                                   const qmomkit::QmomPhaseSpaceEvolutionOptions&);
using CqmomBatchBoundaryMomentSourceBatch = std::vector<qmomkit::CqmomPhaseSpaceEvolutionResult> (*)(
    const std::vector<qmomkit::CqmomState>&, double, double, const qmomkit::BatchCqmomMomentSource&,
    const qmomkit::BatchCqmomMomentBoundaryContribution&, const qmomkit::CqmomPhaseSpaceEvolutionOptions&);
using QmomBatchPhaseSpaceEvolutionBatch = std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> (*)(const std::vector<qmomkit::QmomState>&, std::size_t,
                                                                                                  double, double,
                                                                                                  const qmomkit::BatchPhaseSpaceVelocity&,
                                                                                                  const qmomkit::BatchMomentSource&,
                                                                                                  const qmomkit::QmomPhaseSpaceEvolutionOptions&);
using CqmomBatchPhaseSpaceEvolutionBatch = std::vector<qmomkit::CqmomPhaseSpaceEvolutionResult> (*)(const std::vector<qmomkit::CqmomState>&, double,
                                                                                                    double, const qmomkit::BatchPhaseSpaceVelocity&,
                                                                                                    const qmomkit::BatchCqmomMomentSource&,
                                                                                                    const qmomkit::CqmomPhaseSpaceEvolutionOptions&);
using QmomBatchBoundaryPhaseSpaceEvolutionBatch = std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> (*)(
    const std::vector<qmomkit::QmomState>&, std::size_t, double, double, const qmomkit::BatchPhaseSpaceVelocity&, const qmomkit::BatchMomentSource&,
    const qmomkit::BatchMomentBoundaryContribution&, const qmomkit::QmomPhaseSpaceEvolutionOptions&);
using CqmomBatchBoundaryPhaseSpaceEvolutionBatch = std::vector<qmomkit::CqmomPhaseSpaceEvolutionResult> (*)(
    const std::vector<qmomkit::CqmomState>&, double, double, const qmomkit::BatchPhaseSpaceVelocity&, const qmomkit::BatchCqmomMomentSource&,
    const qmomkit::BatchCqmomMomentBoundaryContribution&, const qmomkit::CqmomPhaseSpaceEvolutionOptions&);

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
    const qmomkit::QmomPhaseSpaceEvolutionOptions qmom_options;
    const qmomkit::CqmomPhaseSpaceEvolutionOptions cqmom_options;

    const auto generic_velocity = [](auto&&...) { return std::vector<double>{}; };
    (void)qmomkit::advance_qmom_advection_batch(qmom_states, 0.0, 0.1, generic_velocity);
    (void)qmomkit::advance_qmom_advection_batch_with_index(qmom_states, 0.0, 0.1, generic_velocity);
    (void)qmomkit::advance_qmom_advection_batch(qmom_states, 0.0, 0.1, DualPhaseSpaceVelocity{});
    (void)qmomkit::advance_qmom_advection_batch_with_index(qmom_states, 0.0, 0.1, DualPhaseSpaceVelocity{});
    (void)qmomkit::advance_qmom_advection_batch(qmom_states, 0.0, 0.1, nullptr);
    (void)qmomkit::advance_qmom_advection_batch_with_index(qmom_states, 0.0, 0.1, nullptr);

    const auto generic_qmom_source = [](auto&&...) { return std::vector<double>{}; };
    const auto generic_qmom_boundary = [](auto&&...) { return qmomkit::MomentBoundaryRate{std::vector<double>{}}; };
    (void)qmomkit::advance_qmom_moment_source_batch(qmom_states, 1, 0.0, 0.1, generic_qmom_source);
    (void)qmomkit::advance_qmom_moment_source_batch_with_index(qmom_states, 1, 0.0, 0.1, generic_qmom_source);
    (void)qmomkit::advance_qmom_moment_source_batch(qmom_states, 1, 0.0, 0.1, generic_qmom_source, generic_qmom_boundary, qmom_options);
    (void)qmomkit::advance_qmom_moment_source_batch_with_index(qmom_states, 1, 0.0, 0.1, generic_qmom_source, generic_qmom_boundary, qmom_options);
    (void)qmomkit::advance_qmom_moment_source_batch(qmom_states, 1, 0.0, 0.1, nullptr);
    (void)qmomkit::advance_qmom_moment_source_batch_with_index(qmom_states, 1, 0.0, 0.1, nullptr);
    (void)qmomkit::advance_qmom_moment_source_batch(qmom_states, 1, 0.0, 0.1, generic_qmom_source, nullptr, qmom_options);
    (void)qmomkit::advance_qmom_moment_source_batch_with_index(qmom_states, 1, 0.0, 0.1, generic_qmom_source, nullptr, qmom_options);

    const auto generic_cqmom_source = [](auto&&...) { return qmomkit::CqmomMomentRate{std::vector<double>{}}; };
    const auto generic_cqmom_boundary = [](auto&&...) { return qmomkit::CqmomMomentRate{std::vector<double>{}}; };
    (void)qmomkit::advance_cqmom_moment_source_batch(cqmom_states, 0.0, 0.1, generic_cqmom_source);
    (void)qmomkit::advance_cqmom_moment_source_batch_with_index(cqmom_states, 0.0, 0.1, generic_cqmom_source);
    (void)qmomkit::advance_cqmom_moment_source_batch(cqmom_states, 0.0, 0.1, generic_cqmom_source, generic_cqmom_boundary, cqmom_options);
    (void)qmomkit::advance_cqmom_moment_source_batch_with_index(cqmom_states, 0.0, 0.1, generic_cqmom_source, generic_cqmom_boundary, cqmom_options);
    (void)qmomkit::advance_cqmom_moment_source_batch(cqmom_states, 0.0, 0.1, nullptr);
    (void)qmomkit::advance_cqmom_moment_source_batch_with_index(cqmom_states, 0.0, 0.1, nullptr);
    (void)qmomkit::advance_cqmom_moment_source_batch(cqmom_states, 0.0, 0.1, generic_cqmom_source, nullptr, cqmom_options);
    (void)qmomkit::advance_cqmom_moment_source_batch_with_index(cqmom_states, 0.0, 0.1, generic_cqmom_source, nullptr, cqmom_options);
    (void)qmomkit::advance_cqmom_moment_source_batch(cqmom_states, 0.0, 0.1, DualCqmomMomentSource{});
    (void)qmomkit::advance_cqmom_moment_source_batch_with_index(cqmom_states, 0.0, 0.1, DualCqmomMomentSource{});
}

[[maybe_unused]] QmomAdvectionStage qmom_advection_stage = qmomkit::advance_qmom_advection;
[[maybe_unused]] QmomMomentSourceStage qmom_moment_source_stage = qmomkit::advance_qmom_moment_source;
[[maybe_unused]] QmomBoundaryMomentSourceStage qmom_boundary_moment_source_stage = qmomkit::advance_qmom_moment_source;
[[maybe_unused]] CqmomAdvectionStage cqmom_advection_stage = qmomkit::advance_cqmom_advection;
[[maybe_unused]] CqmomMomentSourceStage cqmom_moment_source_stage = qmomkit::advance_cqmom_moment_source;
[[maybe_unused]] CqmomBoundaryMomentSourceStage cqmom_boundary_moment_source_stage = qmomkit::advance_cqmom_moment_source;
[[maybe_unused]] QmomAdvectionBatch qmom_advection_batch = qmomkit::advance_qmom_advection_batch;
[[maybe_unused]] CqmomAdvectionBatch cqmom_advection_batch = qmomkit::advance_cqmom_advection_batch;
[[maybe_unused]] QmomMomentSourceBatch qmom_moment_source_batch = qmomkit::advance_qmom_moment_source_batch;
[[maybe_unused]] CqmomMomentSourceBatch cqmom_moment_source_batch = qmomkit::advance_cqmom_moment_source_batch;
[[maybe_unused]] QmomBoundaryMomentSourceBatch qmom_boundary_moment_source_batch = qmomkit::advance_qmom_moment_source_batch;
[[maybe_unused]] CqmomBoundaryMomentSourceBatch cqmom_boundary_moment_source_batch = qmomkit::advance_cqmom_moment_source_batch;
[[maybe_unused]] QmomPhaseSpaceEvolutionBatch qmom_phase_space_evolution_batch = qmomkit::advance_qmom_phase_space_batch;
[[maybe_unused]] CqmomPhaseSpaceEvolutionBatch cqmom_phase_space_evolution_batch = qmomkit::advance_cqmom_phase_space_batch;
[[maybe_unused]] QmomBoundaryPhaseSpaceEvolutionBatch qmom_boundary_phase_space_evolution_batch = qmomkit::advance_qmom_phase_space_batch;
[[maybe_unused]] CqmomBoundaryPhaseSpaceEvolutionBatch cqmom_boundary_phase_space_evolution_batch = qmomkit::advance_cqmom_phase_space_batch;
[[maybe_unused]] QmomBatchAdvectionBatch qmom_batch_advection_batch = qmomkit::advance_qmom_advection_batch_with_index;
[[maybe_unused]] CqmomBatchAdvectionBatch cqmom_batch_advection_batch = qmomkit::advance_cqmom_advection_batch_with_index;
[[maybe_unused]] QmomBatchMomentSourceBatch qmom_batch_moment_source_batch = qmomkit::advance_qmom_moment_source_batch_with_index;
[[maybe_unused]] CqmomBatchMomentSourceBatch cqmom_batch_moment_source_batch = qmomkit::advance_cqmom_moment_source_batch_with_index;
[[maybe_unused]] QmomBatchBoundaryMomentSourceBatch qmom_batch_boundary_moment_source_batch = qmomkit::advance_qmom_moment_source_batch_with_index;
[[maybe_unused]] CqmomBatchBoundaryMomentSourceBatch cqmom_batch_boundary_moment_source_batch = qmomkit::advance_cqmom_moment_source_batch_with_index;
[[maybe_unused]] QmomBatchPhaseSpaceEvolutionBatch qmom_batch_phase_space_evolution_batch = qmomkit::advance_qmom_phase_space_batch_with_index;
[[maybe_unused]] CqmomBatchPhaseSpaceEvolutionBatch cqmom_batch_phase_space_evolution_batch = qmomkit::advance_cqmom_phase_space_batch_with_index;
[[maybe_unused]] QmomBatchBoundaryPhaseSpaceEvolutionBatch qmom_batch_boundary_phase_space_evolution_batch =
    qmomkit::advance_qmom_phase_space_batch_with_index;
[[maybe_unused]] CqmomBatchBoundaryPhaseSpaceEvolutionBatch cqmom_batch_boundary_phase_space_evolution_batch =
    qmomkit::advance_cqmom_phase_space_batch_with_index;

} // namespace
