#include <qmomkit/qmomkit.hpp>
#include <vector>

namespace
{

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

[[maybe_unused]] CqmomAdvectionStage cqmom_advection_stage = qmomkit::advance_cqmom_advection;
[[maybe_unused]] CqmomMomentSourceStage cqmom_moment_source_stage = qmomkit::advance_cqmom_moment_source;
[[maybe_unused]] CqmomBoundaryMomentSourceStage cqmom_boundary_moment_source_stage = qmomkit::advance_cqmom_moment_source;

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
    const qmomkit::QmomAdvectionMomentRateResult advection_rates = qmomkit::evaluate_qmom_advection_moment_rates(moments.values(), 2, 0.0, velocity);
    if (!advection_rates.ok() || advection_rates.rates.size() != moments.values().size())
    {
        return 1;
    }

    const qmomkit::MomentSource source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };
    const qmomkit::QmomMomentSourceRateResult moment_source_rates = qmomkit::evaluate_qmom_moment_source_rates(0.0, moments.values(), source);
    if (!moment_source_rates.ok() || moment_source_rates.rates.size() != moments.values().size())
    {
        return 1;
    }

    const qmomkit::MomentBoundaryContribution boundary = [](double, const std::vector<double>& state)
    { return qmomkit::MomentBoundaryRate{std::vector<double>(state.size(), 0.0)}; };
    const qmomkit::QmomMomentSourceRateResult moment_source_rates_with_boundary =
        qmomkit::evaluate_qmom_moment_source_rates(0.0, moments.values(), source, boundary);
    if (!moment_source_rates_with_boundary.ok() || moment_source_rates_with_boundary.rates.size() != moments.values().size())
    {
        return 1;
    }

    const qmomkit::QmomState qmom_state = qmomkit::QmomCase{moments, 2, {}}.initial_state();
    const std::vector<qmomkit::QmomState> qmom_states{qmom_state};
    qmomkit::QmomPhaseSpaceEvolutionOptions qmom_options;
    qmom_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> qmom_batch =
        qmomkit::advance_qmom_phase_space_batch(qmom_states, 2, 0.0, 0.1, velocity, source, qmom_options);
    if (qmom_batch.size() != 1 || !qmom_batch.front().ok())
    {
        return 1;
    }

    const qmomkit::BatchPhaseSpaceVelocity batch_velocity = [](std::size_t, double, const std::vector<double>& coordinate)
    { return std::vector<double>(coordinate.size(), 0.0); };
    const qmomkit::BatchMomentSource batch_source = [](std::size_t, double, const std::vector<double>& state)
    { return std::vector<double>(state.size(), 0.0); };
    const std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> qmom_batch_aware =
        qmomkit::advance_qmom_phase_space_batch_with_index(qmom_states, 2, 0.0, 0.1, batch_velocity, batch_source, qmom_options);
    if (qmom_batch_aware.size() != 1 || !qmom_batch_aware.front().ok())
    {
        return 1;
    }

    const auto base_case = qmomkit::config::parse_cqmom(R"(
[phase_space]
node_counts = 1, 1

[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 0.0
M[1,0] = 0.0
M[1,1] = 0.0
)");
    if (!base_case.ok())
    {
        return 1;
    }

    const auto replaced_case = qmomkit::config::with_initial_moments(base_case.value(), {{2, 2}, {1.0, 0.5, 1.0, 0.5}});
    if (!replaced_case.ok())
    {
        return 1;
    }
    const qmomkit::CqmomState replaced_state = replaced_case.value().initial_state();
    if (replaced_state.moment_storage() != qmomkit::CqmomMomentStorage::lean || replaced_state.moment_values().empty())
    {
        return 1;
    }

    const std::vector<qmomkit::CqmomState> cqmom_states{replaced_state};
    qmomkit::CqmomPhaseSpaceEvolutionOptions cqmom_options;
    cqmom_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::PhaseSpaceVelocity cqmom_velocity = [](double, const std::vector<double>& coordinate)
    { return std::vector<double>(coordinate.size(), 0.0); };
    const std::vector<qmomkit::CqmomPhaseSpaceEvolutionResult> cqmom_batch =
        qmomkit::advance_cqmom_advection_batch(cqmom_states, 0.0, 0.1, cqmom_velocity, cqmom_options);
    if (cqmom_batch.size() != 1 || !cqmom_batch.front().ok())
    {
        return 1;
    }

    const std::vector<qmomkit::CqmomPhaseSpaceEvolutionResult> cqmom_batch_aware =
        qmomkit::advance_cqmom_advection_batch_with_index(cqmom_states, 0.0, 0.1, batch_velocity, cqmom_options);
    if (cqmom_batch_aware.size() != 1 || !cqmom_batch_aware.front().ok())
    {
        return 1;
    }

    return 0;
}
