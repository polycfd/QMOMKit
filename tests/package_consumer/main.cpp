#include <algorithm>
#include <qmomkit/qmomkit.hpp>
#include <utility>
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

    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity([](double, const std::vector<double>& coordinate)
                                                                                    { return std::vector<double>{0.0 * coordinate.front()}; });
    const qmomkit::QmomMomentSource source = [](double, const qmomkit::QmomSourceContext& state)
    { return qmomkit::QmomMomentRate{std::vector<double>(state.size(), 0.0)}; };
    const qmomkit::QmomState qmom_state = qmomkit::QmomCase{moments, 2, {}}.initial_state();
    if (qmom_state.maximum_node_count() != 2 || qmom_state.active_node_count() != 2 || qmom_state.qmom_options().adaptive)
    {
        return 1;
    }
    const qmomkit::QmomPhaseSpaceEvolutionResult qmom_source = qmomkit::advance_qmom_moment_source(qmom_state, 0.0, 0.1, source);
    if (!qmom_source.ok() || qmom_source.state.maximum_node_count() != qmom_state.maximum_node_count())
    {
        return 1;
    }
    const std::vector<qmomkit::QmomState> qmom_states{qmom_state};
    qmomkit::QmomPhaseSpaceEvolutionOptions qmom_options;
    qmom_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> qmom_batch =
        qmomkit::advance_qmom_phase_space_batch(qmom_states, 0.0, 0.1, velocity, source, qmom_options);
    if (qmom_batch.size() != 1 || !qmom_batch.front().ok())
    {
        return 1;
    }

    const qmomkit::BatchPhaseSpaceVelocity batch_velocity = [](std::size_t, double, const std::vector<double>&, std::vector<double>& velocities,
                                                               std::size_t) { std::fill(velocities.begin(), velocities.end(), 0.0); };
    const qmomkit::BatchQmomMomentSource batch_source = [](std::size_t, double, const qmomkit::QmomSourceContext& state)
    { return qmomkit::QmomMomentRate{std::vector<double>(state.size(), 0.0)}; };
    const std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> qmom_batch_aware =
        qmomkit::advance_qmom_phase_space_batch_with_index(qmom_states, 0.0, 0.1, batch_velocity, batch_source, qmom_options);
    if (qmom_batch_aware.size() != 1 || !qmom_batch_aware.front().ok())
    {
        return 1;
    }

    const qmomkit::PhaseSpaceHierarchy cqmom_hierarchy({1, 1});
    const qmomkit::RawMomentTensor cqmom_moment_tensor{{2, 2}, {1.0, 0.0, 0.0, 0.0}};

    qmomkit::CqmomCase cqmom_case{qmomkit::extract_cqmom_moments(cqmom_moment_tensor, qmomkit::create_cqmom_moment_layout(cqmom_hierarchy))};
    const qmomkit::CqmomResult cqmom_result = qmomkit::compute_cqmom(cqmom_case);
    if (!cqmom_result.ok())
    {
        return 1;
    }
    const qmomkit::CqmomState cqmom_case_state = cqmom_case.initial_state();
    if (cqmom_case_state.quadrature().total_nodes() != 1)
    {
        return 1;
    }
    // Both overloads must store exactly the replacement moments, which differ
    // from the ones the case holds.
    const qmomkit::RawMomentTensor replacement_moment_tensor{{2, 2}, {2.0, 1.0, 0.5, 0.25}};
    const qmomkit::CqmomMomentData replacement_moments = qmomkit::extract_cqmom_moments(replacement_moment_tensor, cqmom_case.moments.layout);
    const auto tensor_replaced_cqmom_case_result = qmomkit::config::with_initial_moments(cqmom_case, replacement_moment_tensor);
    const auto moment_set_replaced_cqmom_case_result = qmomkit::config::with_initial_moments(cqmom_case, replacement_moments);
    if (replacement_moments.values == cqmom_case.moments.values || !tensor_replaced_cqmom_case_result.ok() ||
        tensor_replaced_cqmom_case_result.value().moments.values != replacement_moments.values || !moment_set_replaced_cqmom_case_result.ok() ||
        moment_set_replaced_cqmom_case_result.value().moments.values != replacement_moments.values)
    {
        return 1;
    }
    const qmomkit::CqmomState cqmom_prepared_case_state = qmomkit::CqmomPreparedCase{std::move(cqmom_case)}.initial_state();
    if (cqmom_prepared_case_state.quadrature().total_nodes() != 1 || cqmom_prepared_case_state.quadrature().weights().size() != 1)
    {
        return 1;
    }

    const auto base_cqmom_prepared_case_result = qmomkit::config::parse_cqmom_prepared(R"(
[phase_space]
node_counts = 1, 1

[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 0.0
M[1,0] = 0.0
M[1,1] = 0.0
)");
    if (!base_cqmom_prepared_case_result.ok())
    {
        return 1;
    }

    const auto replaced_cqmom_prepared_case_result =
        qmomkit::config::with_initial_moments(base_cqmom_prepared_case_result.value(), {{2, 2}, {1.0, 0.5, 1.0, 0.5}});
    if (!replaced_cqmom_prepared_case_result.ok())
    {
        return 1;
    }
    const qmomkit::CqmomState replaced_state = replaced_cqmom_prepared_case_result.value().initial_state();
    if (replaced_state.moment_values().empty())
    {
        return 1;
    }

    const std::vector<qmomkit::CqmomState> cqmom_states{replaced_state};
    qmomkit::CqmomPhaseSpaceEvolutionOptions cqmom_options;
    cqmom_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::PhaseSpaceVelocity cqmom_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& coordinate) { return std::vector<double>(coordinate.size(), 0.0); });
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
