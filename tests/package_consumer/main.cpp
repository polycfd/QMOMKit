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

    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>& coordinate)
    { return std::vector<double>{0.0 * coordinate.front()}; };
    const qmomkit::MomentSource source = [](double, const qmomkit::QmomSourceContext& state) { return std::vector<double>(state.size(), 0.0); };
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

    const qmomkit::BatchPhaseSpaceVelocity batch_velocity = [](std::size_t, double, const std::vector<double>& coordinate)
    { return std::vector<double>(coordinate.size(), 0.0); };
    const qmomkit::BatchMomentSource batch_source = [](std::size_t, double, const qmomkit::QmomSourceContext& state)
    { return std::vector<double>(state.size(), 0.0); };
    const std::vector<qmomkit::QmomPhaseSpaceEvolutionResult> qmom_batch_aware =
        qmomkit::advance_qmom_phase_space_batch_with_index(qmom_states, 0.0, 0.1, batch_velocity, batch_source, qmom_options);
    if (qmom_batch_aware.size() != 1 || !qmom_batch_aware.front().ok())
    {
        return 1;
    }

    const qmomkit::PhaseSpaceHierarchy cqmom_hierarchy({1, 1});
    const qmomkit::RawMomentTensor cqmom_dense_moments{{2, 2}, {1.0, 0.0, 0.0, 0.0}};

    qmomkit::CqmomDenseCase dense_case{cqmom_hierarchy, cqmom_dense_moments};
    const qmomkit::CqmomResult dense_static = qmomkit::compute_cqmom(dense_case);
    if (!dense_static.ok())
    {
        return 1;
    }
    const qmomkit::CqmomState dense_state = qmomkit::CqmomCase{std::move(dense_case)}.initial_state();
    if (dense_state.quadrature().total_nodes() != 1 || dense_state.quadrature().weights().size() != 1)
    {
        return 1;
    }

    qmomkit::CqmomLeanCase lean_case{qmomkit::extract_cqmom_moments(cqmom_dense_moments, qmomkit::create_cqmom_moment_layout(cqmom_hierarchy))};
    const qmomkit::CqmomResult lean_static = qmomkit::compute_cqmom(lean_case);
    if (!lean_static.ok())
    {
        return 1;
    }
    const qmomkit::CqmomState lean_state = qmomkit::CqmomCase{std::move(lean_case)}.initial_state();
    if (lean_state.quadrature().total_nodes() != 1 || lean_state.quadrature().weights().size() != 1)
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
