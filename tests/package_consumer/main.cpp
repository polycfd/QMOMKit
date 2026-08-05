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

    const qmomkit::MomentBoundaryContribution boundary = [](double, const std::vector<double>& state)
    { return qmomkit::MomentBoundaryRate{std::vector<double>(state.size(), 0.0)}; };
    const qmomkit::QmomDiscontinuousMomentRateResult discontinuous_rates_with_boundary =
        qmomkit::evaluate_qmom_discontinuous_moment_rates(0.0, moments.values(), source, boundary);
    if (!discontinuous_rates_with_boundary.ok() || discontinuous_rates_with_boundary.rates.size() != moments.values().size())
    {
        return 1;
    }

    const qmomkit::QmomState qmom_state = qmomkit::QmomCase{moments, 2, {}}.initial_state();
    const std::vector<qmomkit::QmomState> qmom_states{qmom_state};
    qmomkit::QmomTransportOptions qmom_options;
    qmom_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const std::vector<qmomkit::QmomTransportResult> qmom_batch =
        qmomkit::advance_transport_qmom_batch(qmom_states, 2, 0.0, 0.1, velocity, source, qmom_options);
    if (qmom_batch.size() != 1 || !qmom_batch.front().ok())
    {
        return 1;
    }

    const qmomkit::BatchPhaseSpaceVelocity batch_velocity = [](std::size_t, double, const std::vector<double>& coordinate)
    { return std::vector<double>(coordinate.size(), 0.0); };
    const qmomkit::BatchMomentSource batch_source = [](std::size_t, double, const std::vector<double>& state)
    { return std::vector<double>(state.size(), 0.0); };
    const std::vector<qmomkit::QmomTransportResult> qmom_batch_aware =
        qmomkit::advance_transport_qmom_batch_with_index(qmom_states, 2, 0.0, 0.1, batch_velocity, batch_source, qmom_options);
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
    qmomkit::CqmomTransportOptions cqmom_options;
    cqmom_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::PhaseSpaceVelocity cqmom_velocity = [](double, const std::vector<double>& coordinate)
    { return std::vector<double>(coordinate.size(), 0.0); };
    const std::vector<qmomkit::CqmomTransportResult> cqmom_batch =
        qmomkit::advance_cqmom_continuous_batch(cqmom_states, 0.0, 0.1, cqmom_velocity, cqmom_options);
    if (cqmom_batch.size() != 1 || !cqmom_batch.front().ok())
    {
        return 1;
    }

    const std::vector<qmomkit::CqmomTransportResult> cqmom_batch_aware =
        qmomkit::advance_cqmom_continuous_batch_with_index(cqmom_states, 0.0, 0.1, batch_velocity, cqmom_options);
    if (cqmom_batch_aware.size() != 1 || !cqmom_batch_aware.front().ok())
    {
        return 1;
    }

    return 0;
}
