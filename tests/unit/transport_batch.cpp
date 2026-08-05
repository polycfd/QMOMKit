#include "qmomkit/qmomkit.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace
{

void assert_near(const double actual, const double expected, const double tolerance = 1.0e-12)
{
    assert(std::abs(actual - expected) <= tolerance);
}

void assert_vector_near(const std::vector<double>& actual, const std::vector<double>& expected, const double tolerance = 1.0e-12)
{
    assert(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        assert_near(actual[i], expected[i], tolerance);
    }
}

void assert_node_diagnostics_equal(const qmomkit::NodeAdvectionDiagnostics& actual, const qmomkit::NodeAdvectionDiagnostics& expected)
{
    assert(actual.scheme == expected.scheme);
    assert(actual.step_count == expected.step_count);
    assert(actual.rejected_step_count == expected.rejected_step_count);
    assert(actual.rhs_evaluations == expected.rhs_evaluations);
    assert_near(actual.final_time, expected.final_time);
    assert_near(actual.last_step, expected.last_step);
    assert(actual.estimated_error == expected.estimated_error);
    assert(actual.adaptive_status == expected.adaptive_status);
    assert(actual.message == expected.message);
}

void assert_source_diagnostics_equal(const qmomkit::MomentSourceUpdateDiagnostics& actual, const qmomkit::MomentSourceUpdateDiagnostics& expected)
{
    assert(actual.ssp_rk3_status == expected.ssp_rk3_status);
    assert(actual.accepted == expected.accepted);
    assert(actual.rejected == expected.rejected);
    assert(actual.retries == expected.retries);
    assert_near(actual.time_step_used, expected.time_step_used);
    assert(actual.error_estimate == expected.error_estimate);
    assert(actual.message == expected.message);
}

void assert_cqmom_diagnostic_equal(const qmomkit::CqmomDiagnostics& actual, const qmomkit::CqmomDiagnostics& expected)
{
    assert(actual.conditioning_level == expected.conditioning_level);
    assert(actual.parent_path == expected.parent_path);
    assert(actual.phase_space_dimension == expected.phase_space_dimension);
    assert(actual.requested_node_count == expected.requested_node_count);
    assert(actual.actual_node_count == expected.actual_node_count);
    assert(actual.adaptive_branch_counts_enabled == expected.adaptive_branch_counts_enabled);
    assert(actual.adaptive_node_count_reduced == expected.adaptive_node_count_reduced);
    assert(actual.qmom_status == expected.qmom_status);
    assert(actual.conditional_moment_status == expected.conditional_moment_status);
    assert(actual.conditional_moment_diagnostics.residual_norm == expected.conditional_moment_diagnostics.residual_norm);
    assert(actual.conditional_moment_diagnostics.condition_estimate == expected.conditional_moment_diagnostics.condition_estimate);
    assert(actual.conditional_moment_diagnostics.conditioning_warning == expected.conditional_moment_diagnostics.conditioning_warning);
    assert(actual.conditional_vandermonde_requested_solver == expected.conditional_vandermonde_requested_solver);
    assert(actual.conditional_vandermonde_selected_solver == expected.conditional_vandermonde_selected_solver);
    assert(actual.dense_vandermonde_fallback_occurred == expected.dense_vandermonde_fallback_occurred);
    assert(actual.message == expected.message);
    assert(actual.bounded_support_checked == expected.bounded_support_checked);
    assert(actual.bounded_support_interval_present == expected.bounded_support_interval_present);
    assert_near(actual.bounded_support_interval.lower, expected.bounded_support_interval.lower);
    assert_near(actual.bounded_support_interval.upper, expected.bounded_support_interval.upper);
    assert(actual.bounded_previous_active_node_count == expected.bounded_previous_active_node_count);
    assert(actual.bounded_candidate_node_count == expected.bounded_candidate_node_count);
    assert(actual.bounded_candidate_attempt_type == expected.bounded_candidate_attempt_type);
    assert_vector_near(actual.bounded_candidate_moments, expected.bounded_candidate_moments);
    assert(actual.bounded_realizability_status == expected.bounded_realizability_status);
    assert_near(actual.bounded_effective_tolerance, expected.bounded_effective_tolerance);
    assert(actual.bounded_failed_matrix == expected.bounded_failed_matrix);
    assert_near(actual.bounded_minimum_margin, expected.bounded_minimum_margin);
    assert(actual.bounded_accepted_fallback_node_count == expected.bounded_accepted_fallback_node_count);
    assert_vector_near(actual.bounded_returned_nodes, expected.bounded_returned_nodes);
    assert_vector_near(actual.bounded_returned_weights, expected.bounded_returned_weights);
    assert(actual.bounded_support_determinants.size() == expected.bounded_support_determinants.size());
    for (std::size_t i = 0; i < actual.bounded_support_determinants.size(); ++i)
    {
        assert(actual.bounded_support_determinants[i].matrix_family == expected.bounded_support_determinants[i].matrix_family);
        assert(actual.bounded_support_determinants[i].matrix_order == expected.bounded_support_determinants[i].matrix_order);
        assert_near(actual.bounded_support_determinants[i].determinant, expected.bounded_support_determinants[i].determinant);
    }
}

void assert_cqmom_diagnostics_equal(const std::vector<qmomkit::CqmomDiagnostics>& actual, const std::vector<qmomkit::CqmomDiagnostics>& expected)
{
    assert(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        assert_cqmom_diagnostic_equal(actual[i], expected[i]);
    }
}

void assert_transport_diagnostics_equal(const qmomkit::TransportDiagnostics& actual, const qmomkit::TransportDiagnostics& expected)
{
    assert(actual.scheme == expected.scheme);
    assert(actual.failed_stage == expected.failed_stage);
    assert(actual.cqmom_moment_storage == expected.cqmom_moment_storage);
    assert(actual.cqmom_diagnostic_level == expected.cqmom_diagnostic_level);
    assert(actual.cqmom_vandermonde_solver == expected.cqmom_vandermonde_solver);
    assert(actual.post_source_qmom_status == expected.post_source_qmom_status);
    assert(actual.initial_cqmom_status == expected.initial_cqmom_status);
    assert(actual.post_source_cqmom_status == expected.post_source_cqmom_status);
    assert_cqmom_diagnostics_equal(actual.initial_cqmom_diagnostics, expected.initial_cqmom_diagnostics);
    assert_cqmom_diagnostics_equal(actual.post_source_cqmom_diagnostics, expected.post_source_cqmom_diagnostics);
    assert(actual.first_continuous_status == expected.first_continuous_status);
    assert(actual.second_continuous_status == expected.second_continuous_status);
    assert_node_diagnostics_equal(actual.first_continuous_diagnostics, expected.first_continuous_diagnostics);
    assert_node_diagnostics_equal(actual.second_continuous_diagnostics, expected.second_continuous_diagnostics);
    assert(actual.first_reconstruction_diagnostics.success == expected.first_reconstruction_diagnostics.success);
    assert(actual.first_reconstruction_diagnostics.moment_count == expected.first_reconstruction_diagnostics.moment_count);
    assert(actual.first_reconstruction_diagnostics.message == expected.first_reconstruction_diagnostics.message);
    assert(actual.final_reconstruction_diagnostics.success == expected.final_reconstruction_diagnostics.success);
    assert(actual.final_reconstruction_diagnostics.moment_count == expected.final_reconstruction_diagnostics.moment_count);
    assert(actual.final_reconstruction_diagnostics.message == expected.final_reconstruction_diagnostics.message);
    assert(actual.source_status == expected.source_status);
    assert_source_diagnostics_equal(actual.source_diagnostics, expected.source_diagnostics);
    assert(actual.cqmom_node_count_transitions.size() == expected.cqmom_node_count_transitions.size());
    for (std::size_t i = 0; i < actual.cqmom_node_count_transitions.size(); ++i)
    {
        const qmomkit::CqmomNodeCountTransition& actual_transition = actual.cqmom_node_count_transitions[i];
        const qmomkit::CqmomNodeCountTransition& expected_transition = expected.cqmom_node_count_transitions[i];
        assert(actual_transition.stage == expected_transition.stage);
        assert(actual_transition.conditioning_level == expected_transition.conditioning_level);
        assert(actual_transition.parent_path == expected_transition.parent_path);
        assert(actual_transition.phase_space_dimension == expected_transition.phase_space_dimension);
        assert(actual_transition.configured_maximum_node_count == expected_transition.configured_maximum_node_count);
        assert(actual_transition.active_node_count_before == expected_transition.active_node_count_before);
        assert(actual_transition.active_node_count_after == expected_transition.active_node_count_after);
        assert(actual_transition.kind == expected_transition.kind);
        assert(actual_transition.rank_reduced == expected_transition.rank_reduced);
        assert(actual_transition.rank_grew == expected_transition.rank_grew);
    }
    assert(actual.total_node_advection_steps == expected.total_node_advection_steps);
    assert(actual.total_node_advection_rejected_steps == expected.total_node_advection_rejected_steps);
    assert(actual.total_node_advection_rhs_evaluations == expected.total_node_advection_rhs_evaluations);
    assert(actual.source_retries == expected.source_retries);
    assert(actual.message == expected.message);
}

void assert_qmom_state_equal(const qmomkit::QmomState& actual, const qmomkit::QmomState& expected)
{
    assert_vector_near(actual.moments(), expected.moments());
    assert(actual.active_node_count() == expected.active_node_count());
    assert(actual.quadrature().status == expected.quadrature().status);
    assert(actual.quadrature().message == expected.quadrature().message);
    assert_vector_near(actual.quadrature().abscissas, expected.quadrature().abscissas);
    assert_vector_near(actual.quadrature().weights, expected.quadrature().weights);
}

void assert_qmom_result_equal(const qmomkit::QmomTransportResult& actual, const qmomkit::QmomTransportResult& expected)
{
    assert_qmom_state_equal(actual.state, expected.state);
    assert(actual.status == expected.status);
    assert_near(actual.final_time, expected.final_time);
    assert(actual.message == expected.message);
    assert_transport_diagnostics_equal(actual.diagnostics, expected.diagnostics);
}

void assert_active_counts_equal(const std::vector<qmomkit::CqmomActiveBranchCount>& actual,
                                const std::vector<qmomkit::CqmomActiveBranchCount>& expected)
{
    assert(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        assert(actual[i].conditioning_level == expected[i].conditioning_level);
        assert(actual[i].parent_path == expected[i].parent_path);
        assert(actual[i].phase_space_dimension == expected[i].phase_space_dimension);
        assert(actual[i].configured_node_count == expected[i].configured_node_count);
        assert(actual[i].active_node_count == expected[i].active_node_count);
    }
}

void assert_cqmom_quadrature_equal(const qmomkit::CqmomQuadrature& actual, const qmomkit::CqmomQuadrature& expected)
{
    assert(actual.hierarchy().node_counts() == expected.hierarchy().node_counts());
    assert(actual.hierarchy().conditioning_order() == expected.hierarchy().conditioning_order());
    assert(actual.adaptive_branch_counts_enabled() == expected.adaptive_branch_counts_enabled());
    assert(actual.has_rectangular_storage() == expected.has_rectangular_storage());
    assert(actual.terminal_paths() == expected.terminal_paths());
    assert_vector_near(actual.weights(), expected.weights());
    assert(actual.abscissas().size() == expected.abscissas().size());
    for (std::size_t dimension = 0; dimension < actual.abscissas().size(); ++dimension)
    {
        assert_vector_near(actual.abscissas()[dimension], expected.abscissas()[dimension]);
    }
}

void assert_cqmom_state_equal(const qmomkit::CqmomState& actual, const qmomkit::CqmomState& expected)
{
    assert(actual.hierarchy().node_counts() == expected.hierarchy().node_counts());
    assert(actual.hierarchy().conditioning_order() == expected.hierarchy().conditioning_order());
    assert(actual.qmom_options().adaptive == expected.qmom_options().adaptive);
    assert_near(actual.qmom_options().rmin, expected.qmom_options().rmin);
    assert_near(actual.qmom_options().eabs, expected.qmom_options().eabs);
    assert_near(actual.qmom_options().cutoff, expected.qmom_options().cutoff);
    assert_near(actual.conditional_options().singular_tolerance, expected.conditional_options().singular_tolerance);
    assert_near(actual.conditional_options().conditioning_warning_threshold, expected.conditional_options().conditioning_warning_threshold);
    assert(actual.adaptive_branch_counts() == expected.adaptive_branch_counts());
    assert(actual.diagnostic_level() == expected.diagnostic_level());
    assert(actual.vandermonde_solver() == expected.vandermonde_solver());
    assert(actual.coordinate_support().size() == expected.coordinate_support().size());
    for (std::size_t dimension = 0; dimension < actual.coordinate_support().size(); ++dimension)
    {
        assert(actual.coordinate_support()[dimension].has_value() == expected.coordinate_support()[dimension].has_value());
        if (actual.coordinate_support()[dimension].has_value())
        {
            assert_near(actual.coordinate_support()[dimension]->lower, expected.coordinate_support()[dimension]->lower);
            assert_near(actual.coordinate_support()[dimension]->upper, expected.coordinate_support()[dimension]->upper);
        }
    }
    assert(actual.moment_storage() == expected.moment_storage());
    assert(actual.moment_indices() == expected.moment_indices());
    assert_vector_near(actual.moment_values(), expected.moment_values());
    assert_active_counts_equal(actual.active_branch_counts(), expected.active_branch_counts());
}

void assert_cqmom_result_equal(const qmomkit::CqmomTransportResult& actual, const qmomkit::CqmomTransportResult& expected)
{
    assert_cqmom_state_equal(actual.state, expected.state);
    assert(actual.status == expected.status);
    assert(actual.quadrature.has_value() == expected.quadrature.has_value());
    if (actual.quadrature.has_value())
    {
        assert_cqmom_quadrature_equal(*actual.quadrature, *expected.quadrature);
    }
    assert_near(actual.final_time, expected.final_time);
    assert(actual.message == expected.message);
    assert_transport_diagnostics_equal(actual.diagnostics, expected.diagnostics);
}

qmomkit::QmomTransportOptions qmom_options()
{
    qmomkit::QmomTransportOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    options.scheme = qmomkit::TransportScheme::strang;
    return options;
}

qmomkit::CqmomTransportOptions cqmom_options()
{
    qmomkit::CqmomTransportOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    options.scheme = qmomkit::TransportScheme::strang;
    return options;
}

std::vector<qmomkit::QmomState> qmom_states()
{
    std::vector<qmomkit::QmomState> states;
    for (const std::pair<double, double>& data : std::vector<std::pair<double, double>>{{1.0, -0.5}, {1.5, 0.25}, {2.0, 1.0}})
    {
        states.push_back(qmomkit::QmomCase{qmomkit::MomentSet({data.first, data.first * data.second}), 1, {}}.initial_state());
    }
    return states;
}

qmomkit::CqmomState cqmom_state(const double weight, const double x, const double y, const std::string& phase_space_settings,
                                const std::string& reconstruction_settings)
{
    const std::string configuration = "[phase_space]\n"
                                      "node_counts = 1, 1\n" +
                                      phase_space_settings + "\n" + reconstruction_settings +
                                      "\n"
                                      "[moments]\n"
                                      "extents = 2, 2\n"
                                      "M[0,0] = " +
                                      std::to_string(weight) +
                                      "\n"
                                      "M[0,1] = " +
                                      std::to_string(weight * y) +
                                      "\n"
                                      "M[1,0] = " +
                                      std::to_string(weight * x) +
                                      "\n"
                                      "M[1,1] = " +
                                      std::to_string(weight * x * y) + "\n";
    const auto parsed = qmomkit::config::parse_cqmom(configuration, "transport-batch");
    assert(parsed.ok());
    return parsed.value().initial_state();
}

std::vector<qmomkit::CqmomState> cqmom_states()
{
    return {cqmom_state(1.0, -0.5, 0.25, "", ""),
            cqmom_state(1.5, 0.2, -0.4,
                        "adaptive_branch_counts = true\n"
                        "support_lower = none, -1.0\n"
                        "support_upper = none, 1.0\n",
                        "[qmom]\n"
                        "adaptive = true\n"
                        "\n"
                        "[cqmom]\n"
                        "diagnostic_level = basic\n"
                        "\n"
                        "[vandermonde]\n"
                        "solver = dense\n"),
            cqmom_state(2.0, 0.8, 1.0, "", "")};
}

qmomkit::PhaseSpaceVelocity ordinary_velocity()
{
    return [](double, const std::vector<double>& location)
    {
        std::vector<double> velocity(location.size(), 0.0);
        velocity[0] = 0.15;
        if (velocity.size() > 1)
        {
            velocity[1] = -0.1;
        }
        return velocity;
    };
}

qmomkit::MomentSource ordinary_qmom_source()
{
    return [](double, const std::vector<double>& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        for (std::size_t i = 0; i < moments.size(); ++i)
        {
            rate[i] = 0.1 * moments[i];
        }
        return rate;
    };
}

qmomkit::MomentBoundaryContribution ordinary_qmom_boundary()
{
    return [](double, const std::vector<double>& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        for (std::size_t i = 0; i < moments.size(); ++i)
        {
            rate[i] = 0.05 * moments[i];
        }
        return qmomkit::MomentBoundaryRate{rate};
    };
}

qmomkit::CqmomMomentSource ordinary_cqmom_source()
{
    return [](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t i = 0; i < context.size(); ++i)
        {
            rate[i] = 0.1 * context.value(i);
        }
        return qmomkit::CqmomMomentRate{rate};
    };
}

qmomkit::CqmomMomentBoundaryContribution ordinary_cqmom_boundary()
{
    return [](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t i = 0; i < context.size(); ++i)
        {
            rate[i] = 0.05 * context.value(i);
        }
        return qmomkit::CqmomMomentRate{rate};
    };
}

void test_qmom_ordinary_scalar_equivalence()
{
    const std::vector<qmomkit::QmomState> states = qmom_states();
    const qmomkit::QmomTransportOptions options = qmom_options();
    const qmomkit::PhaseSpaceVelocity velocity = ordinary_velocity();
    const qmomkit::MomentSource source = ordinary_qmom_source();
    const qmomkit::MomentBoundaryContribution boundary = ordinary_qmom_boundary();

    const auto continuous = qmomkit::advance_qmom_continuous_batch(states, 0.0, 0.2, velocity, options);
    const auto discontinuous = qmomkit::advance_qmom_discontinuous_batch(states, 1, 0.0, 0.2, source, options);
    const auto discontinuous_boundary = qmomkit::advance_qmom_discontinuous_batch(states, 1, 0.0, 0.2, source, boundary, options);
    const auto transport = qmomkit::advance_transport_qmom_batch(states, 1, 0.0, 0.2, velocity, source, options);
    const auto transport_boundary = qmomkit::advance_transport_qmom_batch(states, 1, 0.0, 0.2, velocity, source, boundary, options);

    assert(continuous.size() == states.size());
    assert(discontinuous.size() == states.size());
    assert(discontinuous_boundary.size() == states.size());
    assert(transport.size() == states.size());
    assert(transport_boundary.size() == states.size());
    for (std::size_t i = 0; i < states.size(); ++i)
    {
        assert_qmom_result_equal(continuous[i], qmomkit::advance_qmom_continuous(states[i], 0.0, 0.2, velocity, options));
        assert_qmom_result_equal(discontinuous[i], qmomkit::advance_qmom_discontinuous(states[i], 1, 0.0, 0.2, source, options));
        assert_qmom_result_equal(discontinuous_boundary[i], qmomkit::advance_qmom_discontinuous(states[i], 1, 0.0, 0.2, source, boundary, options));
        assert_qmom_result_equal(transport[i], qmomkit::advance_transport_qmom(states[i], 1, 0.0, 0.2, velocity, source, options));
        assert_qmom_result_equal(transport_boundary[i], qmomkit::advance_transport_qmom(states[i], 1, 0.0, 0.2, velocity, source, boundary, options));
        assert_qmom_state_equal(states[i], qmom_states()[i]);
    }
}

void test_cqmom_ordinary_scalar_equivalence()
{
    const std::vector<qmomkit::CqmomState> states = cqmom_states();
    const qmomkit::CqmomTransportOptions options = cqmom_options();
    const qmomkit::PhaseSpaceVelocity velocity = ordinary_velocity();
    const qmomkit::CqmomMomentSource source = ordinary_cqmom_source();
    const qmomkit::CqmomMomentBoundaryContribution boundary = ordinary_cqmom_boundary();

    const auto continuous = qmomkit::advance_cqmom_continuous_batch(states, 0.0, 0.2, velocity, options);
    const auto discontinuous = qmomkit::advance_cqmom_discontinuous_batch(states, 0.0, 0.2, source, options);
    const auto discontinuous_boundary = qmomkit::advance_cqmom_discontinuous_batch(states, 0.0, 0.2, source, boundary, options);
    const auto transport = qmomkit::advance_transport_cqmom_batch(states, 0.0, 0.2, velocity, source, options);
    const auto transport_boundary = qmomkit::advance_transport_cqmom_batch(states, 0.0, 0.2, velocity, source, boundary, options);

    assert(continuous.size() == states.size());
    assert(discontinuous.size() == states.size());
    assert(discontinuous_boundary.size() == states.size());
    assert(transport.size() == states.size());
    assert(transport_boundary.size() == states.size());
    for (std::size_t i = 0; i < states.size(); ++i)
    {
        assert_cqmom_result_equal(continuous[i], qmomkit::advance_cqmom_continuous(states[i], 0.0, 0.2, velocity, options));
        assert_cqmom_result_equal(discontinuous[i], qmomkit::advance_cqmom_discontinuous(states[i], 0.0, 0.2, source, options));
        assert_cqmom_result_equal(discontinuous_boundary[i], qmomkit::advance_cqmom_discontinuous(states[i], 0.0, 0.2, source, boundary, options));
        assert_cqmom_result_equal(transport[i], qmomkit::advance_transport_cqmom(states[i], 0.0, 0.2, velocity, source, options));
        assert_cqmom_result_equal(transport_boundary[i], qmomkit::advance_transport_cqmom(states[i], 0.0, 0.2, velocity, source, boundary, options));
        assert_cqmom_state_equal(states[i], cqmom_states()[i]);
    }
}

void test_qmom_batch_aware_scalar_equivalence()
{
    const std::vector<qmomkit::QmomState> states = qmom_states();
    const qmomkit::QmomTransportOptions options = qmom_options();
    const std::vector<double> external{0.02, 0.05, 0.09};
    const qmomkit::BatchPhaseSpaceVelocity velocity = [&external](const std::size_t batch_index, double, const std::vector<double>& location)
    { return std::vector<double>(location.size(), external[batch_index]); };
    const qmomkit::BatchMomentSource source = [&external](const std::size_t batch_index, double, const std::vector<double>& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        for (std::size_t i = 0; i < moments.size(); ++i)
        {
            rate[i] = external[batch_index] * moments[i];
        }
        return rate;
    };
    const qmomkit::BatchMomentBoundaryContribution boundary = [&external](const std::size_t batch_index, double, const std::vector<double>& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        for (std::size_t i = 0; i < moments.size(); ++i)
        {
            rate[i] = 0.5 * external[batch_index] * moments[i];
        }
        return qmomkit::MomentBoundaryRate{rate};
    };

    const auto continuous = qmomkit::advance_qmom_continuous_batch_with_index(states, 0.0, 0.2, velocity, options);
    const auto discontinuous = qmomkit::advance_qmom_discontinuous_batch_with_index(states, 1, 0.0, 0.2, source, options);
    const auto discontinuous_boundary = qmomkit::advance_qmom_discontinuous_batch_with_index(states, 1, 0.0, 0.2, source, boundary, options);
    const auto transport = qmomkit::advance_transport_qmom_batch_with_index(states, 1, 0.0, 0.2, velocity, source, options);
    const auto transport_boundary = qmomkit::advance_transport_qmom_batch_with_index(states, 1, 0.0, 0.2, velocity, source, boundary, options);

    assert(continuous.size() == states.size());
    assert(discontinuous.size() == states.size());
    assert(discontinuous_boundary.size() == states.size());
    assert(transport.size() == states.size());
    assert(transport_boundary.size() == states.size());
    for (std::size_t i = 0; i < states.size(); ++i)
    {
        const qmomkit::PhaseSpaceVelocity scalar_velocity = [&velocity, i](double time, const std::vector<double>& location)
        { return velocity(i, time, location); };
        const qmomkit::MomentSource scalar_source = [&source, i](double time, const std::vector<double>& moments)
        { return source(i, time, moments); };
        const qmomkit::MomentBoundaryContribution scalar_boundary = [&boundary, i](double time, const std::vector<double>& moments)
        { return boundary(i, time, moments); };
        assert_qmom_result_equal(continuous[i], qmomkit::advance_qmom_continuous(states[i], 0.0, 0.2, scalar_velocity, options));
        assert_qmom_result_equal(discontinuous[i], qmomkit::advance_qmom_discontinuous(states[i], 1, 0.0, 0.2, scalar_source, options));
        assert_qmom_result_equal(discontinuous_boundary[i],
                                 qmomkit::advance_qmom_discontinuous(states[i], 1, 0.0, 0.2, scalar_source, scalar_boundary, options));
        assert_qmom_result_equal(transport[i], qmomkit::advance_transport_qmom(states[i], 1, 0.0, 0.2, scalar_velocity, scalar_source, options));
        assert_qmom_result_equal(transport_boundary[i],
                                 qmomkit::advance_transport_qmom(states[i], 1, 0.0, 0.2, scalar_velocity, scalar_source, scalar_boundary, options));
    }
}

void test_cqmom_batch_aware_scalar_equivalence()
{
    const std::vector<qmomkit::CqmomState> states = cqmom_states();
    const qmomkit::CqmomTransportOptions options = cqmom_options();
    const std::vector<double> external{0.02, 0.05, 0.09};
    const qmomkit::BatchPhaseSpaceVelocity velocity = [&external](const std::size_t batch_index, double, const std::vector<double>& location)
    {
        std::vector<double> rate(location.size(), 0.0);
        rate[0] = external[batch_index];
        if (rate.size() > 1)
        {
            rate[1] = -0.5 * external[batch_index];
        }
        return rate;
    };
    const qmomkit::BatchCqmomMomentSource source = [&external](const std::size_t batch_index, double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t i = 0; i < context.size(); ++i)
        {
            rate[i] = external[batch_index] * context.value(i);
        }
        return qmomkit::CqmomMomentRate{rate};
    };
    const qmomkit::BatchCqmomMomentBoundaryContribution boundary =
        [&external](const std::size_t batch_index, double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t i = 0; i < context.size(); ++i)
        {
            rate[i] = 0.5 * external[batch_index] * context.value(i);
        }
        return qmomkit::CqmomMomentRate{rate};
    };

    const auto continuous = qmomkit::advance_cqmom_continuous_batch_with_index(states, 0.0, 0.2, velocity, options);
    const auto discontinuous = qmomkit::advance_cqmom_discontinuous_batch_with_index(states, 0.0, 0.2, source, options);
    const auto discontinuous_boundary = qmomkit::advance_cqmom_discontinuous_batch_with_index(states, 0.0, 0.2, source, boundary, options);
    const auto transport = qmomkit::advance_transport_cqmom_batch_with_index(states, 0.0, 0.2, velocity, source, options);
    const auto transport_boundary = qmomkit::advance_transport_cqmom_batch_with_index(states, 0.0, 0.2, velocity, source, boundary, options);

    assert(continuous.size() == states.size());
    assert(discontinuous.size() == states.size());
    assert(discontinuous_boundary.size() == states.size());
    assert(transport.size() == states.size());
    assert(transport_boundary.size() == states.size());
    for (std::size_t i = 0; i < states.size(); ++i)
    {
        const qmomkit::PhaseSpaceVelocity scalar_velocity = [&velocity, i](double time, const std::vector<double>& location)
        { return velocity(i, time, location); };
        const qmomkit::CqmomMomentSource scalar_source = [&source, i](double time, const qmomkit::CqmomSourceContext& context)
        { return source(i, time, context); };
        const qmomkit::CqmomMomentBoundaryContribution scalar_boundary = [&boundary, i](double time, const qmomkit::CqmomSourceContext& context)
        { return boundary(i, time, context); };
        assert_cqmom_result_equal(continuous[i], qmomkit::advance_cqmom_continuous(states[i], 0.0, 0.2, scalar_velocity, options));
        assert_cqmom_result_equal(discontinuous[i], qmomkit::advance_cqmom_discontinuous(states[i], 0.0, 0.2, scalar_source, options));
        assert_cqmom_result_equal(discontinuous_boundary[i],
                                  qmomkit::advance_cqmom_discontinuous(states[i], 0.0, 0.2, scalar_source, scalar_boundary, options));
        assert_cqmom_result_equal(transport[i], qmomkit::advance_transport_cqmom(states[i], 0.0, 0.2, scalar_velocity, scalar_source, options));
        assert_cqmom_result_equal(transport_boundary[i],
                                  qmomkit::advance_transport_cqmom(states[i], 0.0, 0.2, scalar_velocity, scalar_source, scalar_boundary, options));
    }
}

void test_empty_batch_skips_scalar_validation_and_callbacks()
{
    const std::vector<qmomkit::QmomState> states = qmom_states();
    const std::vector<qmomkit::QmomState> empty_states;
    const qmomkit::QmomTransportOptions options = qmom_options();
    const double invalid_time = std::numeric_limits<double>::quiet_NaN();
    std::size_t ordinary_calls = 0;
    std::size_t indexed_calls = 0;
    const qmomkit::PhaseSpaceVelocity ordinary_velocity = [&ordinary_calls](double, const std::vector<double>& location)
    {
        ++ordinary_calls;
        return std::vector<double>(location.size(), 0.0);
    };
    const qmomkit::BatchPhaseSpaceVelocity indexed_velocity = [&indexed_calls](std::size_t, double, const std::vector<double>& location)
    {
        ++indexed_calls;
        return std::vector<double>(location.size(), 0.0);
    };

    // A non-finite time is rejected before velocity evaluation by scalar continuous transport; empty batches must skip that scalar validation.
    const qmomkit::QmomTransportResult scalar = qmomkit::advance_qmom_continuous(states.front(), invalid_time, 0.2, ordinary_velocity, options);
    assert(!scalar.ok());
    assert(scalar.status == qmomkit::TransportStatus::invalid_input);
    assert(ordinary_calls == 0);

    const auto ordinary = qmomkit::advance_qmom_continuous_batch(empty_states, invalid_time, 0.2, ordinary_velocity, options);
    const auto indexed = qmomkit::advance_qmom_continuous_batch_with_index(empty_states, invalid_time, 0.2, indexed_velocity, options);
    assert(ordinary.empty());
    assert(indexed.empty());
    assert(ordinary_calls == 0);
    assert(indexed_calls == 0);
}

void test_normal_failures_are_independent_and_ordered()
{
    const std::vector<qmomkit::QmomState> states = qmom_states();
    const std::vector<qmomkit::QmomState> initial_states = states;
    const qmomkit::QmomTransportOptions options = qmom_options();
    const qmomkit::BatchPhaseSpaceVelocity velocity = [](const std::size_t batch_index, double, const std::vector<double>& location)
    {
        if (batch_index == 1)
        {
            return std::vector<double>(location.size() + 1, 0.0);
        }
        return std::vector<double>(location.size(), 0.0);
    };
    const auto batch = qmomkit::advance_qmom_continuous_batch_with_index(states, 0.0, 0.2, velocity, options);
    assert(batch.size() == states.size());
    for (std::size_t i = 0; i < states.size(); ++i)
    {
        const qmomkit::PhaseSpaceVelocity scalar_velocity = [&velocity, i](double time, const std::vector<double>& location)
        { return velocity(i, time, location); };
        assert_qmom_result_equal(batch[i], qmomkit::advance_qmom_continuous(states[i], 0.0, 0.2, scalar_velocity, options));
    }
    assert(batch[0].ok());
    assert(!batch[1].ok());
    assert(batch[1].status == qmomkit::TransportStatus::continuous_failed);
    assert(batch[2].ok());
    for (std::size_t i = 0; i < states.size(); ++i)
    {
        assert_qmom_state_equal(states[i], initial_states[i]);
    }
}

void test_escaping_exception_stops_at_the_throwing_entry()
{
    const std::vector<qmomkit::QmomState> states = qmom_states();
    const qmomkit::QmomTransportOptions options = qmom_options();
    std::vector<std::size_t> invoked;
    const qmomkit::BatchPhaseSpaceVelocity velocity = [&invoked](const std::size_t batch_index, double, const std::vector<double>& location)
    {
        invoked.push_back(batch_index);
        if (batch_index == 1)
        {
            throw 37;
        }
        return std::vector<double>(location.size(), 0.0);
    };

    bool caught = false;
    try
    {
        (void)qmomkit::advance_qmom_continuous_batch_with_index(states, 0.0, 0.2, velocity, options);
    }
    catch (const int value)
    {
        assert(value == 37);
        caught = true;
    }
    assert(caught);
    assert(!invoked.empty());
    assert(invoked.front() == 0);
    assert(invoked.back() == 1);
    for (const std::size_t batch_index : invoked)
    {
        assert(batch_index != 2);
    }
}

void test_batch_index_remains_stable_for_repeated_callbacks()
{
    const std::vector<qmomkit::QmomState> states = qmom_states();
    const qmomkit::QmomTransportOptions options = qmom_options();
    std::vector<std::vector<std::size_t>> observed(states.size());
    const qmomkit::BatchPhaseSpaceVelocity velocity = [&observed](const std::size_t batch_index, double, const std::vector<double>& location)
    {
        observed[batch_index].push_back(batch_index);
        return std::vector<double>(location.size(), 0.0);
    };
    const auto results = qmomkit::advance_qmom_continuous_batch_with_index(states, 0.0, 0.2, velocity, options);
    assert(results.size() == states.size());
    for (std::size_t i = 0; i < observed.size(); ++i)
    {
        assert(observed[i].size() >= 3);
        for (const std::size_t batch_index : observed[i])
        {
            assert(batch_index == i);
        }
    }
}

void test_cqmom_normal_failures_are_independent_and_ordered()
{
    const std::vector<qmomkit::CqmomState> states = cqmom_states();
    const std::vector<qmomkit::CqmomState> initial_states = states;
    const qmomkit::CqmomTransportOptions options = cqmom_options();
    const qmomkit::BatchPhaseSpaceVelocity velocity = [](const std::size_t batch_index, double, const std::vector<double>& location)
    {
        if (batch_index == 1)
        {
            return std::vector<double>(location.size() + 1, 0.0);
        }
        return std::vector<double>(location.size(), 0.0);
    };

    const auto batch = qmomkit::advance_cqmom_continuous_batch_with_index(states, 0.0, 0.2, velocity, options);
    assert(batch.size() == states.size());
    for (std::size_t i = 0; i < states.size(); ++i)
    {
        const qmomkit::PhaseSpaceVelocity scalar_velocity = [&velocity, i](double time, const std::vector<double>& location)
        { return velocity(i, time, location); };
        assert_cqmom_result_equal(batch[i], qmomkit::advance_cqmom_continuous(states[i], 0.0, 0.2, scalar_velocity, options));
        assert_cqmom_state_equal(states[i], initial_states[i]);
    }
    assert(batch[0].ok());
    assert(!batch[1].ok());
    assert(batch[1].status == qmomkit::TransportStatus::continuous_failed);
    assert(batch[2].ok());
}

struct CqmomCallbackException
{
    int value = 0;
};

void test_cqmom_escaping_exception_stops_at_the_throwing_entry()
{
    const std::vector<qmomkit::CqmomState> states = cqmom_states();
    const qmomkit::CqmomTransportOptions options = cqmom_options();
    std::vector<std::size_t> invoked;
    const qmomkit::BatchPhaseSpaceVelocity velocity = [&invoked](const std::size_t batch_index, double, const std::vector<double>& location)
    {
        invoked.push_back(batch_index);
        if (batch_index == 1)
        {
            throw CqmomCallbackException{91};
        }
        return std::vector<double>(location.size(), 0.0);
    };

    bool caught = false;
    try
    {
        (void)qmomkit::advance_cqmom_continuous_batch_with_index(states, 0.0, 0.2, velocity, options);
    }
    catch (const CqmomCallbackException& error)
    {
        assert(error.value == 91);
        caught = true;
    }
    assert(caught);
    assert(!invoked.empty());
    assert(invoked.front() == 0);
    assert(invoked.back() == 1);
    for (const std::size_t batch_index : invoked)
    {
        assert(batch_index != 2);
    }
}

void test_cqmom_batch_index_remains_stable_for_repeated_callbacks()
{
    const std::vector<qmomkit::CqmomState> states = cqmom_states();
    const qmomkit::CqmomTransportOptions options = cqmom_options();
    std::vector<std::vector<std::size_t>> observed(states.size());
    const qmomkit::BatchPhaseSpaceVelocity velocity = [&observed](const std::size_t batch_index, double, const std::vector<double>& location)
    {
        observed[batch_index].push_back(batch_index);
        return std::vector<double>(location.size(), 0.0);
    };

    const auto results = qmomkit::advance_cqmom_continuous_batch_with_index(states, 0.0, 0.2, velocity, options);
    assert(results.size() == states.size());
    for (std::size_t i = 0; i < observed.size(); ++i)
    {
        assert(observed[i].size() >= 3);
        for (const std::size_t batch_index : observed[i])
        {
            assert(batch_index == i);
        }
    }
}

void test_null_callbacks_preserve_scalar_behavior()
{
    const std::vector<qmomkit::QmomState> states = qmom_states();
    const qmomkit::QmomTransportOptions options = qmom_options();
    const qmomkit::PhaseSpaceVelocity velocity = ordinary_velocity();
    const qmomkit::MomentSource source = ordinary_qmom_source();
    const qmomkit::PhaseSpaceVelocity empty_scalar_velocity;
    const qmomkit::MomentSource empty_scalar_source;
    const qmomkit::MomentBoundaryContribution empty_scalar_boundary;
    const qmomkit::BatchPhaseSpaceVelocity indexed_velocity = [&velocity](std::size_t, const double time, const std::vector<double>& location)
    { return velocity(time, location); };
    const qmomkit::BatchMomentSource indexed_source = [&source](std::size_t, const double time, const std::vector<double>& moments)
    { return source(time, moments); };

    const auto ordinary_continuous = qmomkit::advance_qmom_continuous_batch(states, 0.0, 0.2, nullptr, options);
    const auto indexed_continuous = qmomkit::advance_qmom_continuous_batch_with_index(states, 0.0, 0.2, nullptr, options);
    const auto ordinary_discontinuous = qmomkit::advance_qmom_discontinuous_batch(states, 1, 0.0, 0.2, nullptr, options);
    const auto indexed_discontinuous = qmomkit::advance_qmom_discontinuous_batch_with_index(states, 1, 0.0, 0.2, nullptr, options);
    const auto ordinary_discontinuous_boundary = qmomkit::advance_qmom_discontinuous_batch(states, 1, 0.0, 0.2, source, nullptr, options);
    const auto indexed_discontinuous_boundary =
        qmomkit::advance_qmom_discontinuous_batch_with_index(states, 1, 0.0, 0.2, indexed_source, nullptr, options);
    const auto ordinary_transport_velocity = qmomkit::advance_transport_qmom_batch(states, 1, 0.0, 0.2, nullptr, source, options);
    const auto indexed_transport_velocity = qmomkit::advance_transport_qmom_batch_with_index(states, 1, 0.0, 0.2, nullptr, indexed_source, options);
    const auto ordinary_transport_source = qmomkit::advance_transport_qmom_batch(states, 1, 0.0, 0.2, velocity, nullptr, options);
    const auto indexed_transport_source = qmomkit::advance_transport_qmom_batch_with_index(states, 1, 0.0, 0.2, indexed_velocity, nullptr, options);
    const auto ordinary_transport_boundary = qmomkit::advance_transport_qmom_batch(states, 1, 0.0, 0.2, velocity, source, nullptr, options);
    const auto indexed_transport_boundary =
        qmomkit::advance_transport_qmom_batch_with_index(states, 1, 0.0, 0.2, indexed_velocity, indexed_source, nullptr, options);

    assert(ordinary_continuous.size() == states.size());
    assert(indexed_continuous.size() == states.size());
    assert(ordinary_discontinuous.size() == states.size());
    assert(indexed_discontinuous.size() == states.size());
    assert(ordinary_discontinuous_boundary.size() == states.size());
    assert(indexed_discontinuous_boundary.size() == states.size());
    assert(ordinary_transport_velocity.size() == states.size());
    assert(indexed_transport_velocity.size() == states.size());
    assert(ordinary_transport_source.size() == states.size());
    assert(indexed_transport_source.size() == states.size());
    assert(ordinary_transport_boundary.size() == states.size());
    assert(indexed_transport_boundary.size() == states.size());
    for (std::size_t i = 0; i < states.size(); ++i)
    {
        const qmomkit::QmomTransportResult scalar_continuous = qmomkit::advance_qmom_continuous(states[i], 0.0, 0.2, empty_scalar_velocity, options);
        const qmomkit::QmomTransportResult scalar_discontinuous =
            qmomkit::advance_qmom_discontinuous(states[i], 1, 0.0, 0.2, empty_scalar_source, options);
        const qmomkit::QmomTransportResult scalar_discontinuous_boundary =
            qmomkit::advance_qmom_discontinuous(states[i], 1, 0.0, 0.2, source, empty_scalar_boundary, options);
        const qmomkit::QmomTransportResult scalar_transport_velocity =
            qmomkit::advance_transport_qmom(states[i], 1, 0.0, 0.2, empty_scalar_velocity, source, options);
        const qmomkit::QmomTransportResult scalar_transport_source =
            qmomkit::advance_transport_qmom(states[i], 1, 0.0, 0.2, velocity, empty_scalar_source, options);
        const qmomkit::QmomTransportResult scalar_transport_boundary =
            qmomkit::advance_transport_qmom(states[i], 1, 0.0, 0.2, velocity, source, empty_scalar_boundary, options);

        assert_qmom_result_equal(ordinary_continuous[i], scalar_continuous);
        assert_qmom_result_equal(indexed_continuous[i], scalar_continuous);
        assert_qmom_result_equal(ordinary_discontinuous[i], scalar_discontinuous);
        assert_qmom_result_equal(indexed_discontinuous[i], scalar_discontinuous);
        assert_qmom_result_equal(ordinary_discontinuous_boundary[i], scalar_discontinuous_boundary);
        assert_qmom_result_equal(indexed_discontinuous_boundary[i], scalar_discontinuous_boundary);
        assert_qmom_result_equal(ordinary_transport_velocity[i], scalar_transport_velocity);
        assert_qmom_result_equal(indexed_transport_velocity[i], scalar_transport_velocity);
        assert_qmom_result_equal(ordinary_transport_source[i], scalar_transport_source);
        assert_qmom_result_equal(indexed_transport_source[i], scalar_transport_source);
        assert_qmom_result_equal(ordinary_transport_boundary[i], scalar_transport_boundary);
        assert_qmom_result_equal(indexed_transport_boundary[i], scalar_transport_boundary);
    }
}

} // namespace

int main()
{
    test_qmom_ordinary_scalar_equivalence();
    test_cqmom_ordinary_scalar_equivalence();
    test_qmom_batch_aware_scalar_equivalence();
    test_cqmom_batch_aware_scalar_equivalence();
    test_empty_batch_skips_scalar_validation_and_callbacks();
    test_normal_failures_are_independent_and_ordered();
    test_escaping_exception_stops_at_the_throwing_entry();
    test_batch_index_remains_stable_for_repeated_callbacks();
    test_cqmom_normal_failures_are_independent_and_ordered();
    test_cqmom_escaping_exception_stops_at_the_throwing_entry();
    test_cqmom_batch_index_remains_stable_for_repeated_callbacks();
    test_null_callbacks_preserve_scalar_behavior();
    return 0;
}
