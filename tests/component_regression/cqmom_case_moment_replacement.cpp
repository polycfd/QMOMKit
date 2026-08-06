#include "../verification/verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"

#include <cstddef>
#include <vector>

namespace
{

const char* base_configuration = R"(
[phase_space]
node_counts = 1, 1
conditioning_order = 1, 0
adaptive_branch_counts = true
support_lower = none, 0.0
support_upper = none, 1.0

[qmom]
adaptive = true

[cqmom]
moment_storage = lean
diagnostic_level = detailed

[vandermonde]
solver = dense

[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 0.25
M[1,0] = 0.0
M[1,1] = 0.0
)";

const char* explicit_replacement_configuration = R"(
[phase_space]
node_counts = 1, 1
conditioning_order = 1, 0
adaptive_branch_counts = true
support_lower = none, 0.0
support_upper = none, 1.0

[qmom]
adaptive = true

[cqmom]
moment_storage = lean
diagnostic_level = detailed

[vandermonde]
solver = dense

[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 0.5
M[1,0] = 1.0
M[1,1] = 0.5
)";

qmomkit::RawMomentTensor replacement_tensor()
{
    return {{2, 2}, {1.0, 0.5, 1.0, 0.5}};
}

void assert_cases_match(const qmomkit::CqmomCase& actual, const qmomkit::CqmomCase& expected)
{
    QMOMKIT_VERIFY_CHECK(actual.hierarchy().node_counts() == expected.hierarchy().node_counts());
    QMOMKIT_VERIFY_CHECK(actual.hierarchy().conditioning_order() == expected.hierarchy().conditioning_order());
    QMOMKIT_VERIFY_CHECK(actual.qmom_options().adaptive == expected.qmom_options().adaptive);
    QMOMKIT_VERIFY_NEAR(actual.qmom_options().rmin, expected.qmom_options().rmin, 0.0);
    QMOMKIT_VERIFY_NEAR(actual.qmom_options().eabs, expected.qmom_options().eabs, 0.0);
    QMOMKIT_VERIFY_NEAR(actual.qmom_options().cutoff, expected.qmom_options().cutoff, 0.0);
    QMOMKIT_VERIFY_NEAR(actual.conditional_options().singular_tolerance, expected.conditional_options().singular_tolerance, 0.0);
    QMOMKIT_VERIFY_NEAR(actual.conditional_options().conditioning_warning_threshold, expected.conditional_options().conditioning_warning_threshold,
                        0.0);
    QMOMKIT_VERIFY_CHECK(actual.adaptive_branch_counts() == expected.adaptive_branch_counts());
    QMOMKIT_VERIFY_CHECK(actual.diagnostic_level() == expected.diagnostic_level());
    QMOMKIT_VERIFY_CHECK(actual.vandermonde_solver() == expected.vandermonde_solver());
    QMOMKIT_VERIFY_CHECK(actual.moment_storage() == expected.moment_storage());
    QMOMKIT_VERIFY_CHECK(actual.coordinate_support().size() == expected.coordinate_support().size());
    for (std::size_t dimension = 0; dimension < actual.coordinate_support().size(); ++dimension)
    {
        QMOMKIT_VERIFY_CHECK(actual.coordinate_support()[dimension].has_value() == expected.coordinate_support()[dimension].has_value());
        if (actual.coordinate_support()[dimension].has_value())
        {
            QMOMKIT_VERIFY_NEAR(actual.coordinate_support()[dimension]->lower, expected.coordinate_support()[dimension]->lower, 0.0);
            QMOMKIT_VERIFY_NEAR(actual.coordinate_support()[dimension]->upper, expected.coordinate_support()[dimension]->upper, 0.0);
        }
    }
}

void assert_quadratures_match(const qmomkit::CqmomQuadrature& actual, const qmomkit::CqmomQuadrature& expected)
{
    QMOMKIT_VERIFY_CHECK(actual.weights().size() == expected.weights().size());
    QMOMKIT_VERIFY_CHECK(actual.abscissas().size() == expected.abscissas().size());
    for (std::size_t node = 0; node < actual.weights().size(); ++node)
    {
        QMOMKIT_VERIFY_NEAR(actual.weights()[node], expected.weights()[node], 1.0e-13);
    }
    for (std::size_t dimension = 0; dimension < actual.abscissas().size(); ++dimension)
    {
        QMOMKIT_VERIFY_CHECK(actual.abscissas()[dimension].size() == expected.abscissas()[dimension].size());
        for (std::size_t node = 0; node < actual.abscissas()[dimension].size(); ++node)
        {
            QMOMKIT_VERIFY_NEAR(actual.abscissas()[dimension][node], expected.abscissas()[dimension][node], 1.0e-13);
        }
    }
}

void assert_diagnostics_match(const std::vector<qmomkit::CqmomDiagnostics>& actual, const std::vector<qmomkit::CqmomDiagnostics>& expected)
{
    QMOMKIT_VERIFY_CHECK(actual.size() == expected.size());
    for (std::size_t index = 0; index < actual.size(); ++index)
    {
        const qmomkit::CqmomDiagnostics& lhs = actual[index];
        const qmomkit::CqmomDiagnostics& rhs = expected[index];
        QMOMKIT_VERIFY_CHECK(lhs.conditioning_level == rhs.conditioning_level);
        QMOMKIT_VERIFY_CHECK(lhs.parent_path.size() == rhs.parent_path.size());
        QMOMKIT_VERIFY_CHECK(lhs.phase_space_dimension == rhs.phase_space_dimension);
        QMOMKIT_VERIFY_CHECK(lhs.requested_node_count == rhs.requested_node_count);
        QMOMKIT_VERIFY_CHECK(lhs.actual_node_count == rhs.actual_node_count);
        QMOMKIT_VERIFY_CHECK(lhs.adaptive_branch_counts_enabled == rhs.adaptive_branch_counts_enabled);
        QMOMKIT_VERIFY_CHECK(lhs.adaptive_node_count_reduced == rhs.adaptive_node_count_reduced);
        QMOMKIT_VERIFY_CHECK(lhs.qmom_status == rhs.qmom_status);
        QMOMKIT_VERIFY_CHECK(lhs.conditional_moment_status == rhs.conditional_moment_status);
        QMOMKIT_VERIFY_CHECK(lhs.message == rhs.message);
        QMOMKIT_VERIFY_CHECK(lhs.bounded_support_checked == rhs.bounded_support_checked);
        QMOMKIT_VERIFY_CHECK(lhs.bounded_support_interval_present == rhs.bounded_support_interval_present);
        QMOMKIT_VERIFY_CHECK(lhs.bounded_realizability_status == rhs.bounded_realizability_status);
    }
}

void assert_states_match(const qmomkit::CqmomState& actual, const qmomkit::CqmomState& expected)
{
    QMOMKIT_VERIFY_CHECK(actual.moment_storage() == expected.moment_storage());
    QMOMKIT_VERIFY_CHECK(actual.moment_indices() == expected.moment_indices());
    QMOMKIT_VERIFY_CHECK(actual.moment_values().size() == expected.moment_values().size());
    for (std::size_t index = 0; index < actual.moment_values().size(); ++index)
    {
        QMOMKIT_VERIFY_NEAR(actual.moment_values()[index], expected.moment_values()[index], 1.0e-13);
    }
    QMOMKIT_VERIFY_CHECK(actual.active_branch_counts().size() == expected.active_branch_counts().size());
    for (std::size_t index = 0; index < actual.active_branch_counts().size(); ++index)
    {
        const qmomkit::CqmomActiveBranchCount& lhs = actual.active_branch_counts()[index];
        const qmomkit::CqmomActiveBranchCount& rhs = expected.active_branch_counts()[index];
        QMOMKIT_VERIFY_CHECK(lhs.conditioning_level == rhs.conditioning_level);
        QMOMKIT_VERIFY_CHECK(lhs.parent_path.size() == rhs.parent_path.size());
        QMOMKIT_VERIFY_CHECK(lhs.phase_space_dimension == rhs.phase_space_dimension);
        QMOMKIT_VERIFY_CHECK(lhs.configured_node_count == rhs.configured_node_count);
        QMOMKIT_VERIFY_CHECK(lhs.active_node_count == rhs.active_node_count);
    }
}

void test_parsed_and_programmatic_replacement_match()
{
    const auto base = qmomkit::config::parse_cqmom(base_configuration, "replacement-base");
    const auto parsed = qmomkit::config::parse_cqmom(explicit_replacement_configuration, "replacement-parsed");
    QMOMKIT_VERIFY_CHECK(base.ok());
    QMOMKIT_VERIFY_CHECK(parsed.ok());

    const auto replaced = qmomkit::config::with_initial_moments(base.value(), replacement_tensor());
    QMOMKIT_VERIFY_CHECK(replaced.ok());
    assert_cases_match(replaced.value(), parsed.value());

    const qmomkit::CqmomResult replaced_result = qmomkit::compute_cqmom(replaced.value());
    const qmomkit::CqmomResult parsed_result = qmomkit::compute_cqmom(parsed.value());
    QMOMKIT_VERIFY_CHECK(replaced_result.ok());
    QMOMKIT_VERIFY_CHECK(parsed_result.ok());
    assert_quadratures_match(replaced_result.quadrature, parsed_result.quadrature);
    assert_diagnostics_match(replaced_result.diagnostics, parsed_result.diagnostics);

    const qmomkit::CqmomState replaced_state = replaced.value().initial_state();
    const qmomkit::CqmomState parsed_state = parsed.value().initial_state();
    assert_states_match(replaced_state, parsed_state);

    const qmomkit::PhaseSpaceVelocity zero_velocity = [](double, const std::vector<double>& coordinates)
    { return std::vector<double>(coordinates.size(), 0.0); };
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::CqmomPhaseSpaceEvolutionResult replaced_phase_space_evolution =
        qmomkit::advance_cqmom_phase_space(replaced_state, 0.0, 0.1, zero_velocity, zero_source, options);
    const qmomkit::CqmomPhaseSpaceEvolutionResult parsed_phase_space_evolution =
        qmomkit::advance_cqmom_phase_space(parsed_state, 0.0, 0.1, zero_velocity, zero_source, options);
    QMOMKIT_VERIFY_CHECK(replaced_phase_space_evolution.ok());
    QMOMKIT_VERIFY_CHECK(parsed_phase_space_evolution.ok());
    QMOMKIT_VERIFY_CHECK(replaced_phase_space_evolution.quadrature.has_value());
    QMOMKIT_VERIFY_CHECK(parsed_phase_space_evolution.quadrature.has_value());
    assert_quadratures_match(*replaced_phase_space_evolution.quadrature, *parsed_phase_space_evolution.quadrature);
    assert_states_match(replaced_phase_space_evolution.state, parsed_phase_space_evolution.state);
}

} // namespace

int main()
{
    test_parsed_and_programmatic_replacement_match();
    return 0;
}
