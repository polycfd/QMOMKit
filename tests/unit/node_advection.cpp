#include "qmomkit/qmomkit.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

void assert_contains(const std::string& text, const std::string& needle)
{
    if (text.find(needle) == std::string::npos)
    {
        std::cerr << "Missing expected text:\n" << needle << "\nActual message:\n" << text << '\n';
        assert(false);
    }
}

std::vector<double> raw_moments_1d(const std::vector<double>& abscissas, const std::vector<double>& weights, const std::size_t moment_count)
{
    std::vector<double> moments(moment_count, 0.0);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        for (std::size_t node = 0; node < abscissas.size(); ++node)
        {
            moments[order] += weights[node] * std::pow(abscissas[node], static_cast<int>(order));
        }
    }
    return moments;
}

std::vector<double> raw_moments_from_terminals(const std::vector<std::size_t>& extents, const std::vector<double>& weights,
                                               const std::vector<std::vector<double>>& abscissas)
{
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            double contribution = weights[node];
            for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
            {
                contribution *= std::pow(abscissas[dimension][node], static_cast<int>(moment_index[dimension]));
            }
            moments[flat_moment] += contribution;
        }
    }
    return moments;
}

qmomkit::CqmomResult make_fixed_two_dimensional_cqmom()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    return qmomkit::compute_cqmom({hierarchy, {extents, raw_moments_from_terminals(extents, weights, abscissas)}, {}, {}});
}

qmomkit::CqmomResult make_adaptive_two_dimensional_cqmom()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.1, 0.3, 0.6};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 1.0}, {-2.0, 2.0, 3.0}};
    return qmomkit::compute_cqmom(
        {hierarchy, {extents, raw_moments_from_terminals(extents, weights, abscissas)}, {}, {true, 1.0e-10, 1.0e-10, 1.0e-12}, true});
}

void test_adaptive_quadrature_structure_is_preserved()
{
    const qmomkit::CqmomResult cqmom_result = make_adaptive_two_dimensional_cqmom();
    assert(cqmom_result.ok());
    assert(!cqmom_result.quadrature.has_rectangular_storage());

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::adaptive_rk23;
    options.rk23_options.relative_tolerance = 1.0e-8;
    options.rk23_options.absolute_tolerance = 1.0e-11;

    constexpr double rate_x = 0.2;
    constexpr double rate_y = -0.6;
    constexpr double dt = 0.4;
    const qmomkit::CqmomNodeAdvectionResult result = qmomkit::advance_nodes(
        cqmom_result.quadrature, 0.0, dt,
        [](double, const std::vector<double>& xi)
        {
            assert(xi.size() == 2);
            return std::vector<double>{rate_x * xi[0], rate_y * xi[1]};
        },
        options);

    assert(result.ok());
    assert(result.quadrature.adaptive_branch_counts_enabled());
    assert(!result.quadrature.has_rectangular_storage());
    assert(result.quadrature.terminal_paths() == cqmom_result.quadrature.terminal_paths());
    assert(result.quadrature.weights() == cqmom_result.quadrature.weights());
    assert(result.quadrature.total_nodes() == 3);

    for (std::size_t node = 0; node < result.quadrature.total_nodes(); ++node)
    {
        assert_near(result.quadrature.abscissas()[0][node], cqmom_result.quadrature.abscissas()[0][node] * std::exp(rate_x * dt), 2.0e-8);
        assert_near(result.quadrature.abscissas()[1][node], cqmom_result.quadrature.abscissas()[1][node] * std::exp(rate_y * dt), 2.0e-8);
    }
}

void test_default_scheme_is_ssp_rk3()
{
    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {1.0, 2.0};
    quadrature.weights = {0.25, 0.75};

    const qmomkit::QmomNodeAdvectionResult result =
        qmomkit::advance_nodes(quadrature, 2.0, 0.5, [](double, const std::vector<double>&) { return std::vector<double>{3.0}; });

    assert(result.ok());
    assert(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);
    assert(result.diagnostics.step_count == 1);
    assert(result.diagnostics.rejected_step_count == 0);
    assert(result.diagnostics.rhs_evaluations == 3);
    assert(!result.diagnostics.adaptive_rk23_status.has_value());
    assert_near(result.diagnostics.final_time, 2.5, 1.0e-15);
    assert_near(result.quadrature.abscissas[0], 2.5, 1.0e-15);
    assert_near(result.quadrature.abscissas[1], 3.5, 1.0e-15);
    assert(result.quadrature.weights == quadrature.weights);
}

void test_invalid_input_status_covers_quadrature_and_callback_contracts()
{
    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {1.0};
    quadrature.weights = {1.0};

    const qmomkit::QmomNodeAdvectionResult invalid_qmom =
        qmomkit::advance_nodes(qmomkit::Quadrature1D{{1.0}, {1.0}, qmomkit::QuadratureStatus::not_realizable, "manufactured unsuccessful quadrature"},
                               0.0, 0.1, [](double, const std::vector<double>&) { return std::vector<double>{0.0}; });
    assert(!invalid_qmom.ok());
    assert(invalid_qmom.status == qmomkit::NodeAdvectionStatus::invalid_input);
    assert_contains(invalid_qmom.message, "Input QMOM quadrature is not successful.");

    const qmomkit::CqmomResult valid_cqmom = make_fixed_two_dimensional_cqmom();
    const qmomkit::CqmomNodeAdvectionResult invalid_cqmom_input =
        qmomkit::advance_nodes(valid_cqmom.quadrature, std::numeric_limits<double>::quiet_NaN(), 0.1,
                               [](double, const std::vector<double>&) { return std::vector<double>{0.0, 0.0}; });
    assert(!invalid_cqmom_input.ok());
    assert(invalid_cqmom_input.status == qmomkit::NodeAdvectionStatus::invalid_input);
    assert_contains(invalid_cqmom_input.message, "time input contains non-finite values");

    const qmomkit::PhaseSpaceVelocity missing_velocity;
    const qmomkit::QmomNodeAdvectionResult missing_callback = qmomkit::advance_nodes(quadrature, 0.0, 0.1, missing_velocity);
    assert(!missing_callback.ok());
    assert(missing_callback.status == qmomkit::NodeAdvectionStatus::invalid_input);
    assert_contains(missing_callback.message, "velocity callback must be provided");
}

void test_wrong_size_velocity_diagnostic_reports_expected_and_returned_size()
{
    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {1.0};
    quadrature.weights = {1.0};

    const qmomkit::QmomNodeAdvectionResult result =
        qmomkit::advance_nodes(quadrature, 0.0, 0.25, [](double, const std::vector<double>&) { return std::vector<double>{1.0, 2.0}; });

    assert(!result.ok());
    assert(result.status == qmomkit::NodeAdvectionStatus::velocity_failed);
    assert(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);
    assert(result.diagnostics.message == result.message);
    assert(result.quadrature.abscissas == quadrature.abscissas);
    assert(result.quadrature.weights == quadrature.weights);
    assert_contains(result.message, "Phase-space velocity callback returned the wrong dimension:");
    assert_contains(result.message, "expected size: 1");
    assert_contains(result.message, "returned size: 2");
}

void test_non_finite_velocity_diagnostic_reports_callback_context()
{
    const qmomkit::CqmomResult cqmom_result = make_fixed_two_dimensional_cqmom();
    assert(cqmom_result.ok());

    const qmomkit::CqmomNodeAdvectionResult result =
        qmomkit::advance_nodes(cqmom_result.quadrature, 0.125, 0.25, [](const double time, const std::vector<double>& xi)
                               { return std::vector<double>{xi[0] + time, -std::numeric_limits<double>::infinity()}; });

    assert(!result.ok());
    assert(result.status == qmomkit::NodeAdvectionStatus::velocity_failed);
    assert(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);
    assert(result.diagnostics.message == result.message);
    assert(result.quadrature.weights() == cqmom_result.quadrature.weights());
    assert(result.quadrature.abscissas() == cqmom_result.quadrature.abscissas());
    assert_contains(result.message, "Phase-space velocity callback returned a non-finite value:");
    assert_contains(result.message, "time: 0.125");
    assert_contains(result.message, "dimension: 1");
    assert_contains(result.message, "value: -inf");
    assert_contains(result.message, "state: [-1, -0.49999999999999972]");
    assert_contains(result.message, "velocity: [-0.875, -inf]");
}

void test_rk23_failure_diagnostics_are_reported()
{
    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {1.0};
    quadrature.weights = {1.0};

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::adaptive_rk23;
    options.rk23_options.relative_tolerance = 1.0e-12;
    options.rk23_options.absolute_tolerance = 1.0e-15;
    options.rk23_options.initial_step = 1.0e-3;
    options.rk23_options.max_steps = 1;

    const qmomkit::QmomNodeAdvectionResult result =
        qmomkit::advance_nodes(quadrature, 0.0, 1.0, [](double, const std::vector<double>& xi) { return std::vector<double>{xi[0]}; }, options);

    assert(!result.ok());
    assert(result.status == qmomkit::NodeAdvectionStatus::integration_failed);
    assert(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::adaptive_rk23);
    assert(result.diagnostics.adaptive_rk23_status == qmomkit::AdaptiveRk23Status::max_steps_exceeded);
    assert(result.diagnostics.step_count == 1);
    assert(result.quadrature.abscissas == quadrature.abscissas);
    assert(result.quadrature.weights == quadrature.weights);
}

void test_configured_rk23_options_reach_node_advection()
{
    const auto configured = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection]
scheme = adaptive_rk23

[time_integration.node_advection.rk23]
initial_step = 0.05
max_step = 0.05
relative_tolerance = 1.0e-6
absolute_tolerance = 1.0e-9
)");
    assert(configured.ok());

    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {0.0};
    quadrature.weights = {1.0};

    const qmomkit::QmomNodeAdvectionResult result =
        qmomkit::advance_nodes(quadrature, 0.0, 0.2, [](double, const std::vector<double>&) { return std::vector<double>{1.0}; }, configured.value());

    assert(result.ok());
    assert(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::adaptive_rk23);
    assert(result.diagnostics.adaptive_rk23_status == qmomkit::AdaptiveRk23Status::success);
    assert(result.diagnostics.step_count == 4);
    assert(result.diagnostics.rhs_evaluations == 13);
    assert_near(result.diagnostics.last_step, 0.05, 1.0e-15);
    assert_near(result.quadrature.abscissas[0], 0.2, 1.0e-15);
}

void test_node_advection_options_formatting_reports_semantic_settings()
{
    const std::string fixed = qmomkit::format_node_advection_options(qmomkit::NodeAdvectionOptions{});
    assert_contains(fixed, "scheme = ssp_rk3");
    assert(fixed.find("rk23.relative_tolerance") == std::string::npos);

    qmomkit::NodeAdvectionOptions adaptive_defaults;
    adaptive_defaults.scheme = qmomkit::NodeAdvectionScheme::adaptive_rk23;
    const std::string automatic = qmomkit::format_node_advection_options(adaptive_defaults);
    assert_contains(automatic, "rk23.initial_step = automatic");
    assert_contains(automatic, "rk23.max_step = full integration interval");
    assert_contains(automatic, "rk23.min_step = machine-dependent automatic minimum");

    qmomkit::NodeAdvectionOptions adaptive;
    adaptive.scheme = qmomkit::NodeAdvectionScheme::adaptive_rk23;
    adaptive.rk23_options.relative_tolerance = 0.125;
    adaptive.rk23_options.absolute_tolerance = 0.25;
    adaptive.rk23_options.initial_step = 0.5;
    adaptive.rk23_options.max_step = 0.75;
    adaptive.rk23_options.min_step = 0.0625;
    adaptive.rk23_options.max_steps = 37;

    const std::string formatted = qmomkit::format_node_advection_options(adaptive);
    assert_contains(formatted, "scheme = adaptive_rk23");
    assert_contains(formatted, "rk23.relative_tolerance = 0.125");
    assert_contains(formatted, "rk23.absolute_tolerance = 0.25");
    assert_contains(formatted, "rk23.initial_step = 0.5");
    assert_contains(formatted, "rk23.max_step = 0.75");
    assert_contains(formatted, "rk23.min_step = 0.0625");
    assert_contains(formatted, "rk23.max_steps = 37");
}

} // namespace

int main()
{
    test_adaptive_quadrature_structure_is_preserved();
    test_default_scheme_is_ssp_rk3();
    test_invalid_input_status_covers_quadrature_and_callback_contracts();
    test_wrong_size_velocity_diagnostic_reports_expected_and_returned_size();
    test_non_finite_velocity_diagnostic_reports_callback_context();
    test_rk23_failure_diagnostics_are_reported();
    test_configured_rk23_options_reach_node_advection();
    test_node_advection_options_formatting_reports_semantic_settings();
    return 0;
}
