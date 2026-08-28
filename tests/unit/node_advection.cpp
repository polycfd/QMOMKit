#include "qmomkit/qmomkit.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <sstream>
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

// Diagnostic vectors are printed with the full round-trippable precision of a
// double, so their digits depend on the toolchain. Parse the formatted list and
// compare the values numerically instead of matching the text verbatim.
std::vector<double> parse_diagnostic_vector(const std::string& text, const std::string& label)
{
    const std::size_t label_position = text.find(label);
    assert(label_position != std::string::npos);
    const std::size_t open = text.find('[', label_position);
    assert(open != std::string::npos);
    const std::size_t close = text.find(']', open);
    assert(close != std::string::npos);

    std::vector<double> values;
    std::istringstream entries(text.substr(open + 1, close - open - 1));
    std::string entry;
    while (std::getline(entries, entry, ','))
    {
        values.push_back(std::stod(entry));
    }
    return values;
}

void assert_diagnostic_vector_near(const std::string& text, const std::string& label, const std::vector<double>& expected, const double tolerance)
{
    const std::vector<double> values = parse_diagnostic_vector(text, label);
    assert(values.size() == expected.size());
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (std::isnan(expected[i]))
        {
            assert(std::isnan(values[i]));
        }
        else if (std::isinf(expected[i]))
        {
            assert(values[i] == expected[i]);
        }
        else
        {
            assert_near(values[i], expected[i], tolerance);
        }
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
    options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    options.adaptive_options.relative_tolerance = 1.0e-8;
    options.adaptive_options.absolute_tolerance = 1.0e-11;
    options.adaptive_options.initial_step = 0.01;

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

void test_default_scheme_is_rk23()
{
    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {1.0, 2.0};
    quadrature.weights = {0.25, 0.75};

    qmomkit::NodeAdvectionOptions options;
    options.adaptive_options.initial_step = 0.5;
    const qmomkit::QmomNodeAdvectionResult result =
        qmomkit::advance_nodes(quadrature, 2.0, 0.5, [](double, const std::vector<double>&) { return std::vector<double>{3.0}; }, options);

    assert(result.ok());
    assert(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert(result.diagnostics.step_count == 1);
    assert(result.diagnostics.rejected_step_count == 0);
    assert(result.diagnostics.rhs_evaluations == 4);
    assert(result.diagnostics.adaptive_status == qmomkit::AdaptiveRkStatus::success);
    assert_near(result.diagnostics.final_time, 2.5, 1.0e-15);
    assert_near(result.quadrature.abscissas[0], 2.5, 1.0e-15);
    assert_near(result.quadrature.abscissas[1], 3.5, 1.0e-15);
    assert(result.quadrature.weights == quadrature.weights);
}

void test_default_qmom_rk23_uses_default_adaptive_initial_step()
{
    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {1.0};
    quadrature.weights = {1.0};

    const double time_step = 2.0 * qmomkit::default_adaptive_initial_step;
    const qmomkit::QmomNodeAdvectionResult result =
        qmomkit::advance_nodes(quadrature, 0.0, time_step, [](double, const std::vector<double>&) { return std::vector<double>{3.0}; });

    assert(result.ok());
    assert(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert(result.diagnostics.step_count == 2);
    assert(result.diagnostics.rhs_evaluations == 7);
    assert_near(result.quadrature.abscissas[0], 1.0 + 3.0 * time_step, 1.0e-15);
}

void test_default_cqmom_rk23_uses_default_adaptive_initial_step()
{
    const qmomkit::CqmomResult cqmom_result = make_fixed_two_dimensional_cqmom();
    assert(cqmom_result.ok());

    const double time_step = 2.0 * qmomkit::default_adaptive_initial_step;
    const qmomkit::CqmomNodeAdvectionResult result = qmomkit::advance_nodes(
        cqmom_result.quadrature, 0.0, time_step, [](double, const std::vector<double>&) { return std::vector<double>{2.0, -3.0}; });

    assert(result.ok());
    assert(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert(result.diagnostics.step_count == 2);
    assert(result.diagnostics.rhs_evaluations == 7);
    for (std::size_t node = 0; node < result.quadrature.total_nodes(); ++node)
    {
        assert_near(result.quadrature.abscissas()[0][node], cqmom_result.quadrature.abscissas()[0][node] + 2.0 * time_step, 1.0e-15);
        assert_near(result.quadrature.abscissas()[1][node], cqmom_result.quadrature.abscissas()[1][node] - 3.0 * time_step, 1.0e-15);
    }
}

void test_rk45_uses_default_adaptive_initial_step_when_omitted()
{
    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {1.0};
    quadrature.weights = {1.0};

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::rk45;
    const double time_step = 2.0 * qmomkit::default_adaptive_initial_step;
    const qmomkit::QmomNodeAdvectionResult result =
        qmomkit::advance_nodes(quadrature, 0.0, time_step, [](double, const std::vector<double>&) { return std::vector<double>{3.0}; }, options);

    assert(result.ok());
    assert(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk45);
    assert(result.diagnostics.step_count == 2);
    assert(result.diagnostics.rhs_evaluations == 13);
    assert_near(result.quadrature.abscissas[0], 1.0 + 3.0 * time_step, 1.0e-15);
}

void test_ssp_rk3_ignores_embedded_initial_step()
{
    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {1.0};
    quadrature.weights = {1.0};

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    options.adaptive_options.initial_step = 1.0e-12;
    const qmomkit::QmomNodeAdvectionResult result =
        qmomkit::advance_nodes(quadrature, 0.0, 0.25, [](double, const std::vector<double>&) { return std::vector<double>{2.0}; }, options);

    assert(result.ok());
    assert(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);
    assert(result.diagnostics.rhs_evaluations == 3);
    assert_near(result.quadrature.abscissas[0], 1.5, 1.0e-15);
}

void test_cqmom_rk45_preserves_fixed_topology_and_advects_coordinates()
{
    const qmomkit::CqmomResult cqmom_result = make_fixed_two_dimensional_cqmom();
    assert(cqmom_result.ok());

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::rk45;
    const double time_step = 2.0 * qmomkit::default_adaptive_initial_step;
    std::size_t velocity_calls = 0;
    const qmomkit::CqmomNodeAdvectionResult result = qmomkit::advance_nodes(
        cqmom_result.quadrature, 0.0, time_step,
        [&velocity_calls](double, const std::vector<double>& coordinates)
        {
            ++velocity_calls;
            assert(coordinates.size() == 2);
            return std::vector<double>{2.0, -3.0};
        },
        options);

    assert(result.ok());
    assert(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk45);
    assert(result.diagnostics.adaptive_status == qmomkit::AdaptiveRkStatus::success);
    assert(result.diagnostics.step_count == 2);
    assert(result.diagnostics.rhs_evaluations == 13);
    assert(velocity_calls == result.diagnostics.rhs_evaluations * cqmom_result.quadrature.total_nodes());
    assert(result.quadrature.hierarchy().node_counts() == cqmom_result.quadrature.hierarchy().node_counts());
    assert(result.quadrature.has_rectangular_storage());
    assert(result.quadrature.terminal_paths() == cqmom_result.quadrature.terminal_paths());
    assert(result.quadrature.weights() == cqmom_result.quadrature.weights());
    for (std::size_t node = 0; node < result.quadrature.total_nodes(); ++node)
    {
        assert_near(result.quadrature.abscissas()[0][node], cqmom_result.quadrature.abscissas()[0][node] + 2.0 * time_step, 1.0e-15);
        assert_near(result.quadrature.abscissas()[1][node], cqmom_result.quadrature.abscissas()[1][node] - 3.0 * time_step, 1.0e-15);
    }
}

void test_cqmom_ssp_rk3_preserves_fixed_topology_and_advects_coordinates()
{
    const qmomkit::CqmomResult cqmom_result = make_fixed_two_dimensional_cqmom();
    assert(cqmom_result.ok());

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    options.adaptive_options.initial_step = 1.0e-12;
    constexpr double time_step = 0.25;
    std::size_t velocity_calls = 0;
    const qmomkit::CqmomNodeAdvectionResult result = qmomkit::advance_nodes(
        cqmom_result.quadrature, 0.0, time_step,
        [&velocity_calls](double, const std::vector<double>& coordinates)
        {
            ++velocity_calls;
            assert(coordinates.size() == 2);
            return std::vector<double>{2.0, -3.0};
        },
        options);

    assert(result.ok());
    assert(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);
    assert(result.diagnostics.step_count == 1);
    assert(result.diagnostics.rhs_evaluations == 3);
    assert(!result.diagnostics.adaptive_status.has_value());
    assert(velocity_calls == result.diagnostics.rhs_evaluations * cqmom_result.quadrature.total_nodes());
    assert(result.quadrature.hierarchy().node_counts() == cqmom_result.quadrature.hierarchy().node_counts());
    assert(result.quadrature.has_rectangular_storage());
    assert(result.quadrature.terminal_paths() == cqmom_result.quadrature.terminal_paths());
    assert(result.quadrature.weights() == cqmom_result.quadrature.weights());
    for (std::size_t node = 0; node < result.quadrature.total_nodes(); ++node)
    {
        assert_near(result.quadrature.abscissas()[0][node], cqmom_result.quadrature.abscissas()[0][node] + 2.0 * time_step, 1.0e-15);
        assert_near(result.quadrature.abscissas()[1][node], cqmom_result.quadrature.abscissas()[1][node] - 3.0 * time_step, 1.0e-15);
    }
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

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::QmomNodeAdvectionResult result =
        qmomkit::advance_nodes(quadrature, 0.0, 0.25, [](double, const std::vector<double>&) { return std::vector<double>{1.0, 2.0}; }, options);

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

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomNodeAdvectionResult result = qmomkit::advance_nodes(
        cqmom_result.quadrature, 0.125, 0.25, [](const double time, const std::vector<double>& xi)
        { return std::vector<double>{xi[0] + time, -std::numeric_limits<double>::infinity()}; }, options);

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
    assert_diagnostic_vector_near(result.message, "state:", {-1.0, -0.5}, 1.0e-12);
    assert_diagnostic_vector_near(result.message, "velocity:", {-0.875, -std::numeric_limits<double>::infinity()}, 1.0e-12);
}

void test_rk23_failure_diagnostics_are_reported()
{
    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {1.0};
    quadrature.weights = {1.0};

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    options.adaptive_options.relative_tolerance = 1.0e-12;
    options.adaptive_options.absolute_tolerance = 1.0e-15;
    options.adaptive_options.initial_step = 1.0e-3;
    options.adaptive_options.max_steps = 1;

    const qmomkit::QmomNodeAdvectionResult result =
        qmomkit::advance_nodes(quadrature, 0.0, 1.0, [](double, const std::vector<double>& xi) { return std::vector<double>{xi[0]}; }, options);

    assert(!result.ok());
    assert(result.status == qmomkit::NodeAdvectionStatus::integration_failed);
    assert(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert(result.diagnostics.adaptive_status == qmomkit::AdaptiveRkStatus::max_steps_exceeded);
    assert(result.diagnostics.step_count == 1);
    assert(result.quadrature.abscissas == quadrature.abscissas);
    assert(result.quadrature.weights == quadrature.weights);
}

void test_configured_adaptive_options_reach_node_advection()
{
    const auto configured = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection]
scheme = rk23

[time_integration.node_advection.adaptive]
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
    assert(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert(result.diagnostics.adaptive_status == qmomkit::AdaptiveRkStatus::success);
    assert(result.diagnostics.step_count == 4);
    assert(result.diagnostics.rhs_evaluations == 13);
    assert_near(result.diagnostics.last_step, 0.05, 1.0e-15);
    assert_near(result.quadrature.abscissas[0], 0.2, 1.0e-15);
}

void test_node_advection_options_formatting_reports_semantic_settings()
{
    const std::string fixed = qmomkit::format_node_advection_options(qmomkit::NodeAdvectionOptions{});
    assert_contains(fixed, "scheme = rk23");
    assert_contains(fixed, "adaptive.initial_step = ");
    assert(fixed.find("<required>") == std::string::npos);

    qmomkit::NodeAdvectionOptions adaptive_defaults;
    adaptive_defaults.scheme = qmomkit::NodeAdvectionScheme::rk23;
    const std::string automatic = qmomkit::format_node_advection_options(adaptive_defaults);
    assert_contains(automatic, "adaptive.initial_step = ");
    assert(automatic.find("<required>") == std::string::npos);
    assert_contains(automatic, "adaptive.max_step = full integration interval");
    assert_contains(automatic, "adaptive.min_step = machine-dependent automatic minimum");

    qmomkit::NodeAdvectionOptions adaptive;
    adaptive.scheme = qmomkit::NodeAdvectionScheme::rk23;
    adaptive.adaptive_options.relative_tolerance = 0.125;
    adaptive.adaptive_options.absolute_tolerance = 0.25;
    adaptive.adaptive_options.initial_step = 0.5;
    adaptive.adaptive_options.max_step = 0.75;
    adaptive.adaptive_options.min_step = 0.0625;
    adaptive.adaptive_options.max_steps = 37;

    const std::string formatted = qmomkit::format_node_advection_options(adaptive);
    assert_contains(formatted, "scheme = rk23");
    assert_contains(formatted, "adaptive.relative_tolerance = 0.125");
    assert_contains(formatted, "adaptive.absolute_tolerance = 0.25");
    assert_contains(formatted, "adaptive.initial_step = 0.5");
    assert_contains(formatted, "adaptive.max_step = 0.75");
    assert_contains(formatted, "adaptive.min_step = 0.0625");
    assert_contains(formatted, "adaptive.max_steps = 37");
}

} // namespace

int main()
{
    test_adaptive_quadrature_structure_is_preserved();
    test_default_scheme_is_rk23();
    test_default_qmom_rk23_uses_default_adaptive_initial_step();
    test_default_cqmom_rk23_uses_default_adaptive_initial_step();
    test_rk45_uses_default_adaptive_initial_step_when_omitted();
    test_ssp_rk3_ignores_embedded_initial_step();
    test_cqmom_rk45_preserves_fixed_topology_and_advects_coordinates();
    test_cqmom_ssp_rk3_preserves_fixed_topology_and_advects_coordinates();
    test_invalid_input_status_covers_quadrature_and_callback_contracts();
    test_wrong_size_velocity_diagnostic_reports_expected_and_returned_size();
    test_non_finite_velocity_diagnostic_reports_callback_context();
    test_rk23_failure_diagnostics_are_reported();
    test_configured_adaptive_options_reach_node_advection();
    test_node_advection_options_formatting_reports_semantic_settings();
    return 0;
}
