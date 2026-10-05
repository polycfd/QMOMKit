#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
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
    return qmomkit::compute_cqmom(qmomkit::CqmomCase{qmomkit::extract_cqmom_moments(
        {extents, raw_moments_from_terminals(extents, weights, abscissas)}, qmomkit::create_cqmom_moment_layout(hierarchy))});
}

qmomkit::CqmomResult make_adaptive_two_dimensional_cqmom()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.1, 0.3, 0.6};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 1.0}, {-2.0, 2.0, 3.0}};
    return qmomkit::compute_cqmom(
        qmomkit::CqmomCase{qmomkit::extract_cqmom_moments({extents, raw_moments_from_terminals(extents, weights, abscissas)},
                                                          qmomkit::create_cqmom_moment_layout(hierarchy)),
                           {},
                           {true, 1.0e-10, 1.0e-10, 1.0e-12},
                           true});
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
    const qmomkit::CqmomNodeAdvectionResult result = qmomkit::advance_nodes(cqmom_result.quadrature, 0.0, dt,
                                                                            qmomkit::make_phase_space_velocity(
                                                                                [](double, const std::vector<double>& xi)
                                                                                {
                                                                                    assert(xi.size() == 2);
                                                                                    return std::vector<double>{rate_x * xi[0], rate_y * xi[1]};
                                                                                }),
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
    const qmomkit::QmomNodeAdvectionResult result = qmomkit::advance_nodes(
        quadrature, 2.0, 0.5, qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{3.0}; }),
        options);

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
    const qmomkit::QmomNodeAdvectionResult result = qmomkit::advance_nodes(
        quadrature, 0.0, time_step, qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{3.0}; }));

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
    const qmomkit::CqmomNodeAdvectionResult result =
        qmomkit::advance_nodes(cqmom_result.quadrature, 0.0, time_step,
                               qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{2.0, -3.0}; }));

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
    const qmomkit::QmomNodeAdvectionResult result = qmomkit::advance_nodes(
        quadrature, 0.0, time_step, qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{3.0}; }),
        options);

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
    const qmomkit::QmomNodeAdvectionResult result = qmomkit::advance_nodes(
        quadrature, 0.0, 0.25, qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{2.0}; }),
        options);

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
    const qmomkit::CqmomNodeAdvectionResult result = qmomkit::advance_nodes(cqmom_result.quadrature, 0.0, time_step,
                                                                            qmomkit::make_phase_space_velocity(
                                                                                [&velocity_calls](double, const std::vector<double>& coordinates)
                                                                                {
                                                                                    ++velocity_calls;
                                                                                    assert(coordinates.size() == 2);
                                                                                    return std::vector<double>{2.0, -3.0};
                                                                                }),
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
    const qmomkit::CqmomNodeAdvectionResult result = qmomkit::advance_nodes(cqmom_result.quadrature, 0.0, time_step,
                                                                            qmomkit::make_phase_space_velocity(
                                                                                [&velocity_calls](double, const std::vector<double>& coordinates)
                                                                                {
                                                                                    ++velocity_calls;
                                                                                    assert(coordinates.size() == 2);
                                                                                    return std::vector<double>{2.0, -3.0};
                                                                                }),
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

    const qmomkit::QmomNodeAdvectionResult invalid_qmom = qmomkit::advance_nodes(
        qmomkit::Quadrature1D{{1.0}, {1.0}, qmomkit::QuadratureStatus::not_realizable, "manufactured unsuccessful quadrature"}, 0.0, 0.1,
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{0.0}; }));
    assert(!invalid_qmom.ok());
    assert(invalid_qmom.status == qmomkit::NodeAdvectionStatus::invalid_input);
    assert_contains(invalid_qmom.message, "Input QMOM quadrature is not successful.");

    const qmomkit::CqmomResult valid_cqmom = make_fixed_two_dimensional_cqmom();
    const qmomkit::CqmomNodeAdvectionResult invalid_cqmom_input =
        qmomkit::advance_nodes(valid_cqmom.quadrature, std::numeric_limits<double>::quiet_NaN(), 0.1,
                               qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{0.0, 0.0}; }));
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
    const qmomkit::QmomNodeAdvectionResult result = qmomkit::advance_nodes(
        quadrature, 0.0, 0.25, qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{1.0, 2.0}; }),
        options);

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
        cqmom_result.quadrature, 0.125, 0.25,
        qmomkit::make_phase_space_velocity([](const double time, const std::vector<double>& xi)
                                           { return std::vector<double>{xi[0] + time, -std::numeric_limits<double>::infinity()}; }),
        options);

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

// A velocity written directly in packed form and the same physics reached
// through make_phase_space_velocity must agree exactly for every scheme, not
// merely closely.
void test_adapted_velocity_matches_directly_written_velocity_exactly()
{
    const qmomkit::CqmomResult cqmom_result = make_fixed_two_dimensional_cqmom();
    assert(cqmom_result.ok());

    const qmomkit::PerNodePhaseSpaceVelocity per_node = [](const double time, const std::vector<double>& xi)
    { return std::vector<double>{0.3 * xi[0] - 0.1 * time * xi[1], 0.2 * xi[1] + 0.05 * xi[0]}; };
    const qmomkit::PhaseSpaceVelocity packed =
        [](const double time, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
    {
        for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
        {
            const double x = coordinates[node * dimensions];
            const double y = coordinates[node * dimensions + 1];
            velocities[node * dimensions] = 0.3 * x - 0.1 * time * y;
            velocities[node * dimensions + 1] = 0.2 * y + 0.05 * x;
        }
    };

    for (const qmomkit::NodeAdvectionScheme scheme :
         {qmomkit::NodeAdvectionScheme::rk23, qmomkit::NodeAdvectionScheme::rk45, qmomkit::NodeAdvectionScheme::ssp_rk3})
    {
        qmomkit::NodeAdvectionOptions options;
        options.scheme = scheme;

        const qmomkit::CqmomNodeAdvectionResult from_per_node =
            qmomkit::advance_nodes(cqmom_result.quadrature, 0.125, 0.25, qmomkit::make_phase_space_velocity(per_node), options);
        const qmomkit::CqmomNodeAdvectionResult from_packed = qmomkit::advance_nodes(cqmom_result.quadrature, 0.125, 0.25, packed, options);

        assert(from_per_node.ok());
        assert(from_packed.ok());
        assert(from_packed.quadrature.weights() == from_per_node.quadrature.weights());
        assert(from_packed.quadrature.abscissas() == from_per_node.quadrature.abscissas());
        assert(from_packed.quadrature.terminal_paths() == from_per_node.quadrature.terminal_paths());
        assert(from_packed.diagnostics.rhs_evaluations == from_per_node.diagnostics.rhs_evaluations);
        assert(from_packed.diagnostics.step_count == from_per_node.diagnostics.step_count);
    }
}

void test_velocity_adapter_preserves_results_and_emptiness()
{
    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {0.5, 1.5};
    quadrature.weights = {0.25, 0.75};

    const qmomkit::PerNodePhaseSpaceVelocity per_node = [](double, const std::vector<double>& xi) { return std::vector<double>{-0.4 * xi[0]}; };

    const qmomkit::QmomNodeAdvectionResult direct =
        qmomkit::advance_nodes(quadrature, 0.0, 0.1,
                               [](double, const std::vector<double>& coordinates, std::vector<double>& velocities, std::size_t)
                               {
                                   for (std::size_t entry = 0; entry < coordinates.size(); ++entry)
                                   {
                                       velocities[entry] = -0.4 * coordinates[entry];
                                   }
                               });
    const qmomkit::QmomNodeAdvectionResult adapted = qmomkit::advance_nodes(quadrature, 0.0, 0.1, qmomkit::make_phase_space_velocity(per_node));
    assert(direct.ok());
    assert(adapted.ok());
    assert(adapted.quadrature.abscissas == direct.quadrature.abscissas);
    assert(adapted.quadrature.weights == direct.quadrature.weights);

    // An empty per-node callback must adapt to an empty packed callback so the
    // missing-velocity contract still reports invalid input.
    const qmomkit::PhaseSpaceVelocity adapted_empty = qmomkit::make_phase_space_velocity(qmomkit::PerNodePhaseSpaceVelocity{});
    assert(!adapted_empty);
    const qmomkit::QmomNodeAdvectionResult missing = qmomkit::advance_nodes(quadrature, 0.0, 0.1, adapted_empty);
    assert(!missing.ok());
    assert(missing.status == qmomkit::NodeAdvectionStatus::invalid_input);
    assert_contains(missing.message, "velocity callback must be provided");
}

void test_velocity_diagnostics_report_node_and_dimension_context()
{
    const qmomkit::CqmomResult cqmom_result = make_fixed_two_dimensional_cqmom();
    assert(cqmom_result.ok());

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    // A non-finite entry must still be attributed to its phase-space dimension
    // with the offending node's coordinate and velocity slices.
    const qmomkit::CqmomNodeAdvectionResult non_finite = qmomkit::advance_nodes(
        cqmom_result.quadrature, 0.125, 0.25,
        [](const double time, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
        {
            for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
            {
                velocities[node * dimensions] = coordinates[node * dimensions] + time;
                velocities[node * dimensions + 1] = -std::numeric_limits<double>::infinity();
            }
        },
        options);
    assert(!non_finite.ok());
    assert(non_finite.status == qmomkit::NodeAdvectionStatus::velocity_failed);
    assert_contains(non_finite.message, "Phase-space velocity callback returned a non-finite value:");
    assert_contains(non_finite.message, "time: 0.125");
    assert_contains(non_finite.message, "dimension: 1");
    assert_contains(non_finite.message, "value: -inf");
    assert_diagnostic_vector_near(non_finite.message, "state:", {-1.0, -0.5}, 1.0e-12);
    assert_diagnostic_vector_near(non_finite.message, "velocity:", {-0.875, -std::numeric_limits<double>::infinity()}, 1.0e-12);

    // Resizing the supplied buffer breaks the packed contract and is reported
    // as a dimension mismatch rather than corrupting the integrator state.
    const qmomkit::CqmomNodeAdvectionResult resized = qmomkit::advance_nodes(
        cqmom_result.quadrature, 0.0, 0.25,
        [](double, const std::vector<double>&, std::vector<double>& velocities, std::size_t) { velocities.assign(3, 0.0); }, options);
    assert(!resized.ok());
    assert(resized.status == qmomkit::NodeAdvectionStatus::velocity_failed);
    assert_contains(resized.message, "Phase-space velocity callback returned the wrong dimension:");
    assert(resized.quadrature.abscissas() == cqmom_result.quadrature.abscissas());
}

// The packed callback fills a library-supplied buffer, so a callback that
// writes only part of it cannot be caught by size or finiteness alone. Debug
// builds poison the buffer to detect it; in release builds the check is
// compiled out and, like every assert-based test here, this one is inert.
// The adapter hands the callback a reference to its coordinate buffer. A
// callback that re-enters node advection on the same thread must not be able to
// resize that buffer underneath the reference its caller is still holding: with
// shared storage the re-entrant call reallocates it, so the outer callback's
// argument moves in memory and its contents are replaced mid-call.
void test_reentrant_velocity_callback_keeps_its_coordinate_argument()
{
    // The inner system is wider than the outer one, so a shared buffer has to
    // grow, which is what makes the reallocation observable.
    const qmomkit::PhaseSpaceVelocity inner = qmomkit::make_phase_space_velocity([](const double, const std::vector<double>& coordinate)
                                                                                 { return std::vector<double>(coordinate.size(), 0.0); });

    bool reentered = false;
    bool argument_survived = true;
    const qmomkit::PhaseSpaceVelocity outer = qmomkit::make_phase_space_velocity(
        [&inner, &reentered, &argument_survived](const double time, const std::vector<double>& coordinate)
        {
            const double* const held = coordinate.data();
            const std::size_t held_size = coordinate.size();
            if (!reentered)
            {
                reentered = true;
                std::vector<double> inner_coordinates(3 * 5, 1.0);
                std::vector<double> inner_velocities(3 * 5, 0.0);
                inner(time, inner_coordinates, inner_velocities, 5);
            }
            argument_survived = argument_survived && coordinate.data() == held && coordinate.size() == held_size && coordinate.size() == 2 &&
                                coordinate[0] == 4.0 && coordinate[1] == 9.0;
            return std::vector<double>(coordinate.size(), 0.0);
        });

    std::vector<double> coordinates{4.0, 9.0};
    std::vector<double> velocities(coordinates.size(), 0.0);
    outer(0.0, coordinates, velocities, 2);

    assert(reentered);
    assert(argument_survived);
}

void test_partially_filled_velocity_is_reported_in_debug_builds()
{
    const qmomkit::CqmomResult cqmom_result = make_fixed_two_dimensional_cqmom();
    assert(cqmom_result.ok());

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    // Writes the first phase-space dimension of every node and forgets the second.
    const qmomkit::CqmomNodeAdvectionResult partial = qmomkit::advance_nodes(
        cqmom_result.quadrature, 0.125, 0.25,
        [](const double time, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
        {
            for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
            {
                velocities[node * dimensions] = coordinates[node * dimensions] + time;
            }
        },
        options);

    assert(!partial.ok());
    assert(partial.status == qmomkit::NodeAdvectionStatus::velocity_failed);
    assert(partial.diagnostics.message == partial.message);
    assert_contains(partial.message, "Phase-space velocity callback did not fill every velocity entry:");
    assert_contains(partial.message, "time: 0.125");
    assert_contains(partial.message, "node: 0");
    assert_contains(partial.message, "dimension: 1");
    assert_diagnostic_vector_near(partial.message, "coordinate:", {-1.0, -0.5}, 1.0e-12);
    assert(partial.quadrature.abscissas() == cqmom_result.quadrature.abscissas());
    assert(partial.quadrature.weights() == cqmom_result.quadrature.weights());

    // A callback that fills every entry must remain unaffected by the poison.
    const qmomkit::CqmomNodeAdvectionResult complete = qmomkit::advance_nodes(
        cqmom_result.quadrature, 0.125, 0.25,
        [](const double time, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
        {
            for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
            {
                velocities[node * dimensions] = coordinates[node * dimensions] + time;
                velocities[node * dimensions + 1] = 0.0;
            }
        },
        options);
    assert(complete.ok());
}

// One adapted velocity object may be shared by an application-managed parallel
// loop. The adapter reuses a coordinate buffer across nodes and right-hand-side
// evaluations, so that buffer must not be shared between threads.
void test_adapted_velocity_is_safe_to_share_between_threads()
{
    const qmomkit::CqmomResult cqmom_result = make_fixed_two_dimensional_cqmom();
    assert(cqmom_result.ok());

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    // Dimension-dependent rates, so a buffer shared between threads would mix
    // coordinates between nodes and change the advected result.
    const qmomkit::PhaseSpaceVelocity shared = qmomkit::make_phase_space_velocity(
        [](const double time, const std::vector<double>& coordinate)
        {
            std::vector<double> rate(coordinate.size(), 0.0);
            for (std::size_t dimension = 0; dimension < coordinate.size(); ++dimension)
            {
                rate[dimension] = 0.25 * coordinate[dimension] + 0.1 * static_cast<double>(dimension + 1) - 0.05 * time;
            }
            return rate;
        });

    const qmomkit::CqmomNodeAdvectionResult reference = qmomkit::advance_nodes(cqmom_result.quadrature, 0.125, 0.25, shared, options);
    assert(reference.ok());

    std::vector<int> succeeded(4, 0);
    std::vector<std::thread> workers;
    workers.reserve(succeeded.size());
    for (std::size_t worker = 0; worker < succeeded.size(); ++worker)
    {
        workers.emplace_back(
            [&cqmom_result, &shared, &options, &reference, &succeeded, worker]()
            {
                for (std::size_t iteration = 0; iteration < 32; ++iteration)
                {
                    const qmomkit::CqmomNodeAdvectionResult result = qmomkit::advance_nodes(cqmom_result.quadrature, 0.125, 0.25, shared, options);
                    assert(result.ok());
                    assert(result.quadrature.abscissas() == reference.quadrature.abscissas());
                    assert(result.quadrature.weights() == reference.quadrature.weights());
                }
                succeeded[worker] = 1;
            });
    }
    for (std::thread& worker : workers)
    {
        worker.join();
    }
    assert(std::all_of(succeeded.begin(), succeeded.end(), [](const int value) { return value == 1; }));
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

    const qmomkit::QmomNodeAdvectionResult result = qmomkit::advance_nodes(
        quadrature, 0.0, 1.0, qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>{xi[0]}; }),
        options);

    assert(!result.ok());
    assert(result.status == qmomkit::NodeAdvectionStatus::integration_failed);
    assert(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert(result.diagnostics.adaptive_status == qmomkit::AdaptiveRkStatus::max_steps_exceeded);
    assert(result.diagnostics.step_count == 1);
    assert(result.quadrature.abscissas == quadrature.abscissas);
    assert(result.quadrature.weights == quadrature.weights);
}

// A velocity callback that throws stops the integration before the integrator
// returns, so there is no adaptive status to report. The evaluations made
// before the failure are still known and reported, on both scheme families.
void test_throwing_velocity_reports_its_evaluations_without_an_adaptive_status()
{
    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {1.0, 2.0};
    quadrature.weights = {0.5, 0.5};

    for (const qmomkit::NodeAdvectionScheme scheme :
         {qmomkit::NodeAdvectionScheme::rk23, qmomkit::NodeAdvectionScheme::rk45, qmomkit::NodeAdvectionScheme::ssp_rk3})
    {
        std::size_t calls = 0;
        const qmomkit::PhaseSpaceVelocity failing_velocity =
            [&calls](double, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t)
        {
            ++calls;
            if (calls == 3)
            {
                throw std::runtime_error("velocity failure on the third evaluation");
            }
            for (std::size_t i = 0; i < coordinates.size(); ++i)
            {
                velocities[i] = -coordinates[i];
            }
        };
        qmomkit::NodeAdvectionOptions options;
        options.scheme = scheme;

        const qmomkit::QmomNodeAdvectionResult result = qmomkit::advance_nodes(quadrature, 0.25, 0.5, failing_velocity, options);

        assert(result.status == qmomkit::NodeAdvectionStatus::velocity_failed);
        assert(result.diagnostics.scheme == scheme);
        assert(!result.diagnostics.adaptive_status.has_value());
        assert(result.diagnostics.rhs_evaluations == 3);
        assert(result.diagnostics.final_time == 0.25);
        assert(result.quadrature.abscissas == quadrature.abscissas);
    }
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

    const qmomkit::QmomNodeAdvectionResult result = qmomkit::advance_nodes(
        quadrature, 0.0, 0.2, qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{1.0}; }),
        configured.value());

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

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"adaptive_quadrature_structure_is_preserved", test_adaptive_quadrature_structure_is_preserved},
            {"default_scheme_is_rk23", test_default_scheme_is_rk23},
            {"default_qmom_rk23_uses_default_adaptive_initial_step", test_default_qmom_rk23_uses_default_adaptive_initial_step},
            {"default_cqmom_rk23_uses_default_adaptive_initial_step", test_default_cqmom_rk23_uses_default_adaptive_initial_step},
            {"rk45_uses_default_adaptive_initial_step_when_omitted", test_rk45_uses_default_adaptive_initial_step_when_omitted},
            {"ssp_rk3_ignores_embedded_initial_step", test_ssp_rk3_ignores_embedded_initial_step},
            {"cqmom_rk45_preserves_fixed_topology_and_advects_coordinates", test_cqmom_rk45_preserves_fixed_topology_and_advects_coordinates},
            {"cqmom_ssp_rk3_preserves_fixed_topology_and_advects_coordinates", test_cqmom_ssp_rk3_preserves_fixed_topology_and_advects_coordinates},
            {"invalid_input_status_covers_quadrature_and_callback_contracts", test_invalid_input_status_covers_quadrature_and_callback_contracts},
            {"wrong_size_velocity_diagnostic_reports_expected_and_returned_size",
             test_wrong_size_velocity_diagnostic_reports_expected_and_returned_size},
            {"non_finite_velocity_diagnostic_reports_callback_context", test_non_finite_velocity_diagnostic_reports_callback_context},
            {"adapted_velocity_matches_directly_written_velocity_exactly", test_adapted_velocity_matches_directly_written_velocity_exactly},
            {"velocity_adapter_preserves_results_and_emptiness", test_velocity_adapter_preserves_results_and_emptiness},
            {"velocity_diagnostics_report_node_and_dimension_context", test_velocity_diagnostics_report_node_and_dimension_context},
            {"adapted_velocity_is_safe_to_share_between_threads", test_adapted_velocity_is_safe_to_share_between_threads},
            {"reentrant_velocity_callback_keeps_its_coordinate_argument", test_reentrant_velocity_callback_keeps_its_coordinate_argument},
            {"partially_filled_velocity_is_reported_in_debug_builds", test_partially_filled_velocity_is_reported_in_debug_builds},
            {"rk23_failure_diagnostics_are_reported", test_rk23_failure_diagnostics_are_reported},
            {"throwing_velocity_reports_its_evaluations_without_an_adaptive_status",
             test_throwing_velocity_reports_its_evaluations_without_an_adaptive_status},
            {"configured_adaptive_options_reach_node_advection", test_configured_adaptive_options_reach_node_advection},
            {"node_advection_options_formatting_reports_semantic_settings", test_node_advection_options_formatting_reports_semantic_settings},
        });
}
