#include "../verification/verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

void assert_near(const double actual, const double expected, const double tolerance)
{
    QMOMKIT_VERIFY_NEAR(actual, expected, tolerance);
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

void test_constant_phase_space_velocity_for_qmom()
{
    const std::vector<double> input_abscissas{-1.0, 2.0};
    const std::vector<double> input_weights{0.25, 0.75};
    const qmomkit::QmomResult quadrature = qmomkit::compute_qmom(qmomkit::MomentSet(raw_moments_1d(input_abscissas, input_weights, 4)), 2);
    QMOMKIT_VERIFY_CHECK(quadrature.ok());

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::QmomNodeAdvectionResult result = qmomkit::advance_nodes(quadrature, 0.0, 0.25,
                                                                           qmomkit::make_phase_space_velocity(
                                                                               [](double, const std::vector<double>& xi)
                                                                               {
                                                                                   QMOMKIT_VERIFY_CHECK(xi.size() == 1);
                                                                                   return std::vector<double>{2.0};
                                                                               }),
                                                                           options);

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.step_count == 1);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.rhs_evaluations == 3);
    QMOMKIT_VERIFY_CHECK(!result.diagnostics.estimated_error.has_value());
    for (std::size_t node = 0; node < input_weights.size(); ++node)
    {
        assert_near(result.quadrature.weights[node], quadrature.weights[node], 1.0e-14);
    }
    assert_near(result.quadrature.abscissas[0], -0.5, 1.0e-13);
    assert_near(result.quadrature.abscissas[1], 2.5, 1.0e-13);
}

void test_linear_phase_space_velocity_for_qmom()
{
    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {-0.5, 1.25};
    quadrature.weights = {0.4, 0.6};

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    options.adaptive_options.relative_tolerance = 1.0e-8;
    options.adaptive_options.absolute_tolerance = 1.0e-11;
    options.adaptive_options.initial_step = 0.05;

    constexpr double rate = -0.75;
    constexpr double dt = 0.6;
    const qmomkit::QmomNodeAdvectionResult result = qmomkit::advance_nodes(
        quadrature, 0.0, dt,
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>{rate * xi[0]}; }), options);

    QMOMKIT_VERIFY_CHECK(result.ok());
    assert_near(result.quadrature.abscissas[0], quadrature.abscissas[0] * std::exp(rate * dt), 2.0e-8);
    assert_near(result.quadrature.abscissas[1], quadrature.abscissas[1] * std::exp(rate * dt), 2.0e-8);
    QMOMKIT_VERIFY_CHECK(result.quadrature.weights == quadrature.weights);
}

void test_multi_dimensional_phase_space_evolution()
{
    const qmomkit::CqmomResult cqmom_result = make_fixed_two_dimensional_cqmom();
    QMOMKIT_VERIFY_CHECK(cqmom_result.ok());

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    options.adaptive_options.relative_tolerance = 1.0e-8;
    options.adaptive_options.absolute_tolerance = 1.0e-11;
    options.adaptive_options.initial_step = 0.05;

    constexpr double rate_x = 0.3;
    constexpr double rate_y = -0.4;
    constexpr double dt = 0.5;
    const qmomkit::CqmomNodeAdvectionResult result = qmomkit::advance_nodes(cqmom_result.quadrature, 1.0, dt,
                                                                            qmomkit::make_phase_space_velocity(
                                                                                [](double, const std::vector<double>& xi)
                                                                                {
                                                                                    QMOMKIT_VERIFY_CHECK(xi.size() == 2);
                                                                                    return std::vector<double>{rate_x * xi[0], rate_y * xi[1]};
                                                                                }),
                                                                            options);

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.quadrature.has_rectangular_storage());
    QMOMKIT_VERIFY_CHECK(result.quadrature.weights() == cqmom_result.quadrature.weights());
    for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(cqmom_result.quadrature.hierarchy()))
    {
        assert_near(result.quadrature.abscissa(0, node_index), cqmom_result.quadrature.abscissa(0, node_index) * std::exp(rate_x * dt), 2.0e-8);
        assert_near(result.quadrature.abscissa(1, node_index), cqmom_result.quadrature.abscissa(1, node_index) * std::exp(rate_y * dt), 2.0e-8);
    }
}

void test_ssp_rk3_multi_dimensional_phase_space_evolution()
{
    const qmomkit::CqmomResult cqmom_result = make_fixed_two_dimensional_cqmom();
    QMOMKIT_VERIFY_CHECK(cqmom_result.ok());

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    constexpr double rate_x = 0.3;
    constexpr double rate_y = -0.4;
    constexpr double dt = 0.03;
    const qmomkit::CqmomNodeAdvectionResult result = qmomkit::advance_nodes(cqmom_result.quadrature, 1.0, dt,
                                                                            qmomkit::make_phase_space_velocity(
                                                                                [](double, const std::vector<double>& xi)
                                                                                {
                                                                                    QMOMKIT_VERIFY_CHECK(xi.size() == 2);
                                                                                    return std::vector<double>{rate_x * xi[0], rate_y * xi[1]};
                                                                                }),
                                                                            options);

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.step_count == 1);
    QMOMKIT_VERIFY_CHECK(result.quadrature.has_rectangular_storage());
    QMOMKIT_VERIFY_CHECK(result.quadrature.weights() == cqmom_result.quadrature.weights());
    for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(cqmom_result.quadrature.hierarchy()))
    {
        assert_near(result.quadrature.abscissa(0, node_index), cqmom_result.quadrature.abscissa(0, node_index) * std::exp(rate_x * dt), 3.0e-9);
        assert_near(result.quadrature.abscissa(1, node_index), cqmom_result.quadrature.abscissa(1, node_index) * std::exp(rate_y * dt), 3.0e-9);
    }
}

void test_rk23_convergence_sanity_check()
{
    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {1.0};
    quadrature.weights = {1.0};

    qmomkit::NodeAdvectionOptions loose_options;
    loose_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    loose_options.adaptive_options.relative_tolerance = 1.0e-3;
    loose_options.adaptive_options.absolute_tolerance = 1.0e-6;
    loose_options.adaptive_options.initial_step = 0.1;

    qmomkit::NodeAdvectionOptions tight_options;
    tight_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    tight_options.adaptive_options.relative_tolerance = 1.0e-8;
    tight_options.adaptive_options.absolute_tolerance = 1.0e-11;
    tight_options.adaptive_options.initial_step = 0.1;

    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>{xi[0]}; });
    const qmomkit::QmomNodeAdvectionResult loose = qmomkit::advance_nodes(quadrature, 0.0, 1.0, velocity, loose_options);
    const qmomkit::QmomNodeAdvectionResult tight = qmomkit::advance_nodes(quadrature, 0.0, 1.0, velocity, tight_options);

    QMOMKIT_VERIFY_CHECK(loose.ok());
    QMOMKIT_VERIFY_CHECK(tight.ok());
    const double exact = std::exp(1.0);
    const double loose_error = std::abs(loose.quadrature.abscissas[0] - exact);
    const double tight_error = std::abs(tight.quadrature.abscissas[0] - exact);

    QMOMKIT_VERIFY_CHECK(tight.diagnostics.step_count >= loose.diagnostics.step_count);
    QMOMKIT_VERIFY_CHECK(tight_error < loose_error);
    QMOMKIT_VERIFY_CHECK(tight_error < 1.0e-7);
}

void test_ssp_rk3_linear_phase_space_velocity_for_qmom()
{
    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {-0.5, 1.25};
    quadrature.weights = {0.4, 0.6};

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    constexpr double rate = -0.75;
    constexpr double dt = 0.05;
    const qmomkit::QmomNodeAdvectionResult result = qmomkit::advance_nodes(
        quadrature, 0.0, dt,
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>{rate * xi[0]}; }), options);

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.step_count == 1);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.rejected_step_count == 0);
    assert_near(result.quadrature.abscissas[0], quadrature.abscissas[0] * std::exp(rate * dt), 2.0e-7);
    assert_near(result.quadrature.abscissas[1], quadrature.abscissas[1] * std::exp(rate * dt), 2.0e-7);
    QMOMKIT_VERIFY_CHECK(result.quadrature.weights == quadrature.weights);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"constant_phase_space_velocity_for_qmom", test_constant_phase_space_velocity_for_qmom},
            {"linear_phase_space_velocity_for_qmom", test_linear_phase_space_velocity_for_qmom},
            {"multi_dimensional_phase_space_evolution", test_multi_dimensional_phase_space_evolution},
            {"ssp_rk3_multi_dimensional_phase_space_evolution", test_ssp_rk3_multi_dimensional_phase_space_evolution},
            {"rk23_convergence_sanity_check", test_rk23_convergence_sanity_check},
            {"ssp_rk3_linear_phase_space_velocity_for_qmom", test_ssp_rk3_linear_phase_space_velocity_for_qmom},
        });
}
