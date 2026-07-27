#include "qmomkit/qmomkit.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
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

} // namespace

int main()
{
    test_adaptive_quadrature_structure_is_preserved();
    test_default_scheme_is_ssp_rk3();
    test_rk23_failure_diagnostics_are_reported();
    return 0;
}
