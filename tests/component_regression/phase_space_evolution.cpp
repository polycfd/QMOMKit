#include "../verification/verification_helpers.hpp"
#include "internal/cqmom_case.hpp"
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

void assert_vector_near(const std::vector<double>& actual, const std::vector<double>& expected, const double tolerance)
{
    QMOMKIT_VERIFY_CHECK(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        assert_near(actual[i], expected[i], tolerance);
    }
}

std::vector<double> raw_moments_1d(const std::vector<double>& abscissas, const std::vector<double>& weights, const std::size_t moment_count)
{
    std::vector<double> moments(moment_count, 0.0);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        for (std::size_t node = 0; node < weights.size(); ++node)
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

qmomkit::RawMomentTensor fixed_cqmom_moments()
{
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    return {extents, raw_moments_from_terminals(extents, weights, abscissas)};
}

qmomkit::NodeAdvectionOptions tight_node_options()
{
    qmomkit::NodeAdvectionOptions options;
    options.adaptive_options.relative_tolerance = 1.0e-9;
    options.adaptive_options.absolute_tolerance = 1.0e-12;
    options.adaptive_options.initial_step = 0.01;
    return options;
}

qmomkit::QmomState qmom_state_from_moments(const std::vector<double>& moments, const std::size_t node_count, const qmomkit::QmomOptions& options = {})
{
    return qmomkit::QmomCase{qmomkit::MomentSet(moments), node_count, options}.initial_state();
}

void test_zero_advection_and_zero_source_preserve_qmom()
{
    const std::vector<double> moments = raw_moments_1d({-1.0, 2.0}, {0.25, 0.75}, 4);
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 2);
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_phase_space(
        state, 0.0, 0.4,
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); }),
        [](double, const qmomkit::QmomSourceContext& context) { return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)}; }, options);

    QMOMKIT_VERIFY_CHECK(result.ok());
    assert_vector_near(result.state.moments(), moments, 1.0e-12);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::none);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.first_reconstruction_diagnostics.success);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.final_reconstruction_diagnostics.success);
}

void test_combined_advection_and_source_qmom_analytical()
{
    const std::vector<double> moments{2.0, 3.0};
    constexpr double shift_rate = 1.25;
    constexpr double source_rate = 0.4;
    constexpr double dt = 0.2;
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>& xi)
        {
            QMOMKIT_VERIFY_CHECK(xi.size() == 1);
            return std::vector<double>{shift_rate};
        });
    const auto source = [](double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> rhs(context.size(), 0.0);
        for (std::size_t i = 0; i < context.size(); ++i)
        {
            rhs[i] = source_rate * context.value(i);
        }
        return qmomkit::QmomMomentRate{std::move(rhs)};
    };

    const qmomkit::QmomState state = qmom_state_from_moments(moments, 1);
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_phase_space(state, 0.0, dt, velocity, source, options);

    QMOMKIT_VERIFY_CHECK(result.ok());
    const double initial_weight = moments[0];
    const double initial_x = moments[1] / moments[0];
    const double source_factor =
        1.0 + source_rate * dt + 0.5 * source_rate * source_rate * dt * dt + source_rate * source_rate * source_rate * dt * dt * dt / 6.0;
    const double final_weight = source_factor * initial_weight;
    const double final_x = initial_x + shift_rate * dt;
    assert_near(result.state.moments()[0], final_weight, 1.0e-13);
    assert_near(result.state.moments()[1], final_weight * final_x, 1.0e-13);
}

void test_cqmom_case_preserves_phase_space_ordering()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& state)
    { return qmomkit::CqmomMomentRate{std::vector<double>(state.size(), 0.0)}; };
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>& xi)
        {
            QMOMKIT_VERIFY_CHECK(xi.size() == 2);
            return std::vector<double>{0.25, -0.5};
        });

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.node_options = tight_node_options();
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomCase cqmom_case{qmomkit::extract_cqmom_moments(moments, layout), {}, {}};
    const qmomkit::CqmomState state = qmomkit::internal::CqmomConstructionAccess::create_state(cqmom_case);
    const qmomkit::CqmomPhaseSpaceEvolutionResult result = qmomkit::advance_cqmom_phase_space(state, 0.0, 0.2, velocity, zero_source, options);
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.quadrature.has_value());
    QMOMKIT_VERIFY_CHECK(result.quadrature->has_rectangular_storage());

    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_case);
    QMOMKIT_VERIFY_CHECK(initial.ok());
    const qmomkit::CqmomNodeAdvectionResult direct = qmomkit::advance_nodes(initial.quadrature, 0.0, 0.2, velocity, options.node_options);
    QMOMKIT_VERIFY_CHECK(direct.ok());
    // The state carries the moments its hierarchy reconstructs from, so the
    // independently advected reference is reduced the same way before comparison.
    const qmomkit::CqmomMomentData expected = qmomkit::reconstruct_cqmom_moments(direct.quadrature, layout);
    assert_vector_near(result.state.moment_values(), expected.values, 5.0e-9);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"zero_advection_and_zero_source_preserve_qmom", test_zero_advection_and_zero_source_preserve_qmom},
            {"combined_advection_and_source_qmom_analytical", test_combined_advection_and_source_qmom_analytical},
            {"cqmom_case_preserves_phase_space_ordering", test_cqmom_case_preserves_phase_space_ordering},
        });
}
