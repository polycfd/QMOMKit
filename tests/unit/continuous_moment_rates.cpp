#include "internal/continuous_moment_rates.hpp"

#include "internal/cqmom_case.hpp"
#include "internal/moment_reconstruction.hpp"
#include "internal/performance.hpp"
#include "qmomkit/qmomkit.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

void assert_vector_near(const std::vector<double>& actual, const std::vector<double>& expected, const double tolerance)
{
    assert(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        assert_near(actual[i], expected[i], tolerance);
    }
}

bool same_layout(const qmomkit::CqmomMomentLayout& lhs, const qmomkit::CqmomMomentLayout& rhs)
{
    return lhs.hierarchy.node_counts() == rhs.hierarchy.node_counts() && lhs.hierarchy.conditioning_order() == rhs.hierarchy.conditioning_order() &&
           lhs.required_indices == rhs.required_indices;
}

double integer_power(const double base, const std::size_t exponent)
{
    double value = 1.0;
    for (std::size_t i = 0; i < exponent; ++i)
    {
        value *= base;
    }
    return value;
}

std::vector<double> qmom_moments_from_terminals(const std::vector<double>& weights, const std::vector<double>& abscissas,
                                                const std::size_t moment_count)
{
    std::vector<double> moments(moment_count, 0.0);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            moments[order] += weights[node] * integer_power(abscissas[node], order);
        }
    }
    return moments;
}

std::vector<qmomkit::MultiIndex> qmom_indices(const std::size_t moment_count)
{
    std::vector<qmomkit::MultiIndex> indices;
    indices.reserve(moment_count);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        indices.push_back(qmomkit::MultiIndex{order});
    }
    return indices;
}

std::vector<double> selected_terminal_moments(const std::vector<double>& weights, const std::vector<std::vector<double>>& abscissas,
                                              const std::vector<qmomkit::MultiIndex>& indices)
{
    std::vector<double> moments;
    moments.reserve(indices.size());
    for (const qmomkit::MultiIndex& index : indices)
    {
        double value = 0.0;
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            double contribution = weights[node];
            for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
            {
                contribution *= integer_power(abscissas[dimension][node], index[dimension]);
            }
            value += contribution;
        }
        moments.push_back(value);
    }
    return moments;
}

qmomkit::CqmomMomentData cqmom_moment_data_from_terminals(const qmomkit::CqmomMomentLayout& layout, const std::vector<double>& weights,
                                                          const std::vector<std::vector<double>>& abscissas)
{
    return qmomkit::CqmomMomentData{layout, selected_terminal_moments(weights, abscissas, layout.required_indices)};
}

std::vector<std::vector<double>> manual_terminal_velocities(const qmomkit::PhaseSpaceVelocity& velocity, const double time,
                                                            const std::vector<std::vector<double>>& abscissas)
{
    const std::size_t dimensions = abscissas.size();
    const std::size_t node_count = abscissas.front().size();
    std::vector<std::vector<double>> velocities(node_count, std::vector<double>(dimensions, 0.0));
    std::vector<double> coordinate(dimensions, 0.0);
    for (std::size_t node = 0; node < node_count; ++node)
    {
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
        {
            coordinate[dimension] = abscissas[dimension][node];
        }
        velocities[node] = velocity(time, coordinate);
    }
    return velocities;
}

std::vector<double> independently_composed_qmom_rates(const std::vector<double>& moments, const std::size_t node_count, const double time,
                                                      const qmomkit::PhaseSpaceVelocity& velocity, const qmomkit::QmomOptions& options = {})
{
    const qmomkit::QmomResult quadrature = qmomkit::compute_qmom(qmomkit::MomentSet(moments), node_count, options);
    assert(quadrature.ok());
    const std::vector<std::vector<double>> abscissas{quadrature.abscissas};
    return qmomkit::internal::reconstruct_continuous_moment_rates(
        quadrature.weights, abscissas, manual_terminal_velocities(velocity, time, abscissas), qmom_indices(moments.size()));
}

std::vector<double> independently_composed_cqmom_rates(const qmomkit::CqmomLeanCase& lean_case, const double time,
                                                       const qmomkit::PhaseSpaceVelocity& velocity)
{
    const qmomkit::CqmomResult quadrature = qmomkit::internal::compute_trusted_cqmom(lean_case);
    assert(quadrature.ok());
    return qmomkit::internal::reconstruct_continuous_moment_rates(quadrature.quadrature.weights(), quadrature.quadrature.abscissas(),
                                                                  manual_terminal_velocities(velocity, time, quadrature.quadrature.abscissas()),
                                                                  lean_case.moments.layout.required_indices);
}

void test_qmom_one_node_constant_velocity()
{
    const std::vector<double> moments = qmom_moments_from_terminals({2.0}, {3.0}, 4);
    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>&) { return std::vector<double>{0.5}; };

    const qmomkit::internal::QmomContinuousMomentRateResult result =
        qmomkit::internal::evaluate_qmom_continuous_moment_rates(moments, 1, 0.25, velocity);

    assert(result.ok());
    assert(result.inversion_status == qmomkit::QmomStatus::success);
    assert(result.rates[0] == 0.0);
    assert_vector_near(result.rates, {0.0, 1.0, 6.0, 27.0}, 1.0e-13);
}

void test_qmom_two_node_distinct_velocities_and_independent_composition()
{
    const std::vector<double> weights{0.25, 0.75};
    const std::vector<double> abscissas{-1.0, 2.0};
    const std::vector<double> moments = qmom_moments_from_terminals(weights, abscissas, 4);
    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>& coordinate)
    { return std::vector<double>{coordinate[0] < 0.0 ? 0.5 : -0.25}; };

    const qmomkit::internal::QmomContinuousMomentRateResult result =
        qmomkit::internal::evaluate_qmom_continuous_moment_rates(moments, 2, 0.5, velocity);

    assert(result.ok());
    assert(result.rates[0] == 0.0);
    assert_vector_near(result.rates, {0.0, -0.0625, -1.0, -1.875}, 1.0e-12);
    assert_vector_near(result.rates, independently_composed_qmom_rates(moments, 2, 0.5, velocity), 1.0e-12);
}

void test_qmom_zero_and_linear_velocity()
{
    {
        const std::vector<double> moments = qmom_moments_from_terminals({1.0, 2.0, 3.0}, {-2.0, 0.0, 4.0}, 6);
        const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>&) { return std::vector<double>{0.0}; };

        const qmomkit::internal::QmomContinuousMomentRateResult result =
            qmomkit::internal::evaluate_qmom_continuous_moment_rates(moments, 3, 0.0, velocity);

        assert(result.ok());
        assert_vector_near(result.rates, std::vector<double>(moments.size(), 0.0), 0.0);
    }
    {
        const std::vector<double> moments = qmom_moments_from_terminals({0.5, 1.5}, {1.0, 3.0}, 4);
        const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>& coordinate)
        { return std::vector<double>{0.4 - 0.2 * coordinate[0]}; };

        const qmomkit::internal::QmomContinuousMomentRateResult result =
            qmomkit::internal::evaluate_qmom_continuous_moment_rates(moments, 2, 0.0, velocity);

        assert(result.ok());
        assert(result.rates[0] == 0.0);
        assert_vector_near(result.rates, {0.0, -0.2, -1.6, -7.8}, 1.0e-12);
    }
}

void test_qmom_failure_paths()
{
    const qmomkit::internal::QmomContinuousMomentRateResult inversion_failure = qmomkit::internal::evaluate_qmom_continuous_moment_rates(
        {1.0, 0.0, -1.0, 0.0}, 2, 0.0, [](double, const std::vector<double>&) { return std::vector<double>{0.0}; });
    assert(!inversion_failure.ok());
    assert(inversion_failure.status == qmomkit::internal::ContinuousMomentRateStatus::inversion_failed);
    assert(inversion_failure.inversion_status == qmomkit::QmomStatus::not_realizable);
    assert(inversion_failure.rates.empty());

    const qmomkit::internal::QmomContinuousMomentRateResult missing_velocity =
        qmomkit::internal::evaluate_qmom_continuous_moment_rates({1.0, 0.0}, 1, 0.0, {});
    assert(!missing_velocity.ok());
    assert(missing_velocity.status == qmomkit::internal::ContinuousMomentRateStatus::invalid_input);

    const std::vector<double> moments = qmom_moments_from_terminals({1.0}, {2.0}, 3);
    const qmomkit::internal::QmomContinuousMomentRateResult wrong_size = qmomkit::internal::evaluate_qmom_continuous_moment_rates(
        moments, 1, 0.0, [](double, const std::vector<double>&) { return std::vector<double>{1.0, 2.0}; });
    assert(!wrong_size.ok());
    assert(wrong_size.status == qmomkit::internal::ContinuousMomentRateStatus::velocity_failed);
    assert(wrong_size.message == "Phase-space velocity callback returned the wrong dimension.");

    const qmomkit::internal::QmomContinuousMomentRateResult non_finite = qmomkit::internal::evaluate_qmom_continuous_moment_rates(
        moments, 1, 0.0, [](double, const std::vector<double>&) { return std::vector<double>{std::numeric_limits<double>::infinity()}; });
    assert(!non_finite.ok());
    assert(non_finite.status == qmomkit::internal::ContinuousMomentRateStatus::velocity_failed);
    assert(non_finite.message == "Phase-space velocity callback returned a non-finite value.");

    const qmomkit::internal::QmomContinuousMomentRateResult thrown = qmomkit::internal::evaluate_qmom_continuous_moment_rates(
        moments, 1, 0.0, [](double, const std::vector<double>&) -> std::vector<double> { throw std::runtime_error("velocity model failed"); });
    assert(!thrown.ok());
    assert(thrown.status == qmomkit::internal::ContinuousMomentRateStatus::velocity_failed);
    assert(thrown.message.find("velocity model failed") != std::string::npos);
}

void test_cqmom_identity_ordering_and_independent_composition()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    const qmomkit::CqmomLeanCase lean_case{cqmom_moment_data_from_terminals(layout, weights, abscissas)};
    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>& coordinate)
    { return std::vector<double>{0.25 * coordinate[0], -0.5 * coordinate[1]}; };

    const qmomkit::internal::CqmomContinuousMomentRateResult result =
        qmomkit::internal::evaluate_cqmom_continuous_moment_rates(lean_case, 0.0, velocity);

    assert(result.ok());
    assert(result.rates.has_value());
    assert(same_layout(result.rates->layout, layout));
    assert(result.rates->values.size() == layout.required_indices.size());
    assert_vector_near(result.rates->values, independently_composed_cqmom_rates(lean_case, 0.0, velocity), 1.0e-11);
    assert(result.rates->values.front() == 0.0);
}

void test_cqmom_non_identity_ordering_and_global_terminal_weights()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    const qmomkit::CqmomLeanCase lean_case{cqmom_moment_data_from_terminals(layout, weights, abscissas)};
    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>& coordinate)
    { return std::vector<double>{0.2 + coordinate[1], -0.1 * coordinate[0]}; };

    const qmomkit::internal::CqmomContinuousMomentRateResult result =
        qmomkit::internal::evaluate_cqmom_continuous_moment_rates(lean_case, 0.5, velocity);
    const qmomkit::CqmomResult quadrature = qmomkit::internal::compute_trusted_cqmom(lean_case);
    assert(quadrature.ok());
    const std::vector<double> expected = qmomkit::internal::reconstruct_continuous_moment_rates(
        quadrature.quadrature.weights(), quadrature.quadrature.abscissas(),
        manual_terminal_velocities(velocity, 0.5, quadrature.quadrature.abscissas()), layout.required_indices);

    assert(result.ok());
    assert(result.rates.has_value());
    assert(same_layout(result.rates->layout, layout));
    assert_vector_near(result.rates->values, expected, 1.0e-11);
}

void test_cqmom_zero_velocity_and_adaptive_terminal_arrays()
{
    {
        const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
        const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
        const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
        const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
        const qmomkit::CqmomLeanCase lean_case{cqmom_moment_data_from_terminals(layout, weights, abscissas)};
        const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>& coordinate)
        { return std::vector<double>(coordinate.size(), 0.0); };

        const qmomkit::internal::CqmomContinuousMomentRateResult result =
            qmomkit::internal::evaluate_cqmom_continuous_moment_rates(lean_case, 0.0, velocity);

        assert(result.ok());
        assert(result.rates.has_value());
        assert_vector_near(result.rates->values, std::vector<double>(layout.required_indices.size(), 0.0), 0.0);
    }
    {
        const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
        const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
        const std::vector<double> weights{0.1, 0.3, 0.6};
        const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 1.0}, {-2.0, 2.0, 3.0}};
        qmomkit::CqmomLeanCase lean_case{cqmom_moment_data_from_terminals(layout, weights, abscissas)};
        lean_case.qmom_options = {true, 1.0e-10, 1.0e-10, 1.0e-12};
        lean_case.adaptive_branch_counts = true;
        const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>& coordinate)
        { return std::vector<double>{0.25 * coordinate[0], -0.5 * coordinate[1]}; };

        const qmomkit::CqmomResult quadrature = qmomkit::internal::compute_trusted_cqmom(lean_case);
        assert(quadrature.ok());
        assert(!quadrature.quadrature.has_rectangular_storage());
        assert(quadrature.quadrature.weights().size() == 3);

        const qmomkit::internal::CqmomContinuousMomentRateResult result =
            qmomkit::internal::evaluate_cqmom_continuous_moment_rates(lean_case, 0.0, velocity);

        assert(result.ok());
        assert(result.rates.has_value());
        assert_vector_near(result.rates->values, independently_composed_cqmom_rates(lean_case, 0.0, velocity), 1.0e-11);
        assert(result.rates->values.front() == 0.0);
    }
}

void test_cqmom_failure_paths()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    std::vector<double> invalid_values(layout.required_indices.size(), 0.0);
    invalid_values.front() = -1.0;
    const qmomkit::CqmomLeanCase invalid_case{qmomkit::CqmomMomentData{layout, invalid_values}};

    const qmomkit::internal::CqmomContinuousMomentRateResult inversion_failure = qmomkit::internal::evaluate_cqmom_continuous_moment_rates(
        invalid_case, 0.0, [](double, const std::vector<double>& coordinate) { return std::vector<double>(coordinate.size(), 0.0); });
    assert(!inversion_failure.ok());
    assert(inversion_failure.status == qmomkit::internal::ContinuousMomentRateStatus::inversion_failed);
    assert(inversion_failure.inversion_status == qmomkit::CqmomStatus::qmom_failed);
    assert(!inversion_failure.rates.has_value());
    assert(!inversion_failure.diagnostics.empty());

    const qmomkit::CqmomLeanCase valid_case{
        cqmom_moment_data_from_terminals(layout, {0.2, 0.3, 0.1, 0.4}, {{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}})};

    const qmomkit::internal::CqmomContinuousMomentRateResult wrong_size = qmomkit::internal::evaluate_cqmom_continuous_moment_rates(
        valid_case, 0.0, [](double, const std::vector<double>&) { return std::vector<double>{1.0}; });
    assert(!wrong_size.ok());
    assert(wrong_size.status == qmomkit::internal::ContinuousMomentRateStatus::velocity_failed);
    assert(wrong_size.message == "Phase-space velocity callback returned the wrong dimension.");

    const qmomkit::internal::CqmomContinuousMomentRateResult non_finite =
        qmomkit::internal::evaluate_cqmom_continuous_moment_rates(valid_case, 0.0,
                                                                  [](double, const std::vector<double>& coordinate)
                                                                  {
                                                                      std::vector<double> value(coordinate.size(), 0.0);
                                                                      value.front() = std::numeric_limits<double>::quiet_NaN();
                                                                      return value;
                                                                  });
    assert(!non_finite.ok());
    assert(non_finite.status == qmomkit::internal::ContinuousMomentRateStatus::velocity_failed);
    assert(non_finite.message == "Phase-space velocity callback returned a non-finite value.");

    const qmomkit::internal::CqmomContinuousMomentRateResult thrown = qmomkit::internal::evaluate_cqmom_continuous_moment_rates(
        valid_case, 0.0, [](double, const std::vector<double>&) -> std::vector<double> { throw std::runtime_error("CQMOM velocity model failed"); });
    assert(!thrown.ok());
    assert(thrown.status == qmomkit::internal::ContinuousMomentRateStatus::velocity_failed);
    assert(thrown.message.find("CQMOM velocity model failed") != std::string::npos);
}

void test_new_evaluators_do_not_record_native_node_rhs_counters()
{
    const std::vector<double> qmom_moments = qmom_moments_from_terminals({0.25, 0.75}, {-1.0, 2.0}, 4);
    const qmomkit::PhaseSpaceVelocity qmom_velocity = [](double, const std::vector<double>& coordinate)
    { return std::vector<double>{0.1 * coordinate[0]}; };

    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomLeanCase cqmom_case{
        cqmom_moment_data_from_terminals(layout, {0.2, 0.3, 0.1, 0.4}, {{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}})};
    const qmomkit::PhaseSpaceVelocity cqmom_velocity = [](double, const std::vector<double>& coordinate)
    { return std::vector<double>{coordinate[0], coordinate[1]}; };

    qmomkit::internal::performance::CounterSession session;
    {
        qmomkit::internal::performance::CounterScope scope(session);
        const qmomkit::internal::QmomContinuousMomentRateResult qmom_result =
            qmomkit::internal::evaluate_qmom_continuous_moment_rates(qmom_moments, 2, 0.0, qmom_velocity);
        const qmomkit::internal::CqmomContinuousMomentRateResult cqmom_result =
            qmomkit::internal::evaluate_cqmom_continuous_moment_rates(cqmom_case, 0.0, cqmom_velocity);
        assert(qmom_result.ok());
        assert(cqmom_result.ok());
    }

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
    assert(counters.time_integration.node_rhs_evaluations == 0);
}

} // namespace

int main()
{
    test_qmom_one_node_constant_velocity();
    test_qmom_two_node_distinct_velocities_and_independent_composition();
    test_qmom_zero_and_linear_velocity();
    test_qmom_failure_paths();
    test_cqmom_identity_ordering_and_independent_composition();
    test_cqmom_non_identity_ordering_and_global_terminal_weights();
    test_cqmom_zero_velocity_and_adaptive_terminal_arrays();
    test_cqmom_failure_paths();
    test_new_evaluators_do_not_record_native_node_rhs_counters();
    return 0;
}
