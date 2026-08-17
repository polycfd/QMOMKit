#include "qmomkit/qmomkit.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
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

void assert_contains(const std::string& text, const std::string& needle)
{
    if (text.find(needle) == std::string::npos)
    {
        std::cerr << "Missing expected text:\n" << needle << "\nActual message:\n" << text << '\n';
        assert(false);
    }
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

void test_public_qmom_advection_rates_and_status_mapping()
{
    const std::vector<double> moments = qmom_moments_from_terminals({2.0}, {3.0}, 4);
    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>&) { return std::vector<double>{0.5}; };

    const qmomkit::QmomAdvectionMomentRateResult success = qmomkit::evaluate_qmom_advection_moment_rates(moments, 1, 0.25, velocity);
    assert(success.ok());
    assert(success.status == qmomkit::AdvectionMomentRateStatus::success);
    assert(success.inversion_status == qmomkit::QmomStatus::success);
    assert(success.rates[0] == 0.0);
    assert_vector_near(success.rates, {0.0, 1.0, 6.0, 27.0}, 1.0e-13);

    const qmomkit::QmomAdvectionMomentRateResult inversion_failed = qmomkit::evaluate_qmom_advection_moment_rates(
        {1.0, 0.0, -1.0, 0.0}, 2, 0.0, [](double, const std::vector<double>&) { return std::vector<double>{0.0}; });
    assert(!inversion_failed.ok());
    assert(inversion_failed.status == qmomkit::AdvectionMomentRateStatus::inversion_failed);
    assert(inversion_failed.inversion_status == qmomkit::QmomStatus::not_realizable);
    assert(inversion_failed.rates.empty());

    const qmomkit::QmomAdvectionMomentRateResult velocity_failed = qmomkit::evaluate_qmom_advection_moment_rates(
        moments, 1, 0.0, [](double, const std::vector<double>&) { return std::vector<double>{1.0, 2.0}; });
    assert(!velocity_failed.ok());
    assert(velocity_failed.status == qmomkit::AdvectionMomentRateStatus::velocity_failed);
    assert(velocity_failed.inversion_status == qmomkit::QmomStatus::success);
    assert_contains(velocity_failed.message, "Phase-space velocity callback returned the wrong dimension:");
    assert_contains(velocity_failed.message, "expected size: 1");
    assert_contains(velocity_failed.message, "returned size: 2");

    const qmomkit::QmomAdvectionMomentRateResult invalid_input = qmomkit::evaluate_qmom_advection_moment_rates({}, 1, 0.0, velocity);
    assert(!invalid_input.ok());
    assert(invalid_input.status == qmomkit::AdvectionMomentRateStatus::invalid_input);
    assert(invalid_input.inversion_status == qmomkit::QmomStatus::invalid_input);
    assert(invalid_input.rates.empty());
}

void test_public_cqmom_advection_rates_preserve_layout_and_status()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomLeanCase lean_case{
        cqmom_moment_data_from_terminals(layout, {0.2, 0.3, 0.1, 0.4}, {{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}})};
    const qmomkit::PhaseSpaceVelocity zero_velocity = [](double, const std::vector<double>& coordinate)
    { return std::vector<double>(coordinate.size(), 0.0); };

    const qmomkit::CqmomAdvectionMomentRateResult success = qmomkit::evaluate_cqmom_advection_moment_rates(lean_case, 0.5, zero_velocity);
    assert(success.ok());
    assert(success.status == qmomkit::AdvectionMomentRateStatus::success);
    assert(success.inversion_status == qmomkit::CqmomStatus::success);
    assert(success.rates.has_value());
    assert(same_layout(success.rates->layout, layout));
    assert_vector_near(success.rates->values, std::vector<double>(layout.required_indices.size(), 0.0), 0.0);
    assert(!success.diagnostics.empty());

    std::vector<double> invalid_values(layout.required_indices.size(), 0.0);
    invalid_values.front() = -1.0;
    const qmomkit::CqmomLeanCase invalid_case{qmomkit::CqmomMomentData{layout, invalid_values}};
    const qmomkit::CqmomAdvectionMomentRateResult inversion_failed = qmomkit::evaluate_cqmom_advection_moment_rates(invalid_case, 0.0, zero_velocity);
    assert(!inversion_failed.ok());
    assert(inversion_failed.status == qmomkit::AdvectionMomentRateStatus::inversion_failed);
    assert(inversion_failed.inversion_status == qmomkit::CqmomStatus::qmom_failed);
    assert(!inversion_failed.rates.has_value());
    assert(!inversion_failed.diagnostics.empty());

    const qmomkit::CqmomAdvectionMomentRateResult velocity_failed = qmomkit::evaluate_cqmom_advection_moment_rates(
        lean_case, 0.0, [](double, const std::vector<double>&) { return std::vector<double>{std::numeric_limits<double>::quiet_NaN(), 0.0}; });
    assert(!velocity_failed.ok());
    assert(velocity_failed.status == qmomkit::AdvectionMomentRateStatus::velocity_failed);
    assert(velocity_failed.inversion_status == qmomkit::CqmomStatus::success);
    assert_contains(velocity_failed.message, "Phase-space velocity callback returned a non-finite value:");
    assert_contains(velocity_failed.message, "time: 0");
    assert_contains(velocity_failed.message, "dimension: 0");
    assert_contains(velocity_failed.message, "value: nan");
    assert_contains(velocity_failed.message, "state: [-0.99999999999999967, -0.49999999999999961]");
    assert_contains(velocity_failed.message, "velocity: [nan, 0]");

    qmomkit::CqmomLeanCase malformed_case = lean_case;
    malformed_case.moments.values.pop_back();
    const qmomkit::CqmomAdvectionMomentRateResult invalid_input = qmomkit::evaluate_cqmom_advection_moment_rates(malformed_case, 0.0, zero_velocity);
    assert(!invalid_input.ok());
    assert(invalid_input.status == qmomkit::AdvectionMomentRateStatus::invalid_input);
    assert(invalid_input.inversion_status == qmomkit::CqmomStatus::invalid_input);
    assert(!invalid_input.rates.has_value());
}

void test_public_cqmom_advection_rates_preserve_adaptive_layout()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::CqmomLeanCase lean_case{cqmom_moment_data_from_terminals(layout, {0.1, 0.3, 0.6}, {{-1.0, -1.0, 1.0}, {-2.0, 2.0, 3.0}})};
    lean_case.qmom_options = {true, 1.0e-10, 1.0e-10, 1.0e-12};
    lean_case.adaptive_branch_counts = true;

    const qmomkit::CqmomAdvectionMomentRateResult result = qmomkit::evaluate_cqmom_advection_moment_rates(
        lean_case, 0.0, [](double, const std::vector<double>& coordinate) { return std::vector<double>(coordinate.size(), 0.0); });
    assert(result.ok());
    assert(result.rates.has_value());
    assert(same_layout(result.rates->layout, layout));
    assert_vector_near(result.rates->values, std::vector<double>(layout.required_indices.size(), 0.0), 0.0);
}

void test_public_qmom_moment_source_rates_and_status_mapping()
{
    const std::vector<double> moments{1.0, 2.0, 3.0};
    std::size_t calls = 0;
    const qmomkit::QmomMomentSourceRateResult success =
        qmomkit::evaluate_qmom_moment_source_rates(0.25, moments,
                                                   [&calls](const double time, const qmomkit::QmomSourceContext& state)
                                                   {
                                                       ++calls;
                                                       assert_near(time, 0.25, 0.0);
                                                       assert_vector_near(state.values(), {1.0, 2.0, 3.0}, 0.0);
                                                       return std::vector<double>{state.value(2), state.value(0), state.value(1)};
                                                   });
    assert(success.ok());
    assert(success.status == qmomkit::MomentSourceRateStatus::success);
    assert_vector_near(success.rates, {3.0, 1.0, 2.0}, 0.0);
    assert(calls == 1);

    const qmomkit::QmomMomentSourceRateResult invalid_rate = qmomkit::evaluate_qmom_moment_source_rates(
        0.0, moments, [](double, const qmomkit::QmomSourceContext& state) { return std::vector<double>(state.size() + 1, 0.0); });
    assert(!invalid_rate.ok());
    assert(invalid_rate.status == qmomkit::MomentSourceRateStatus::invalid_rate);
    assert(invalid_rate.message == "Moment-source callback rate values must match the moment state size.");

    const qmomkit::QmomMomentSourceRateResult non_finite_rate =
        qmomkit::evaluate_qmom_moment_source_rates(0.0, moments,
                                                   [](double, const qmomkit::QmomSourceContext& state)
                                                   {
                                                       std::vector<double> rates(state.size(), 0.0);
                                                       rates.back() = std::numeric_limits<double>::infinity();
                                                       return rates;
                                                   });
    assert(!non_finite_rate.ok());
    assert(non_finite_rate.status == qmomkit::MomentSourceRateStatus::invalid_rate);
    assert(non_finite_rate.message == "Moment-source callback rate values must be finite.");

    const qmomkit::QmomMomentSourceRateResult invalid_input = qmomkit::evaluate_qmom_moment_source_rates(0.0, {}, {});
    assert(!invalid_input.ok());
    assert(invalid_input.status == qmomkit::MomentSourceRateStatus::invalid_input);
    assert(invalid_input.rates.empty());

    const qmomkit::QmomMomentSourceRateResult source_failed = qmomkit::evaluate_qmom_moment_source_rates(
        0.0, moments, [](double, const qmomkit::QmomSourceContext&) -> std::vector<double> { throw std::runtime_error("source failed"); });
    assert(!source_failed.ok());
    assert(source_failed.status == qmomkit::MomentSourceRateStatus::source_failed);
    assert(source_failed.message.find("source failed") != std::string::npos);
}

void test_public_qmom_moment_source_rates_with_boundary()
{
    const std::vector<double> moments{1.0, 2.0, 3.0};
    const std::vector<double> original = moments;
    std::size_t source_calls = 0;
    std::size_t boundary_calls = 0;
    int callback_order = 0;

    const qmomkit::QmomMomentSourceRateResult result = qmomkit::evaluate_qmom_moment_source_rates(
        0.25, moments,
        [&moments, &source_calls, &callback_order](const double time, const qmomkit::QmomSourceContext& state)
        {
            ++source_calls;
            assert(callback_order == 0);
            callback_order = 1;
            assert_near(time, 0.25, 0.0);
            assert(&state.values() == &moments);
            return std::vector<double>{1.0, 2.0, 3.0};
        },
        [&moments, &boundary_calls, &callback_order](const double time, const std::vector<double>& state)
        {
            ++boundary_calls;
            assert(callback_order == 1);
            callback_order = 2;
            assert_near(time, 0.25, 0.0);
            assert(&state == &moments);
            return qmomkit::MomentBoundaryRate{std::vector<double>{0.5, -1.0, 2.0}};
        });

    assert(result.ok());
    assert(result.status == qmomkit::MomentSourceRateStatus::success);
    assert_vector_near(result.rates, {1.5, 1.0, 5.0}, 0.0);
    assert(source_calls == 1);
    assert(boundary_calls == 1);
    assert(callback_order == 2);
    assert(moments == original);
}

void test_public_qmom_moment_source_rates_with_zero_boundary()
{
    const std::vector<double> moments{1.0, 2.0};
    const qmomkit::QmomMomentSourceRateResult result = qmomkit::evaluate_qmom_moment_source_rates(
        0.0, moments, [](double, const qmomkit::QmomSourceContext&) { return std::vector<double>{3.0, -4.0}; },
        [](double, const std::vector<double>& state) { return qmomkit::MomentBoundaryRate{std::vector<double>(state.size(), 0.0)}; });

    assert(result.ok());
    assert_vector_near(result.rates, {3.0, -4.0}, 0.0);
}

void test_public_qmom_moment_source_rates_validate_inputs_before_callbacks()
{
    bool source_called = false;
    bool boundary_called = false;
    const qmomkit::MomentSource source = [&source_called](double, const qmomkit::QmomSourceContext&)
    {
        source_called = true;
        return std::vector<double>{0.0};
    };
    const qmomkit::MomentBoundaryContribution boundary = [&boundary_called](double, const std::vector<double>&)
    {
        boundary_called = true;
        return qmomkit::MomentBoundaryRate{std::vector<double>{0.0}};
    };

    const qmomkit::QmomMomentSourceRateResult invalid_time =
        qmomkit::evaluate_qmom_moment_source_rates(std::numeric_limits<double>::quiet_NaN(), {1.0}, source, boundary);
    assert(!invalid_time.ok());
    assert(invalid_time.status == qmomkit::MomentSourceRateStatus::invalid_input);
    assert(invalid_time.message == "Moment-source rate evaluation time must be finite.");
    assert(!source_called);
    assert(!boundary_called);

    const qmomkit::QmomMomentSourceRateResult invalid_moments = qmomkit::evaluate_qmom_moment_source_rates(0.0, {}, source, boundary);
    assert(!invalid_moments.ok());
    assert(invalid_moments.status == qmomkit::MomentSourceRateStatus::invalid_input);
    assert(invalid_moments.message == "Moment-source rate evaluation requires at least one moment.");
    assert(!source_called);
    assert(!boundary_called);

    const qmomkit::QmomMomentSourceRateResult missing_boundary = qmomkit::evaluate_qmom_moment_source_rates(0.0, {1.0}, source, {});
    assert(!missing_boundary.ok());
    assert(missing_boundary.status == qmomkit::MomentSourceRateStatus::invalid_input);
    assert(missing_boundary.message == "Moment-boundary contribution callback must be provided.");
    assert(!source_called);
    assert(!boundary_called);
}

void test_public_qmom_moment_source_rates_source_failure_prevents_boundary()
{
    bool boundary_called = false;
    const qmomkit::QmomMomentSourceRateResult wrong_size = qmomkit::evaluate_qmom_moment_source_rates(
        0.0, {1.0, 2.0}, [](double, const qmomkit::QmomSourceContext&) { return std::vector<double>{0.0}; },
        [&boundary_called](double, const std::vector<double>& state)
        {
            boundary_called = true;
            return qmomkit::MomentBoundaryRate{std::vector<double>(state.size(), 0.0)};
        });
    assert(!wrong_size.ok());
    assert(wrong_size.status == qmomkit::MomentSourceRateStatus::invalid_rate);
    assert(wrong_size.message == "Moment-source callback rate values must match the moment state size.");
    assert(wrong_size.rates.empty());
    assert(!boundary_called);

    const qmomkit::QmomMomentSourceRateResult non_finite = qmomkit::evaluate_qmom_moment_source_rates(
        0.0, {1.0, 2.0},
        [](double, const qmomkit::QmomSourceContext& state)
        {
            std::vector<double> rate(state.size(), 0.0);
            rate.front() = std::numeric_limits<double>::infinity();
            return rate;
        },
        [&boundary_called](double, const std::vector<double>& state)
        {
            boundary_called = true;
            return qmomkit::MomentBoundaryRate{std::vector<double>(state.size(), 0.0)};
        });
    assert(!non_finite.ok());
    assert(non_finite.status == qmomkit::MomentSourceRateStatus::invalid_rate);
    assert(non_finite.message == "Moment-source callback rate values must be finite.");
    assert(non_finite.rates.empty());
    assert(!boundary_called);

    const qmomkit::QmomMomentSourceRateResult thrown = qmomkit::evaluate_qmom_moment_source_rates(
        0.0, {1.0}, [](double, const qmomkit::QmomSourceContext&) -> std::vector<double> { throw std::runtime_error("source failed"); },
        [&boundary_called](double, const std::vector<double>& state)
        {
            boundary_called = true;
            return qmomkit::MomentBoundaryRate{std::vector<double>(state.size(), 0.0)};
        });
    assert(!thrown.ok());
    assert(thrown.status == qmomkit::MomentSourceRateStatus::source_failed);
    assert(thrown.message.find("source failed") != std::string::npos);
    assert(thrown.rates.empty());
    assert(!boundary_called);
}

void test_public_qmom_moment_source_rates_boundary_failures()
{
    const std::vector<double> moments{1.0, 2.0};

    const qmomkit::QmomMomentSourceRateResult wrong_size = qmomkit::evaluate_qmom_moment_source_rates(
        0.0, moments, [](double, const qmomkit::QmomSourceContext& state) { return std::vector<double>(state.size(), 0.0); },
        [](double, const std::vector<double>& state) { return qmomkit::MomentBoundaryRate{std::vector<double>(state.size() + 1, 0.0)}; });
    assert(!wrong_size.ok());
    assert(wrong_size.status == qmomkit::MomentSourceRateStatus::invalid_rate);
    assert(wrong_size.message == "Moment-boundary contribution rate values must match the moment state size.");
    assert(wrong_size.rates.empty());

    const qmomkit::QmomMomentSourceRateResult non_finite = qmomkit::evaluate_qmom_moment_source_rates(
        0.0, moments, [](double, const qmomkit::QmomSourceContext& state) { return std::vector<double>(state.size(), 0.0); },
        [](double, const std::vector<double>& state)
        {
            std::vector<double> rate(state.size(), 0.0);
            rate.front() = std::numeric_limits<double>::quiet_NaN();
            return qmomkit::MomentBoundaryRate{rate};
        });
    assert(!non_finite.ok());
    assert(non_finite.status == qmomkit::MomentSourceRateStatus::invalid_rate);
    assert(non_finite.message == "Moment-boundary contribution rate values must be finite.");
    assert(non_finite.rates.empty());

    const qmomkit::QmomMomentSourceRateResult thrown = qmomkit::evaluate_qmom_moment_source_rates(
        0.0, {1.0}, [](double, const qmomkit::QmomSourceContext&) { return std::vector<double>{0.0}; },
        [](double, const std::vector<double>&) -> qmomkit::MomentBoundaryRate { throw std::runtime_error("boundary failed"); });
    assert(!thrown.ok());
    assert(thrown.status == qmomkit::MomentSourceRateStatus::source_failed);
    assert(thrown.message == "Moment-boundary contribution callback failed: boundary failed");
    assert(thrown.rates.empty());
}

void test_public_qmom_moment_source_rates_detect_combined_overflow()
{
    const double maximum = std::numeric_limits<double>::max();
    const qmomkit::QmomMomentSourceRateResult result = qmomkit::evaluate_qmom_moment_source_rates(
        0.0, {1.0}, [maximum](double, const qmomkit::QmomSourceContext&) { return std::vector<double>{maximum}; },
        [maximum](double, const std::vector<double>&) { return qmomkit::MomentBoundaryRate{std::vector<double>{maximum}}; });

    assert(!result.ok());
    assert(result.status == qmomkit::MomentSourceRateStatus::invalid_rate);
    assert(result.message == "Combined moment-source rate values must be finite.");
    assert(result.rates.empty());
}

void test_public_cqmom_moment_source_rates_preserve_layout_and_status()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData moments =
        cqmom_moment_data_from_terminals(layout, {0.2, 0.3, 0.1, 0.4}, {{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}});
    constexpr double evaluation_time = 1.2345;
    std::size_t calls = 0;

    const qmomkit::CqmomMomentSourceRateResult success =
        qmomkit::evaluate_cqmom_moment_source_rates(evaluation_time, moments,
                                                    [&moments, &calls](const double time, const qmomkit::CqmomSourceContext& context)
                                                    {
                                                        ++calls;
                                                        assert_near(time, evaluation_time, 0.0);
                                                        assert_near(context.time(), evaluation_time, 0.0);
                                                        assert(context.size() == moments.values.size());
                                                        std::vector<double> rates(context.size(), 0.0);
                                                        for (std::size_t i = 0; i < context.size(); ++i)
                                                        {
                                                            assert(context.index(i) == moments.layout.required_indices[i]);
                                                            assert_near(context.value(i), moments.values[i], 0.0);
                                                            rates[i] = static_cast<double>(i) + 0.125;
                                                        }
                                                        return qmomkit::CqmomMomentRate{rates};
                                                    });
    assert(success.ok());
    assert(success.status == qmomkit::MomentSourceRateStatus::success);
    assert(success.rates.has_value());
    assert(same_layout(success.rates->layout, layout));
    for (std::size_t i = 0; i < success.rates->values.size(); ++i)
    {
        assert_near(success.rates->values[i], static_cast<double>(i) + 0.125, 0.0);
    }
    assert(calls == 1);

    const qmomkit::CqmomMomentSourceRateResult invalid_rate =
        qmomkit::evaluate_cqmom_moment_source_rates(0.0, moments, [](double, const qmomkit::CqmomSourceContext& context)
                                                    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size() + 1, 0.0)}; });
    assert(!invalid_rate.ok());
    assert(invalid_rate.status == qmomkit::MomentSourceRateStatus::invalid_rate);
    assert(!invalid_rate.rates.has_value());

    const qmomkit::CqmomMomentSourceRateResult non_finite_rate =
        qmomkit::evaluate_cqmom_moment_source_rates(0.0, moments,
                                                    [](double, const qmomkit::CqmomSourceContext& context)
                                                    {
                                                        std::vector<double> rates(context.size(), 0.0);
                                                        rates.front() = std::numeric_limits<double>::quiet_NaN();
                                                        return qmomkit::CqmomMomentRate{rates};
                                                    });
    assert(!non_finite_rate.ok());
    assert(non_finite_rate.status == qmomkit::MomentSourceRateStatus::invalid_rate);
    assert(!non_finite_rate.rates.has_value());

    const qmomkit::CqmomMomentSourceRateResult invalid_input = qmomkit::evaluate_cqmom_moment_source_rates(0.0, moments, {});
    assert(!invalid_input.ok());
    assert(invalid_input.status == qmomkit::MomentSourceRateStatus::invalid_input);
    assert(!invalid_input.rates.has_value());

    const qmomkit::CqmomMomentSourceRateResult source_failed =
        qmomkit::evaluate_cqmom_moment_source_rates(0.0, moments, [](double, const qmomkit::CqmomSourceContext&) -> qmomkit::CqmomMomentRate
                                                    { throw std::runtime_error("CQMOM source failed"); });
    assert(!source_failed.ok());
    assert(source_failed.status == qmomkit::MomentSourceRateStatus::source_failed);
    assert(source_failed.message.find("CQMOM source failed") != std::string::npos);
}

} // namespace

int main()
{
    test_public_qmom_advection_rates_and_status_mapping();
    test_public_cqmom_advection_rates_preserve_layout_and_status();
    test_public_cqmom_advection_rates_preserve_adaptive_layout();
    test_public_qmom_moment_source_rates_and_status_mapping();
    test_public_qmom_moment_source_rates_with_boundary();
    test_public_qmom_moment_source_rates_with_zero_boundary();
    test_public_qmom_moment_source_rates_validate_inputs_before_callbacks();
    test_public_qmom_moment_source_rates_source_failure_prevents_boundary();
    test_public_qmom_moment_source_rates_boundary_failures();
    test_public_qmom_moment_source_rates_detect_combined_overflow();
    test_public_cqmom_moment_source_rates_preserve_layout_and_status();
    return 0;
}
