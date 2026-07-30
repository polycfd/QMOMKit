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

void test_public_qmom_continuous_rates_and_status_mapping()
{
    const std::vector<double> moments = qmom_moments_from_terminals({2.0}, {3.0}, 4);
    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>&) { return std::vector<double>{0.5}; };

    const qmomkit::QmomContinuousMomentRateResult success = qmomkit::evaluate_qmom_continuous_moment_rates(moments, 1, 0.25, velocity);
    assert(success.ok());
    assert(success.status == qmomkit::ContinuousMomentRateStatus::success);
    assert(success.inversion_status == qmomkit::QmomStatus::success);
    assert(success.rates[0] == 0.0);
    assert_vector_near(success.rates, {0.0, 1.0, 6.0, 27.0}, 1.0e-13);

    const qmomkit::QmomContinuousMomentRateResult inversion_failed = qmomkit::evaluate_qmom_continuous_moment_rates(
        {1.0, 0.0, -1.0, 0.0}, 2, 0.0, [](double, const std::vector<double>&) { return std::vector<double>{0.0}; });
    assert(!inversion_failed.ok());
    assert(inversion_failed.status == qmomkit::ContinuousMomentRateStatus::inversion_failed);
    assert(inversion_failed.inversion_status == qmomkit::QmomStatus::not_realizable);
    assert(inversion_failed.rates.empty());

    const qmomkit::QmomContinuousMomentRateResult velocity_failed = qmomkit::evaluate_qmom_continuous_moment_rates(
        moments, 1, 0.0, [](double, const std::vector<double>&) { return std::vector<double>{1.0, 2.0}; });
    assert(!velocity_failed.ok());
    assert(velocity_failed.status == qmomkit::ContinuousMomentRateStatus::velocity_failed);
    assert(velocity_failed.inversion_status == qmomkit::QmomStatus::success);
    assert(velocity_failed.message == "Phase-space velocity callback returned the wrong dimension.");

    const qmomkit::QmomContinuousMomentRateResult invalid_input = qmomkit::evaluate_qmom_continuous_moment_rates({}, 1, 0.0, velocity);
    assert(!invalid_input.ok());
    assert(invalid_input.status == qmomkit::ContinuousMomentRateStatus::invalid_input);
    assert(invalid_input.inversion_status == qmomkit::QmomStatus::invalid_input);
    assert(invalid_input.rates.empty());
}

void test_public_cqmom_continuous_rates_preserve_layout_and_status()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomLeanCase lean_case{
        cqmom_moment_data_from_terminals(layout, {0.2, 0.3, 0.1, 0.4}, {{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}})};
    const qmomkit::PhaseSpaceVelocity zero_velocity = [](double, const std::vector<double>& coordinate)
    { return std::vector<double>(coordinate.size(), 0.0); };

    const qmomkit::CqmomContinuousMomentRateResult success = qmomkit::evaluate_cqmom_continuous_moment_rates(lean_case, 0.5, zero_velocity);
    assert(success.ok());
    assert(success.status == qmomkit::ContinuousMomentRateStatus::success);
    assert(success.inversion_status == qmomkit::CqmomStatus::success);
    assert(success.rates.has_value());
    assert(same_layout(success.rates->layout, layout));
    assert_vector_near(success.rates->values, std::vector<double>(layout.required_indices.size(), 0.0), 0.0);
    assert(!success.diagnostics.empty());

    std::vector<double> invalid_values(layout.required_indices.size(), 0.0);
    invalid_values.front() = -1.0;
    const qmomkit::CqmomLeanCase invalid_case{qmomkit::CqmomMomentData{layout, invalid_values}};
    const qmomkit::CqmomContinuousMomentRateResult inversion_failed =
        qmomkit::evaluate_cqmom_continuous_moment_rates(invalid_case, 0.0, zero_velocity);
    assert(!inversion_failed.ok());
    assert(inversion_failed.status == qmomkit::ContinuousMomentRateStatus::inversion_failed);
    assert(inversion_failed.inversion_status == qmomkit::CqmomStatus::qmom_failed);
    assert(!inversion_failed.rates.has_value());
    assert(!inversion_failed.diagnostics.empty());

    const qmomkit::CqmomContinuousMomentRateResult velocity_failed = qmomkit::evaluate_cqmom_continuous_moment_rates(
        lean_case, 0.0, [](double, const std::vector<double>&) { return std::vector<double>{std::numeric_limits<double>::quiet_NaN(), 0.0}; });
    assert(!velocity_failed.ok());
    assert(velocity_failed.status == qmomkit::ContinuousMomentRateStatus::velocity_failed);
    assert(velocity_failed.inversion_status == qmomkit::CqmomStatus::success);
    assert(velocity_failed.message == "Phase-space velocity callback returned a non-finite value.");

    qmomkit::CqmomLeanCase malformed_case = lean_case;
    malformed_case.moments.values.pop_back();
    const qmomkit::CqmomContinuousMomentRateResult invalid_input =
        qmomkit::evaluate_cqmom_continuous_moment_rates(malformed_case, 0.0, zero_velocity);
    assert(!invalid_input.ok());
    assert(invalid_input.status == qmomkit::ContinuousMomentRateStatus::invalid_input);
    assert(invalid_input.inversion_status == qmomkit::CqmomStatus::invalid_input);
    assert(!invalid_input.rates.has_value());
}

void test_public_cqmom_continuous_rates_preserve_adaptive_layout()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::CqmomLeanCase lean_case{cqmom_moment_data_from_terminals(layout, {0.1, 0.3, 0.6}, {{-1.0, -1.0, 1.0}, {-2.0, 2.0, 3.0}})};
    lean_case.qmom_options = {true, 1.0e-10, 1.0e-10, 1.0e-12};
    lean_case.adaptive_branch_counts = true;

    const qmomkit::CqmomContinuousMomentRateResult result = qmomkit::evaluate_cqmom_continuous_moment_rates(
        lean_case, 0.0, [](double, const std::vector<double>& coordinate) { return std::vector<double>(coordinate.size(), 0.0); });
    assert(result.ok());
    assert(result.rates.has_value());
    assert(same_layout(result.rates->layout, layout));
    assert_vector_near(result.rates->values, std::vector<double>(layout.required_indices.size(), 0.0), 0.0);
}

void test_public_qmom_discontinuous_rates_and_status_mapping()
{
    const std::vector<double> moments{1.0, 2.0, 3.0};
    std::size_t calls = 0;
    const qmomkit::QmomDiscontinuousMomentRateResult success =
        qmomkit::evaluate_qmom_discontinuous_moment_rates(0.25, moments,
                                                          [&calls](const double time, const std::vector<double>& state)
                                                          {
                                                              ++calls;
                                                              assert_near(time, 0.25, 0.0);
                                                              assert_vector_near(state, {1.0, 2.0, 3.0}, 0.0);
                                                              return std::vector<double>{state[2], state[0], state[1]};
                                                          });
    assert(success.ok());
    assert(success.status == qmomkit::DiscontinuousMomentRateStatus::success);
    assert_vector_near(success.rates, {3.0, 1.0, 2.0}, 0.0);
    assert(calls == 1);

    const qmomkit::QmomDiscontinuousMomentRateResult invalid_rate = qmomkit::evaluate_qmom_discontinuous_moment_rates(
        0.0, moments, [](double, const std::vector<double>& state) { return std::vector<double>(state.size() + 1, 0.0); });
    assert(!invalid_rate.ok());
    assert(invalid_rate.status == qmomkit::DiscontinuousMomentRateStatus::invalid_rate);
    assert(invalid_rate.message == "Moment-source callback rate values must match the moment state size.");

    const qmomkit::QmomDiscontinuousMomentRateResult non_finite_rate =
        qmomkit::evaluate_qmom_discontinuous_moment_rates(0.0, moments,
                                                          [](double, const std::vector<double>& state)
                                                          {
                                                              std::vector<double> rates(state.size(), 0.0);
                                                              rates.back() = std::numeric_limits<double>::infinity();
                                                              return rates;
                                                          });
    assert(!non_finite_rate.ok());
    assert(non_finite_rate.status == qmomkit::DiscontinuousMomentRateStatus::invalid_rate);
    assert(non_finite_rate.message == "Moment-source callback rate values must be finite.");

    const qmomkit::QmomDiscontinuousMomentRateResult invalid_input = qmomkit::evaluate_qmom_discontinuous_moment_rates(0.0, {}, {});
    assert(!invalid_input.ok());
    assert(invalid_input.status == qmomkit::DiscontinuousMomentRateStatus::invalid_input);
    assert(invalid_input.rates.empty());

    const qmomkit::QmomDiscontinuousMomentRateResult source_failed = qmomkit::evaluate_qmom_discontinuous_moment_rates(
        0.0, moments, [](double, const std::vector<double>&) -> std::vector<double> { throw std::runtime_error("source failed"); });
    assert(!source_failed.ok());
    assert(source_failed.status == qmomkit::DiscontinuousMomentRateStatus::source_failed);
    assert(source_failed.message.find("source failed") != std::string::npos);
}

void test_public_cqmom_discontinuous_rates_preserve_layout_and_status()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData moments =
        cqmom_moment_data_from_terminals(layout, {0.2, 0.3, 0.1, 0.4}, {{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}});
    std::size_t calls = 0;

    const qmomkit::CqmomDiscontinuousMomentRateResult success =
        qmomkit::evaluate_cqmom_discontinuous_moment_rates(0.5, moments,
                                                           [&moments, &calls](const double time, const qmomkit::CqmomSourceContext& context)
                                                           {
                                                               ++calls;
                                                               assert_near(time, 0.5, 0.0);
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
    assert(success.status == qmomkit::DiscontinuousMomentRateStatus::success);
    assert(success.rates.has_value());
    assert(same_layout(success.rates->layout, layout));
    for (std::size_t i = 0; i < success.rates->values.size(); ++i)
    {
        assert_near(success.rates->values[i], static_cast<double>(i) + 0.125, 0.0);
    }
    assert(calls == 1);

    const qmomkit::CqmomDiscontinuousMomentRateResult invalid_rate =
        qmomkit::evaluate_cqmom_discontinuous_moment_rates(0.0, moments, [](double, const qmomkit::CqmomSourceContext& context)
                                                           { return qmomkit::CqmomMomentRate{std::vector<double>(context.size() + 1, 0.0)}; });
    assert(!invalid_rate.ok());
    assert(invalid_rate.status == qmomkit::DiscontinuousMomentRateStatus::invalid_rate);
    assert(!invalid_rate.rates.has_value());

    const qmomkit::CqmomDiscontinuousMomentRateResult non_finite_rate =
        qmomkit::evaluate_cqmom_discontinuous_moment_rates(0.0, moments,
                                                           [](double, const qmomkit::CqmomSourceContext& context)
                                                           {
                                                               std::vector<double> rates(context.size(), 0.0);
                                                               rates.front() = std::numeric_limits<double>::quiet_NaN();
                                                               return qmomkit::CqmomMomentRate{rates};
                                                           });
    assert(!non_finite_rate.ok());
    assert(non_finite_rate.status == qmomkit::DiscontinuousMomentRateStatus::invalid_rate);
    assert(!non_finite_rate.rates.has_value());

    const qmomkit::CqmomDiscontinuousMomentRateResult invalid_input = qmomkit::evaluate_cqmom_discontinuous_moment_rates(0.0, moments, {});
    assert(!invalid_input.ok());
    assert(invalid_input.status == qmomkit::DiscontinuousMomentRateStatus::invalid_input);
    assert(!invalid_input.rates.has_value());

    const qmomkit::CqmomDiscontinuousMomentRateResult source_failed =
        qmomkit::evaluate_cqmom_discontinuous_moment_rates(0.0, moments, [](double, const qmomkit::CqmomSourceContext&) -> qmomkit::CqmomMomentRate
                                                           { throw std::runtime_error("CQMOM source failed"); });
    assert(!source_failed.ok());
    assert(source_failed.status == qmomkit::DiscontinuousMomentRateStatus::source_failed);
    assert(source_failed.message.find("CQMOM source failed") != std::string::npos);
}

} // namespace

int main()
{
    test_public_qmom_continuous_rates_and_status_mapping();
    test_public_cqmom_continuous_rates_preserve_layout_and_status();
    test_public_cqmom_continuous_rates_preserve_adaptive_layout();
    test_public_qmom_discontinuous_rates_and_status_mapping();
    test_public_cqmom_discontinuous_rates_preserve_layout_and_status();
    return 0;
}
