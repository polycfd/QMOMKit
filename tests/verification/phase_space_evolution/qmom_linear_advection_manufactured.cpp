#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr std::size_t node_count = 3;
constexpr std::size_t highest_order = 10;
constexpr std::size_t phase_space_order_count = 6;
constexpr double final_time = 0.1;
constexpr double source_lambda = 0.3;
constexpr double injection_rate = 0.1;
constexpr double injection_coordinate = 1.0;

struct RefinementRun
{
    double time_step = 0.0;
    double moment_rate_error = 0.0;
    double node_advection_error = 0.0;
    double moment_rate_order = std::numeric_limits<double>::quiet_NaN();
    double node_advection_order = std::numeric_limits<double>::quiet_NaN();
};

struct CombinedSourceRun
{
    double time_step = 0.0;
    double moment_rate_error = 0.0;
    std::array<double, 3> phase_space_errors{};
    double moment_rate_order = std::numeric_limits<double>::quiet_NaN();
    std::array<double, 3> phase_space_orders{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                                             std::numeric_limits<double>::quiet_NaN()};
};

double integer_power(const double value, const std::size_t order)
{
    double result = 1.0;
    for (std::size_t i = 0; i < order; ++i)
    {
        result *= value;
    }
    return result;
}

std::vector<double> manufactured_moments(const std::size_t count)
{
    constexpr std::array<double, node_count> nodes{0.5, 1.5, 3.0};
    constexpr std::array<double, node_count> weights{0.2, 0.3, 0.5};
    std::vector<double> moments(count, 0.0);
    for (std::size_t order = 0; order < count; ++order)
    {
        for (std::size_t node = 0; node < node_count; ++node)
        {
            moments[order] += weights[node] * integer_power(nodes[node], order);
        }
    }
    return moments;
}

qmomkit::QmomOptions fixed_qmom_options()
{
    qmomkit::QmomOptions options;
    options.adaptive = false;
    options.automatic_scaling = false;
    return options;
}

qmomkit::QmomState manufactured_state(const std::size_t moment_count)
{
    return qmomkit::QmomCase{qmomkit::MomentSet(manufactured_moments(moment_count)), node_count, fixed_qmom_options()}.initial_state();
}

qmomkit::Quadrature1D quadrature_from_state(const qmomkit::QmomState& state)
{
    const qmomkit::QmomResult& quadrature = state.quadrature();
    return qmomkit::Quadrature1D{quadrature.abscissas, quadrature.weights, qmomkit::QuadratureStatus::success, quadrature.message};
}

qmomkit::PhaseSpaceVelocity linear_velocity()
{
    return [](const double, const std::vector<double>& coordinate) { return std::vector<double>{coordinate[0]}; };
}

std::vector<double> exact_advection_moments(const std::vector<double>& initial, const double time)
{
    std::vector<double> exact(initial.size(), 0.0);
    for (std::size_t order = 0; order < initial.size(); ++order)
    {
        exact[order] = initial[order] * std::exp(static_cast<double>(order) * time);
    }
    return exact;
}

double scaled_infinity_error(const std::vector<double>& actual, const std::vector<double>& expected)
{
    QMOMKIT_VERIFY_CHECK(actual.size() == expected.size());
    double error = 0.0;
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        error = std::max(error, std::abs(actual[i] - expected[i]) / std::max(1.0, std::abs(expected[i])));
    }
    return error;
}

double observed_order(const double coarse_error, const double fine_error)
{
    QMOMKIT_VERIFY_CHECK(coarse_error > 0.0);
    QMOMKIT_VERIFY_CHECK(fine_error > 0.0);
    return std::log(coarse_error / fine_error) / std::log(2.0);
}

const char* scheme_name(const qmomkit::PhaseSpaceEvolutionScheme scheme)
{
    switch (scheme)
    {
    case qmomkit::PhaseSpaceEvolutionScheme::strang:
        return "strang";
    case qmomkit::PhaseSpaceEvolutionScheme::advsrc:
        return "advsrc";
    case qmomkit::PhaseSpaceEvolutionScheme::srcadv:
        return "srcadv";
    }
    return "unknown";
}

qmomkit::QmomPhaseSpaceEvolutionOptions ssp_rk3_advection_options()
{
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.qmom_options = fixed_qmom_options();
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    return options;
}

qmomkit::QmomPhaseSpaceEvolutionOptions accurate_phase_space_options(const qmomkit::PhaseSpaceEvolutionScheme scheme)
{
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.qmom_options = fixed_qmom_options();
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk45;
    options.node_options.adaptive_options.relative_tolerance = 1.0e-12;
    options.node_options.adaptive_options.absolute_tolerance = 1.0e-14;
    options.node_options.adaptive_options.initial_step = 0.0025;
    options.node_options.adaptive_options.max_step = 0.0025;
    options.scheme = scheme;
    return options;
}

void verify_instantaneous_rate_identity()
{
    const std::vector<double> moments = manufactured_moments(highest_order + 1);
    const qmomkit::QmomAdvectionMomentRateResult result =
        qmomkit::evaluate_qmom_advection_moment_rates(moments, node_count, 0.0, linear_velocity(), fixed_qmom_options());
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.rates.size() == moments.size());
    QMOMKIT_VERIFY_CHECK(result.rates.front() == 0.0);

    double maximum_absolute_error = 0.0;
    double maximum_scaled_error = 0.0;
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        const double exact = static_cast<double>(order) * moments[order];
        const double error = std::abs(result.rates[order] - exact);
        maximum_absolute_error = std::max(maximum_absolute_error, error);
        maximum_scaled_error = std::max(maximum_scaled_error, error / std::max(1.0, std::abs(exact)));
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(result.rates[order], exact, 1.0e-10, 2.0e-12);
    }
    std::cout << std::setprecision(17) << "instantaneous-rate orders=0.." << highest_order << " maximum_absolute_error=" << maximum_absolute_error
              << " maximum_scaled_error=" << maximum_scaled_error << '\n';
}

void verify_direct_phase_space_advection()
{
    qmomkit::QmomState state = manufactured_state(highest_order + 1);
    const qmomkit::Quadrature1D initial_quadrature = quadrature_from_state(state);
    const double time_step = 0.00025;
    const std::size_t steps = static_cast<std::size_t>(std::lround(final_time / time_step));
    const qmomkit::QmomPhaseSpaceEvolutionOptions options = ssp_rk3_advection_options();
    double time = 0.0;
    for (std::size_t step = 0; step < steps; ++step)
    {
        const qmomkit::QmomPhaseSpaceEvolutionResult advanced = qmomkit::advance_qmom_advection(state, time, time_step, linear_velocity(), options);
        QMOMKIT_VERIFY_CHECK(advanced.ok());
        QMOMKIT_VERIFY_NEAR(advanced.final_time, time + time_step, 1.0e-14);
        state = advanced.state;
        time += time_step;
    }

    double maximum_node_error = 0.0;
    double maximum_weight_change = 0.0;
    for (std::size_t node = 0; node < node_count; ++node)
    {
        const double exact_node = initial_quadrature.abscissas[node] * std::exp(final_time);
        maximum_node_error = std::max(maximum_node_error, std::abs(state.quadrature().abscissas[node] - exact_node));
        maximum_weight_change = std::max(maximum_weight_change, std::abs(state.quadrature().weights[node] - initial_quadrature.weights[node]));
    }
    const std::vector<double> exact = exact_advection_moments(manufactured_moments(highest_order + 1), final_time);
    const std::vector<double> reconstructed = qmomkit::reconstruct_raw_moments(quadrature_from_state(state), exact.size());
    const double maximum_moment_error = scaled_infinity_error(reconstructed, exact);

    QMOMKIT_VERIFY_CHECK(maximum_node_error < 2.0e-12);
    QMOMKIT_VERIFY_CHECK(maximum_weight_change == 0.0);
    QMOMKIT_VERIFY_CHECK(maximum_moment_error < 2.0e-11);
    std::cout << std::setprecision(17) << "direct-phase-space end_time=" << final_time << " dt=" << time_step
              << " maximum_node_error=" << maximum_node_error << " maximum_weight_change=" << maximum_weight_change
              << " maximum_scaled_raw_moment_error=" << maximum_moment_error << '\n';
}

std::vector<double> advance_moment_rates_ssp_rk3(const double time_step)
{
    std::vector<double> moments = manufactured_moments(highest_order + 1);
    const std::size_t steps = static_cast<std::size_t>(std::lround(final_time / time_step));
    double time = 0.0;
    for (std::size_t step = 0; step < steps; ++step)
    {
        const qmomkit::SspRk3Result updated =
            qmomkit::ssp_rk3_step(moments, time, time_step,
                                  [](const double stage_time, const std::vector<double>& stage_moments)
                                  {
                                      const qmomkit::QmomAdvectionMomentRateResult rate = qmomkit::evaluate_qmom_advection_moment_rates(
                                          stage_moments, node_count, stage_time, linear_velocity(), fixed_qmom_options());
                                      if (!rate.ok())
                                      {
                                          throw std::runtime_error("public QMOM advection moment-rate evaluation failed: " + rate.message);
                                      }
                                      return rate.rates;
                                  },
                                  {});
        QMOMKIT_VERIFY_CHECK(updated.ok());
        QMOMKIT_VERIFY_NEAR(updated.time_step_used, time_step, 0.0);
        moments = updated.state;
        time += time_step;
    }
    QMOMKIT_VERIFY_NEAR(time, final_time, 1.0e-14);
    return moments;
}

std::vector<double> advance_nodes_ssp_rk3(const double time_step)
{
    qmomkit::QmomState state = manufactured_state(highest_order + 1);
    const std::size_t steps = static_cast<std::size_t>(std::lround(final_time / time_step));
    const qmomkit::QmomPhaseSpaceEvolutionOptions options = ssp_rk3_advection_options();
    double time = 0.0;
    for (std::size_t step = 0; step < steps; ++step)
    {
        const qmomkit::QmomPhaseSpaceEvolutionResult advanced = qmomkit::advance_qmom_advection(state, time, time_step, linear_velocity(), options);
        QMOMKIT_VERIFY_CHECK(advanced.ok());
        state = advanced.state;
        time += time_step;
    }
    return qmomkit::reconstruct_raw_moments(quadrature_from_state(state), highest_order + 1);
}

void verify_cross_route_ssp_rk3_convergence()
{
    // The moment-rate system has eigenvalues k while every node follows x'=x,
    // so equal finite-step error constants are neither expected nor asserted.
    const std::vector<double> exact = exact_advection_moments(manufactured_moments(highest_order + 1), final_time);
    const std::vector<double> time_steps{0.02, 0.01, 0.005, 0.0025};
    std::vector<RefinementRun> runs;
    for (const double time_step : time_steps)
    {
        const double moment_rate_error = scaled_infinity_error(advance_moment_rates_ssp_rk3(time_step), exact);
        const double node_advection_error = scaled_infinity_error(advance_nodes_ssp_rk3(time_step), exact);
        runs.push_back(RefinementRun{time_step, moment_rate_error, node_advection_error});
    }
    for (std::size_t run = 1; run < runs.size(); ++run)
    {
        QMOMKIT_VERIFY_CHECK(runs[run].moment_rate_error < runs[run - 1].moment_rate_error);
        QMOMKIT_VERIFY_CHECK(runs[run].node_advection_error < runs[run - 1].node_advection_error);
        runs[run].moment_rate_order = observed_order(runs[run - 1].moment_rate_error, runs[run].moment_rate_error);
        runs[run].node_advection_order = observed_order(runs[run - 1].node_advection_error, runs[run].node_advection_error);
    }
    QMOMKIT_VERIFY_CHECK(runs[runs.size() - 2].moment_rate_order > 2.8);
    QMOMKIT_VERIFY_CHECK(runs.back().moment_rate_order > 2.8);
    QMOMKIT_VERIFY_CHECK(runs[runs.size() - 2].node_advection_order > 2.8);
    QMOMKIT_VERIFY_CHECK(runs.back().node_advection_order > 2.8);

    std::cout << std::setprecision(17) << "cross-route SSP-RK3, T=" << final_time << '\n';
    std::cout << "dt moment_rate_error node_advection_error moment_rate_order node_advection_order\n";
    for (const RefinementRun& run : runs)
    {
        std::cout << run.time_step << ' ' << run.moment_rate_error << ' ' << run.node_advection_error << ' ' << run.moment_rate_order << ' '
                  << run.node_advection_order << '\n';
    }
}

qmomkit::MomentSource commuting_source()
{
    return [](const double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> rates(context.size(), 0.0);
        for (std::size_t order = 0; order < rates.size(); ++order)
        {
            rates[order] = source_lambda * context.value(order);
        }
        return rates;
    };
}

qmomkit::MomentSource injection_source()
{
    return [](const double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> rates(context.size(), 0.0);
        for (std::size_t order = 0; order < rates.size(); ++order)
        {
            rates[order] = injection_rate * integer_power(injection_coordinate, order);
        }
        return rates;
    };
}

std::vector<double> exact_commuting_solution(const std::vector<double>& initial, const double time)
{
    std::vector<double> exact(initial.size(), 0.0);
    for (std::size_t order = 0; order < initial.size(); ++order)
    {
        exact[order] = initial[order] * std::exp((static_cast<double>(order) + source_lambda) * time);
    }
    return exact;
}

std::vector<double> exact_injection_solution(const std::vector<double>& initial, const double time)
{
    std::vector<double> exact(initial.size(), 0.0);
    exact.front() = initial.front() + injection_rate * time;
    for (std::size_t order = 1; order < initial.size(); ++order)
    {
        const double k = static_cast<double>(order);
        exact[order] =
            initial[order] * std::exp(k * time) + injection_rate * integer_power(injection_coordinate, order) * (std::exp(k * time) - 1.0) / k;
    }
    return exact;
}

std::vector<double> evaluate_combined_moment_rates(const double time, const std::vector<double>& moments, const qmomkit::MomentSource& source)
{
    const qmomkit::QmomAdvectionMomentRateResult advection =
        qmomkit::evaluate_qmom_advection_moment_rates(moments, node_count, time, linear_velocity(), fixed_qmom_options());
    if (!advection.ok())
    {
        throw std::runtime_error("public QMOM advection moment-rate evaluation failed: " + advection.message);
    }

    const qmomkit::QmomMomentSourceRateResult moment_source = qmomkit::evaluate_qmom_moment_source_rates(time, moments, source);
    if (!moment_source.ok())
    {
        throw std::runtime_error("public QMOM moment-source-rate evaluation failed: " + moment_source.message);
    }

    QMOMKIT_VERIFY_CHECK(advection.rates.size() == moment_source.rates.size());
    std::vector<double> total_rates(advection.rates.size(), 0.0);
    for (std::size_t order = 0; order < total_rates.size(); ++order)
    {
        total_rates[order] = advection.rates[order] + moment_source.rates[order];
    }
    return total_rates;
}

std::vector<double> advance_combined_moment_rates_ssp_rk3(const double time_step, const qmomkit::MomentSource& source)
{
    std::vector<double> moments = manufactured_moments(phase_space_order_count);
    const std::size_t steps = static_cast<std::size_t>(std::lround(final_time / time_step));
    double time = 0.0;
    for (std::size_t step = 0; step < steps; ++step)
    {
        const qmomkit::SspRk3Result updated =
            qmomkit::ssp_rk3_step(moments, time, time_step, [&source](const double stage_time, const std::vector<double>& stage_moments)
                                  { return evaluate_combined_moment_rates(stage_time, stage_moments, source); }, {});
        QMOMKIT_VERIFY_CHECK(updated.ok());
        QMOMKIT_VERIFY_NEAR(updated.time_step_used, time_step, 0.0);
        moments = updated.state;
        time += time_step;
    }
    QMOMKIT_VERIFY_NEAR(time, final_time, 1.0e-14);
    return moments;
}

std::vector<double> advance_native_phase_space(const double time_step, const qmomkit::PhaseSpaceEvolutionScheme scheme,
                                               const qmomkit::MomentSource& source, const double end_time = final_time)
{
    qmomkit::QmomState state = manufactured_state(phase_space_order_count);
    const qmomkit::QmomPhaseSpaceEvolutionOptions options = accurate_phase_space_options(scheme);
    const std::size_t steps = static_cast<std::size_t>(std::lround(end_time / time_step));
    double time = 0.0;
    for (std::size_t step = 0; step < steps; ++step)
    {
        const qmomkit::QmomPhaseSpaceEvolutionResult advanced =
            qmomkit::advance_qmom_phase_space(state, node_count, time, time_step, linear_velocity(), source, options);
        QMOMKIT_VERIFY_CHECK(advanced.ok());
        QMOMKIT_VERIFY_CHECK(advanced.diagnostics.moment_source_diagnostics.accepted);
        QMOMKIT_VERIFY_NEAR(advanced.final_time, time + time_step, 1.0e-14);
        state = advanced.state;
        time += time_step;
    }
    QMOMKIT_VERIFY_NEAR(time, end_time, 1.0e-14);
    return state.moments();
}

std::vector<CombinedSourceRun> run_combined_source_refinement(const qmomkit::MomentSource& source, const std::vector<double>& exact)
{
    const std::vector<double> time_steps{0.02, 0.01, 0.005, 0.0025};
    const std::vector<qmomkit::PhaseSpaceEvolutionScheme> schemes{
        qmomkit::PhaseSpaceEvolutionScheme::strang, qmomkit::PhaseSpaceEvolutionScheme::advsrc, qmomkit::PhaseSpaceEvolutionScheme::srcadv};
    std::vector<CombinedSourceRun> runs;
    for (const double time_step : time_steps)
    {
        CombinedSourceRun run;
        run.time_step = time_step;
        run.moment_rate_error = scaled_infinity_error(advance_combined_moment_rates_ssp_rk3(time_step, source), exact);
        for (std::size_t scheme_index = 0; scheme_index < schemes.size(); ++scheme_index)
        {
            run.phase_space_errors[scheme_index] = scaled_infinity_error(advance_native_phase_space(time_step, schemes[scheme_index], source), exact);
        }
        runs.push_back(run);
    }
    for (std::size_t run = 1; run < runs.size(); ++run)
    {
        QMOMKIT_VERIFY_CHECK(runs[run].moment_rate_error < runs[run - 1].moment_rate_error);
        runs[run].moment_rate_order = observed_order(runs[run - 1].moment_rate_error, runs[run].moment_rate_error);
        for (std::size_t scheme_index = 0; scheme_index < schemes.size(); ++scheme_index)
        {
            QMOMKIT_VERIFY_CHECK(runs[run].phase_space_errors[scheme_index] < runs[run - 1].phase_space_errors[scheme_index]);
            runs[run].phase_space_orders[scheme_index] =
                observed_order(runs[run - 1].phase_space_errors[scheme_index], runs[run].phase_space_errors[scheme_index]);
        }
    }
    return runs;
}

void print_combined_source_runs(const std::string& title, const std::vector<CombinedSourceRun>& runs)
{
    std::cout << std::setprecision(17) << title << '\n';
    std::cout << "dt rate_error strang_error advsrc_error srcadv_error rate_order strang_order advsrc_order srcadv_order\n";
    for (const CombinedSourceRun& run : runs)
    {
        std::cout << run.time_step << ' ' << run.moment_rate_error;
        for (const double error : run.phase_space_errors)
        {
            std::cout << ' ' << error;
        }
        std::cout << ' ' << run.moment_rate_order;
        for (const double order : run.phase_space_orders)
        {
            std::cout << ' ' << order;
        }
        std::cout << '\n';
    }
}

void verify_commuting_orders(const std::vector<CombinedSourceRun>& runs)
{
    for (const CombinedSourceRun& run : runs)
    {
        for (const double error : run.phase_space_errors)
        {
            QMOMKIT_VERIFY_CHECK(error < 2.0e-7);
        }
    }
    for (std::size_t run = runs.size() - 2; run < runs.size(); ++run)
    {
        QMOMKIT_VERIFY_CHECK(runs[run].moment_rate_order > 2.8);
        for (const double order : runs[run].phase_space_orders)
        {
            QMOMKIT_VERIFY_CHECK(order > 2.8);
        }
    }
}

void verify_noncommuting_orders(const std::vector<CombinedSourceRun>& runs)
{
    for (std::size_t run = runs.size() - 2; run < runs.size(); ++run)
    {
        QMOMKIT_VERIFY_CHECK(runs[run].moment_rate_order > 2.8);
        QMOMKIT_VERIFY_CHECK(runs[run].phase_space_orders[0] > 1.7);
        QMOMKIT_VERIFY_CHECK(runs[run].phase_space_orders[1] > 0.7);
        QMOMKIT_VERIFY_CHECK(runs[run].phase_space_orders[1] < 1.3);
        QMOMKIT_VERIFY_CHECK(runs[run].phase_space_orders[2] > 0.7);
        QMOMKIT_VERIFY_CHECK(runs[run].phase_space_orders[2] < 1.3);
    }
}

std::vector<double> exact_lie_composition(const std::vector<double>& initial, const double time_step, const bool source_then_advection)
{
    std::vector<double> exact(initial.size(), 0.0);
    for (std::size_t order = 0; order < initial.size(); ++order)
    {
        const double advection_factor = std::exp(static_cast<double>(order) * time_step);
        const double injection = time_step * injection_rate * integer_power(injection_coordinate, order);
        exact[order] = source_then_advection ? advection_factor * (initial[order] + injection) : advection_factor * initial[order] + injection;
    }
    return exact;
}

void verify_lie_composition_semantics()
{
    constexpr double time_step = 0.02;
    const std::vector<double> initial = manufactured_moments(phase_space_order_count);
    const std::vector<double> advsrc =
        advance_native_phase_space(time_step, qmomkit::PhaseSpaceEvolutionScheme::advsrc, injection_source(), time_step);
    const std::vector<double> srcadv =
        advance_native_phase_space(time_step, qmomkit::PhaseSpaceEvolutionScheme::srcadv, injection_source(), time_step);
    const std::vector<double> advsrc_expected = exact_lie_composition(initial, time_step, false);
    const std::vector<double> srcadv_expected = exact_lie_composition(initial, time_step, true);
    const double advsrc_error = scaled_infinity_error(advsrc, advsrc_expected);
    const double srcadv_error = scaled_infinity_error(srcadv, srcadv_expected);
    const double advsrc_swapped_error = scaled_infinity_error(advsrc, srcadv_expected);
    const double srcadv_swapped_error = scaled_infinity_error(srcadv, advsrc_expected);

    QMOMKIT_VERIFY_CHECK(advsrc_error < 5.0e-11);
    QMOMKIT_VERIFY_CHECK(srcadv_error < 5.0e-11);
    QMOMKIT_VERIFY_CHECK(advsrc_swapped_error > 1.0e-5);
    QMOMKIT_VERIFY_CHECK(srcadv_swapped_error > 1.0e-5);
    std::cout << std::setprecision(17) << "Lie one-step h=" << time_step << " advsrc=B(A(M)) error=" << advsrc_error
              << " swapped_error=" << advsrc_swapped_error << " srcadv=A(B(M)) error=" << srcadv_error << " swapped_error=" << srcadv_swapped_error
              << '\n';
}

void verify_combined_moment_rates_and_phase_space_source()
{
    const std::vector<double> initial = manufactured_moments(phase_space_order_count);
    const std::vector<CombinedSourceRun> commuting =
        run_combined_source_refinement(commuting_source(), exact_commuting_solution(initial, final_time));
    verify_commuting_orders(commuting);

    const std::vector<CombinedSourceRun> noncommuting =
        run_combined_source_refinement(injection_source(), exact_injection_solution(initial, final_time));
    verify_noncommuting_orders(noncommuting);
    verify_lie_composition_semantics();

    print_combined_source_runs("combined commuting source S_k=lambda M_k", commuting);
    print_combined_source_runs("combined noncommuting injection S_k=J chi_s^k", noncommuting);
}

} // namespace

int main(const int argc, char** argv)
{
    if (argc != 2)
    {
        std::cerr << "Usage: qmom_linear_advection_manufactured <rate_identity|direct_advection|cross_route|combined_source>\n";
        return 2;
    }
    const std::string selector = argv[1];
    if (selector == "rate_identity")
    {
        verify_instantaneous_rate_identity();
    }
    else if (selector == "direct_advection")
    {
        verify_direct_phase_space_advection();
    }
    else if (selector == "cross_route")
    {
        verify_cross_route_ssp_rk3_convergence();
    }
    else if (selector == "combined_source")
    {
        verify_combined_moment_rates_and_phase_space_source();
    }
    else
    {
        std::cerr << "Unknown selector '" << selector << "'.\n";
        return 2;
    }
    return 0;
}
