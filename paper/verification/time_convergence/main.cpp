#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr double final_time = 1.0;
constexpr double splitting_final_time = 1.2;
constexpr double noncommuting_b = 0.45;
constexpr double noncommuting_gamma = 0.17;

struct RkRun
{
    std::string scheme;
    std::size_t level = 0;
    double step_size = 0.0;
    double absolute_error = 0.0;
    double observed_order = 0.0;
};

struct SplitRun
{
    std::string mode;
    std::size_t level = 0;
    double step_size = 0.0;
    double maximum_scaled_moment_error = 0.0;
    double observed_order = 0.0;
};

qmomkit::OdeRhs exponential_decay_rhs()
{
    return [](double, const std::vector<double>& state, std::vector<double>& derivative) { derivative[0] = -state[0]; };
}

double observed_order(const double coarse_error, const double fine_error, const double coarse_step, const double fine_step)
{
    if (!(coarse_error > 0.0) || !(fine_error > 0.0) || !(coarse_step > fine_step))
    {
        throw std::runtime_error("Observed-order inputs must be strictly positive and refined.");
    }
    return std::log(coarse_error / fine_error) / std::log(coarse_step / fine_step);
}

void verify_finite(const double value, const std::string& description)
{
    if (!std::isfinite(value))
    {
        throw std::runtime_error(description + " is not finite.");
    }
}

std::vector<RkRun> run_ssp_rk3()
{
    // y'=-y has the exact solution exp(-t), so the reported final error is an
    // independent measure of the native SSP-RK3 integration error.
    const std::vector<std::size_t> step_counts{8, 16, 32, 64, 128};
    const double exact = std::exp(-final_time);
    std::vector<RkRun> runs;
    runs.reserve(step_counts.size());
    for (std::size_t level = 0; level < step_counts.size(); ++level)
    {
        const std::size_t step_count = step_counts[level];
        const double step_size = final_time / static_cast<double>(step_count);
        double time = 0.0;
        std::vector<double> state{1.0};
        for (std::size_t step = 0; step < step_count; ++step)
        {
            const qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step(state, time, step_size, exponential_decay_rhs(), {});
            if (!result.ok() || result.state.size() != 1 || !std::isfinite(result.state[0]))
            {
                throw std::runtime_error("SSP-RK3 fixed-step advancement failed.");
            }
            state = result.state;
            time += step_size;
        }
        if (std::abs(time - final_time) > 1.0e-14)
        {
            throw std::runtime_error("SSP-RK3 did not reach the requested final time.");
        }
        const double error = std::abs(state[0] - exact);
        verify_finite(error, "SSP-RK3 final error");
        runs.push_back({"ssp_rk3", level, step_size, error, 0.0});
    }
    for (std::size_t level = 1; level < runs.size(); ++level)
    {
        runs[level].observed_order =
            observed_order(runs[level - 1].absolute_error, runs[level].absolute_error, runs[level - 1].step_size, runs[level].step_size);
    }
    if (runs[runs.size() - 2].observed_order <= 2.85 || runs.back().observed_order <= 2.85)
    {
        throw std::runtime_error("SSP-RK3 did not demonstrate third-order fixed-step convergence.");
    }
    return runs;
}

std::vector<RkRun> run_embedded_method(const qmomkit::EmbeddedRkMethod method, const std::string& name, const double expected_order)
{
    const std::vector<std::size_t> step_counts{8, 16, 32, 64, 128};
    const double exact = std::exp(-final_time);
    std::vector<RkRun> runs;
    runs.reserve(step_counts.size());
    for (std::size_t level = 0; level < step_counts.size(); ++level)
    {
        const std::size_t step_count = step_counts[level];
        const double step_size = final_time / static_cast<double>(step_count);
        qmomkit::AdaptiveRkOptions options;
        options.initial_step = step_size;
        options.min_step = step_size;
        options.max_step = step_size;
        options.relative_tolerance = 1.0e-15;
        options.absolute_tolerance = 1.0e-15;
        const qmomkit::AdaptiveRkResult result = qmomkit::integrate_embedded_rk({1.0}, 0.0, final_time, exponential_decay_rhs(), method, options);
        if (result.status != qmomkit::AdaptiveRkStatus::success || result.state.size() != 1 || result.accepted_steps != step_count ||
            !std::isfinite(result.state[0]) || std::abs(result.t - final_time) > 1.0e-14)
        {
            throw std::runtime_error(name + " fixed-step advancement failed.");
        }
        const double error = std::abs(result.state[0] - exact);
        verify_finite(error, name + " final error");
        runs.push_back({name, level, step_size, error, 0.0});
    }
    for (std::size_t level = 1; level < runs.size(); ++level)
    {
        runs[level].observed_order =
            observed_order(runs[level - 1].absolute_error, runs[level].absolute_error, runs[level - 1].step_size, runs[level].step_size);
    }
    if (runs[runs.size() - 2].observed_order <= expected_order - 0.2 || runs.back().observed_order <= expected_order - 0.2)
    {
        throw std::runtime_error(name + " did not demonstrate its expected fixed-step convergence order.");
    }
    return runs;
}

std::vector<double> noncommuting_exact(const std::vector<double>& initial, const double elapsed)
{
    // The triangular exact solution is independent of the split evolution
    // implementation and supplies the moment reference for all three modes.
    const std::size_t count = initial.size();
    std::vector<double> diagonal(count, 0.0);
    for (std::size_t order = 0; order < count; ++order)
    {
        diagonal[order] = -noncommuting_gamma * static_cast<double>(order);
    }
    std::vector<std::vector<double>> coefficients(count, std::vector<double>(count, 0.0));
    coefficients[0][0] = initial[0];
    for (std::size_t order = 1; order < count; ++order)
    {
        double lower_sum = 0.0;
        for (std::size_t basis = 0; basis < order; ++basis)
        {
            coefficients[order][basis] =
                static_cast<double>(order) * noncommuting_b * coefficients[order - 1][basis] / (diagonal[basis] - diagonal[order]);
            lower_sum += coefficients[order][basis];
        }
        coefficients[order][order] = initial[order] - lower_sum;
    }

    std::vector<double> exact(count, 0.0);
    for (std::size_t order = 0; order < count; ++order)
    {
        for (std::size_t basis = 0; basis <= order; ++basis)
        {
            exact[order] += coefficients[order][basis] * std::exp(diagonal[basis] * elapsed);
        }
    }
    return exact;
}

qmomkit::PhaseSpaceVelocity noncommuting_velocity()
{
    return qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{noncommuting_b}; });
}

qmomkit::MomentSource noncommuting_source()
{
    return [](double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t order = 1; order < context.size(); ++order)
        {
            rate[order] = -noncommuting_gamma * static_cast<double>(order) * context.value(order);
        }
        return rate;
    };
}

const char* split_mode_name(const qmomkit::PhaseSpaceEvolutionScheme scheme)
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

SplitRun run_split_case(const qmomkit::PhaseSpaceEvolutionScheme scheme, const std::size_t level, const double step_size)
{
    const std::vector<double> initial{1.0, 0.5, 1.2399999999999998, 1.4119999999999997};
    const std::size_t step_count = static_cast<std::size_t>(std::lround(splitting_final_time / step_size));
    const qmomkit::QmomCase case_{qmomkit::MomentSet{initial}, 2, {}};
    qmomkit::QmomState state = case_.initial_state();
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk45;
    options.node_options.adaptive_options.relative_tolerance = 1.0e-13;
    options.node_options.adaptive_options.absolute_tolerance = 1.0e-15;
    options.node_options.adaptive_options.initial_step = 0.0025;
    options.node_options.adaptive_options.max_step = 0.0025;
    options.scheme = scheme;

    double time = 0.0;
    for (std::size_t step = 0; step < step_count; ++step)
    {
        const qmomkit::QmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_qmom_phase_space(state, time, step_size, noncommuting_velocity(), noncommuting_source(), options);
        if (!result.ok() || result.state.active_node_count() != 2 ||
            !std::all_of(result.state.moments().begin(), result.state.moments().end(), [](const double value) { return std::isfinite(value); }))
        {
            throw std::runtime_error(std::string("Phase-space evolution failed for split mode ") + split_mode_name(scheme) + ".");
        }
        state = result.state;
        time = result.final_time;
    }
    if (std::abs(time - splitting_final_time) > 1.0e-12)
    {
        throw std::runtime_error("Phase-space evolution did not reach the requested final time.");
    }

    const std::vector<double> exact = noncommuting_exact(initial, splitting_final_time);
    SplitRun run{split_mode_name(scheme), level, step_size};
    for (std::size_t moment = 0; moment < exact.size(); ++moment)
    {
        const double absolute = std::abs(state.moments()[moment] - exact[moment]);
        run.maximum_scaled_moment_error = std::max(run.maximum_scaled_moment_error, absolute / std::max(1.0, std::abs(exact[moment])));
    }
    verify_finite(run.maximum_scaled_moment_error, "Maximum scaled moment error");
    return run;
}

std::vector<SplitRun> run_splitting_convergence()
{
    const std::vector<double> step_sizes{0.3, 0.15, 0.075, 0.0375, 0.01875};
    std::vector<SplitRun> runs;
    for (const qmomkit::PhaseSpaceEvolutionScheme scheme :
         {qmomkit::PhaseSpaceEvolutionScheme::advsrc, qmomkit::PhaseSpaceEvolutionScheme::srcadv, qmomkit::PhaseSpaceEvolutionScheme::strang})
    {
        const std::size_t first = runs.size();
        for (std::size_t level = 0; level < step_sizes.size(); ++level)
        {
            runs.push_back(run_split_case(scheme, level, step_sizes[level]));
            if (level > 0)
            {
                runs.back().observed_order = observed_order(runs[runs.size() - 2].maximum_scaled_moment_error,
                                                            runs.back().maximum_scaled_moment_error, step_sizes[level - 1], step_sizes[level]);
            }
        }
        const double threshold = scheme == qmomkit::PhaseSpaceEvolutionScheme::strang ? 1.85 : 0.85;
        if (runs[first + 3].observed_order <= threshold - 0.1 || runs[first + 4].observed_order <= threshold)
        {
            throw std::runtime_error(std::string("Split mode ") + split_mode_name(scheme) + " did not show its expected asymptotic order.");
        }
    }
    return runs;
}

void write_rk_csv(const std::vector<RkRun>& runs)
{
    std::ofstream output("rk_convergence.csv");
    if (!output)
    {
        throw std::runtime_error("Could not open rk_convergence.csv for writing.");
    }
    output << "scheme,level,step_size,absolute_error,observed_order\n";
    output << std::scientific << std::setprecision(17);
    for (const RkRun& run : runs)
    {
        output << run.scheme << ',' << run.level << ',' << run.step_size << ',' << run.absolute_error << ',' << run.observed_order << '\n';
    }
}

void write_splitting_csv(const std::vector<SplitRun>& runs)
{
    std::ofstream output("splitting_convergence.csv");
    if (!output)
    {
        throw std::runtime_error("Could not open splitting_convergence.csv for writing.");
    }
    output << "split_mode,level,outer_step_size,max_scaled_moment_error,observed_order\n";
    output << std::scientific << std::setprecision(17);
    for (const SplitRun& run : runs)
    {
        output << run.mode << ',' << run.level << ',' << run.step_size << ',' << run.maximum_scaled_moment_error << ',' << run.observed_order << '\n';
    }
}

void print_rk_table(const std::vector<RkRun>& runs)
{
    std::cout << "Runge--Kutta fixed-step convergence\n";
    std::cout << std::left << std::setw(10) << "scheme" << std::right << std::setw(7) << "level" << std::setw(16) << "step" << std::setw(18)
              << "absolute error" << std::setw(14) << "order" << '\n';
    std::cout << std::scientific << std::setprecision(6);
    for (const RkRun& run : runs)
    {
        std::cout << std::left << std::setw(10) << run.scheme << std::right << std::setw(7) << run.level << std::setw(16) << run.step_size
                  << std::setw(18) << run.absolute_error << std::setw(14) << run.observed_order << '\n';
    }
}

void print_splitting_table(const std::vector<SplitRun>& runs)
{
    std::cout << "Operator-splitting convergence\n";
    std::cout << std::left << std::setw(10) << "split mode" << std::right << std::setw(7) << "level" << std::setw(16) << "outer step" << std::setw(22)
              << "max scaled moment" << std::setw(14) << "order" << '\n';
    std::cout << std::scientific << std::setprecision(6);
    for (const SplitRun& run : runs)
    {
        std::cout << std::left << std::setw(10) << run.mode << std::right << std::setw(7) << run.level << std::setw(16) << run.step_size
                  << std::setw(22) << run.maximum_scaled_moment_error << std::setw(14) << run.observed_order << '\n';
    }
}

} // namespace

int main()
{
    try
    {
        std::vector<RkRun> rk_runs = run_ssp_rk3();
        const std::vector<RkRun> rk23_runs = run_embedded_method(qmomkit::EmbeddedRkMethod::rk23, "rk23", 3.0);
        const std::vector<RkRun> rk45_runs = run_embedded_method(qmomkit::EmbeddedRkMethod::rk45, "rk45", 5.0);
        rk_runs.insert(rk_runs.end(), rk23_runs.begin(), rk23_runs.end());
        rk_runs.insert(rk_runs.end(), rk45_runs.begin(), rk45_runs.end());
        const std::vector<SplitRun> splitting_runs = run_splitting_convergence();

        write_rk_csv(rk_runs);
        write_splitting_csv(splitting_runs);
        print_rk_table(rk_runs);
        print_splitting_table(splitting_runs);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
