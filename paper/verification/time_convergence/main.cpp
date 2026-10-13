#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

constexpr double splitting_final_time = 1.2;
constexpr double noncommuting_b = 0.45;
constexpr double noncommuting_gamma = 0.17;

struct SplitRun
{
    std::string mode;
    std::size_t level = 0;
    double step_size = 0.0;
    double maximum_scaled_moment_error = 0.0;
    double observed_order = 0.0;
};

double observed_order(const double coarse_error, const double fine_error, const double coarse_step, const double fine_step)
{
    if (!(coarse_error > 0.0) || !(fine_error > 0.0) || !(coarse_step > fine_step))
    {
        throw std::runtime_error("Observed-order inputs must be strictly positive and refined.");
    }
    return std::log(coarse_error / fine_error) / std::log(coarse_step / fine_step);
}

// Splitting problem: two weighted points move under the affine node equation
// dxi/dt = b - gamma xi, split into the advection dxi/dt = b and the moment
// source dm_k/dt = -gamma k m_k. Two nodes close m_0 to m_3 exactly, so the only
// errors are those of the splitting and of the SSP-RK3 source step.
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

qmomkit::QmomMomentSource noncommuting_source()
{
    return [](double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t order = 1; order < context.size(); ++order)
        {
            rate[order] = -noncommuting_gamma * static_cast<double>(order) * context.value(order);
        }
        return qmomkit::QmomMomentRate{std::move(rate)};
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
    const qmomkit::QmomCase qmom_case{qmomkit::MomentSet{initial}, 2, {}};
    qmomkit::QmomState state = qmom_case.initial_state();
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
        if (!result.ok())
        {
            throw std::runtime_error(std::string("Phase-space evolution failed for split mode ") + split_mode_name(scheme) + ": " + result.message);
        }
        state = result.state;
        time = result.final_time;
    }

    const std::vector<double> exact = noncommuting_exact(initial, splitting_final_time);
    SplitRun run{split_mode_name(scheme), level, step_size};
    for (std::size_t moment = 0; moment < exact.size(); ++moment)
    {
        const double absolute = std::abs(state.moment_values()[moment] - exact[moment]);
        run.maximum_scaled_moment_error = std::max(run.maximum_scaled_moment_error, absolute / std::max(1.0, std::abs(exact[moment])));
    }
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
        // The observed orders on the two finest refinements, levels 3 and 4.
        const double threshold = scheme == qmomkit::PhaseSpaceEvolutionScheme::strang ? 1.85 : 0.85;
        if (runs[first + 3].observed_order <= threshold - 0.1 || runs[first + 4].observed_order <= threshold)
        {
            throw std::runtime_error(std::string("Split mode ") + split_mode_name(scheme) + " did not show its expected asymptotic order.");
        }
    }
    return runs;
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

} // namespace

int main()
{
    try
    {
        write_splitting_csv(run_splitting_convergence());
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
