#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "qmomkit/multi_index.hpp"
#include "test_case_dispatch.hpp"
#include "test_checks.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// A population of a few weighted points is sampled, and the samples are moved
// with evolve_ensemble(...). Every sample must follow the analytical trajectory
// of its point, and the moments estimated from the moved samples must agree,
// within six standard errors, with the exact moments of the weighted points
// moved along those trajectories.
namespace
{

using qmomkit::monte_carlo::Ensemble;
using qmomkit::monte_carlo::Sample;

struct WeightedPoints
{
    std::vector<Sample> points;
    std::vector<double> weights;
};

using Trajectory = std::function<Sample(const Sample&)>;

qmomkit::AdaptiveRkOptions accurate_options()
{
    qmomkit::AdaptiveRkOptions options;
    options.relative_tolerance = 1.0e-10;
    options.absolute_tolerance = 1.0e-12;
    options.initial_step = 0.01;
    options.max_step = 0.05;
    return options;
}

double total_weight(const WeightedPoints& population)
{
    double total = 0.0;
    for (const double weight : population.weights)
    {
        total += weight;
    }
    return total;
}

// Draws each sample from the points with probability proportional to its
// weight, so every sample carries the weight total / sample_count.
Ensemble sample_points(const WeightedPoints& population, const std::size_t sample_count, const unsigned long long seed)
{
    std::mt19937_64 engine(seed);
    std::discrete_distribution<std::size_t> point(population.weights.begin(), population.weights.end());
    Ensemble ensemble(sample_count);
    for (Sample& sample : ensemble)
    {
        sample = population.points[point(engine)];
    }
    return ensemble;
}

double exact_moment(const WeightedPoints& population, const qmomkit::MultiIndex& order)
{
    double moment = 0.0;
    for (std::size_t index = 0; index < population.points.size(); ++index)
    {
        double monomial = population.weights[index];
        for (std::size_t dimension = 0; dimension < order.size(); ++dimension)
        {
            monomial *= std::pow(population.points[index][dimension], static_cast<double>(order[dimension]));
        }
        moment += monomial;
    }
    return moment;
}

void check_sampled_evolution(const WeightedPoints& initial, const Trajectory& trajectory, const qmomkit::PhaseSpaceVelocity& velocity,
                             const double t0, const double t1, const std::size_t sample_count, const unsigned long long seed,
                             const std::vector<std::size_t>& maximum_orders)
{
    const Ensemble samples = sample_points(initial, sample_count, seed);
    const qmomkit::monte_carlo::EnsembleEvolutionResult evolved =
        qmomkit::monte_carlo::evolve_ensemble(samples, velocity, t0, t1, accurate_options());
    CHECK(evolved.ok());
    CHECK(evolved.ensemble.size() == sample_count);
    CHECK_CLOSE(evolved.final_time, t1, 0.0);
    for (std::size_t sample = 0; sample < sample_count; ++sample)
    {
        const Sample expected = trajectory(samples[sample]);
        for (std::size_t dimension = 0; dimension < expected.size(); ++dimension)
        {
            CHECK_CLOSE(evolved.ensemble[sample][dimension], expected[dimension], 1.0e-9);
        }
    }

    WeightedPoints final_points{{}, initial.weights};
    for (const Sample& point : initial.points)
    {
        final_points.points.push_back(trajectory(point));
    }
    const qmomkit::monte_carlo::TensorMomentEstimate estimate =
        qmomkit::monte_carlo::compute_ensemble_moments(evolved.ensemble, total_weight(initial), maximum_orders);
    CHECK(estimate.ok());
    const std::vector<qmomkit::MultiIndex> orders = qmomkit::enumerate_multi_indices(estimate.moments.extents);
    for (std::size_t flat = 0; flat < orders.size(); ++flat)
    {
        const double reference = exact_moment(final_points, orders[flat]);
        const double standard_error = estimate.standard_errors[flat];
        std::ostringstream context;
        context << "flat_index=" << flat << " reference=" << reference << " estimate=" << estimate.moments.values[flat]
                << " standard_error=" << standard_error;
        CHECK_CLOSE_CONTEXT(estimate.moments.values[flat], reference, std::max(6.0 * standard_error, 1.0e-12), 0.0, context.str());
    }
}

void test_constant_velocity_comparison()
{
    const double t0 = 0.1;
    const double t1 = 0.6;
    const double velocity = -0.35;
    check_sampled_evolution(
        WeightedPoints{{{-0.2}, {0.1}, {0.35}}, {0.8, 1.1, 0.9}}, [&](const Sample& x) { return Sample{x[0] + velocity * (t1 - t0)}; },
        qmomkit::make_phase_space_velocity([velocity](double, const std::vector<double>&) { return std::vector<double>{velocity}; }), t0, t1, 6000,
        4101, {2});
}

void test_affine_evolution_comparison()
{
    const double t0 = -0.2;
    const double t1 = 0.4;
    const double rate = -0.25;
    const double shift = 0.8;
    check_sampled_evolution(
        WeightedPoints{{{-0.3}, {0.2}, {0.55}}, {0.7, 1.0, 0.8}},
        [&](const Sample& x) { return Sample{(x[0] + shift / rate) * std::exp(rate * (t1 - t0)) - shift / rate}; },
        qmomkit::make_phase_space_velocity([rate, shift](double, const std::vector<double>& x) { return std::vector<double>{rate * x[0] + shift}; }),
        t0, t1, 8000, 4201, {2});
}

void test_two_dimensional_constant_velocity_comparison()
{
    const double t0 = 0.0;
    const double t1 = 0.25;
    const double cx = 0.12;
    const double cy = -0.08;
    check_sampled_evolution(WeightedPoints{{{-0.2, 0.1}, {0.25, -0.15}, {0.4, 0.3}}, {0.8, 1.0, 0.9}},
                            [&](const Sample& x) { return Sample{x[0] + cx * (t1 - t0), x[1] + cy * (t1 - t0)}; },
                            qmomkit::make_phase_space_velocity([cx, cy](double, const std::vector<double>&) { return std::vector<double>{cx, cy}; }),
                            t0, t1, 7000, 5201, {2, 2});
}

void test_two_dimensional_affine_comparison()
{
    const double t0 = -0.1;
    const double t1 = 0.2;
    const double ax = -0.2;
    const double bx = 0.35;
    const double ay = 0.15;
    const double by = -0.2;
    check_sampled_evolution(
        WeightedPoints{{{-0.15, 0.2}, {0.25, -0.1}, {0.45, 0.35}}, {0.7, 1.1, 0.8}}, [&](const Sample& x)
        { return Sample{(x[0] + bx / ax) * std::exp(ax * (t1 - t0)) - bx / ax, (x[1] + by / ay) * std::exp(ay * (t1 - t0)) - by / ay}; },
        qmomkit::make_phase_space_velocity([ax, bx, ay, by](double, const std::vector<double>& x)
                                           { return std::vector<double>{ax * x[0] + bx, ay * x[1] + by}; }),
        t0, t1, 8000, 5301, {2, 2});
}

void test_two_dimensional_coupled_linear_comparison()
{
    const double t0 = 0.0;
    const double t1 = 0.2;
    const double a = -0.3;
    const double b = 0.4;
    const double c = -0.1;
    check_sampled_evolution(WeightedPoints{{{-0.25, 0.1}, {0.15, -0.2}, {0.35, 0.25}}, {0.9, 0.8, 1.0}},
                            [&](const Sample& x)
                            {
                                const double elapsed = t1 - t0;
                                return Sample{x[0] * std::exp(a * elapsed),
                                              x[1] * std::exp(c * elapsed) + b * x[0] * (std::exp(a * elapsed) - std::exp(c * elapsed)) / (a - c)};
                            },
                            qmomkit::make_phase_space_velocity([a, b, c](double, const std::vector<double>& x)
                                                               { return std::vector<double>{a * x[0], b * x[0] + c * x[1]}; }),
                            t0, t1, 9000, 5401, {2, 2});
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::monte_carlo::test::dispatch_test_case(
        argc, argv,
        {
            {"constant_velocity_comparison", test_constant_velocity_comparison},
            {"affine_evolution_comparison", test_affine_evolution_comparison},
            {"two_dimensional_constant_velocity_comparison", test_two_dimensional_constant_velocity_comparison},
            {"two_dimensional_affine_comparison", test_two_dimensional_affine_comparison},
            {"two_dimensional_coupled_linear_comparison", test_two_dimensional_coupled_linear_comparison},
        });
}
