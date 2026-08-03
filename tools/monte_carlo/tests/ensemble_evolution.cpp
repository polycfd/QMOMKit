#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "monte_carlo/sampling.hpp"
#include "qmomkit/distributions.hpp"
#include "test_checks.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using qmomkit::monte_carlo::Ensemble;
using qmomkit::monte_carlo::EnsembleEvolutionStatus;
using qmomkit::monte_carlo::Sample;
using qmomkit::monte_carlo::SamplingOptions;

qmomkit::AdaptiveRkOptions accurate_options()
{
    qmomkit::AdaptiveRkOptions options;
    options.relative_tolerance = 1.0e-10;
    options.absolute_tolerance = 1.0e-12;
    options.initial_step = 0.01;
    options.max_step = 0.01;
    return options;
}

void check_close(const double actual, const double expected, const double tolerance)
{
    CHECK_CLOSE(actual, expected, tolerance);
}

void check_contains(const std::string& text, const std::string& needle)
{
    CHECK(text.find(needle) != std::string::npos);
}

void check_ensemble_equal(const Ensemble& lhs, const Ensemble& rhs)
{
    CHECK(lhs.size() == rhs.size());
    for (std::size_t sample = 0; sample < lhs.size(); ++sample)
    {
        CHECK(lhs[sample] == rhs[sample]);
    }
}

void test_single_sample_constant_velocity()
{
    const Ensemble ensemble{{1.25}};
    const double t0 = 0.0;
    const double t1 = 0.75;
    const double velocity = -0.4;
    const qmomkit::PhaseSpaceVelocity callback = [velocity](double, const std::vector<double>& coordinate)
    {
        CHECK(coordinate.size() == 1);
        return std::vector<double>{velocity};
    };

    const qmomkit::monte_carlo::EnsembleEvolutionResult result =
        qmomkit::monte_carlo::evolve_ensemble(ensemble, callback, t0, t1, accurate_options());

    CHECK(result.ok());
    CHECK(result.ensemble.size() == 1);
    CHECK(ensemble[0][0] == 1.25);
    check_close(result.ensemble[0][0], 1.25 + velocity * (t1 - t0), 1.0e-12);
    check_close(result.final_time, t1, 0.0);
}

void test_convenience_evolution_uses_default_adaptive_initial_step()
{
    const double time_step = 2.0 * qmomkit::default_adaptive_initial_step;
    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>&) { return std::vector<double>{-0.4}; };

    const qmomkit::monte_carlo::EnsembleEvolutionResult ensemble = qmomkit::monte_carlo::evolve_ensemble(Ensemble{{1.25}}, velocity, 0.0, time_step);
    CHECK(ensemble.ok());
    check_close(ensemble.ensemble[0][0], 1.25 - 0.4 * time_step, 1.0e-15);
    CHECK(ensemble.total_accepted_steps == 2);

    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {1.25};
    quadrature.weights = {1.0};
    const qmomkit::QmomNodeAdvectionResult nodes = qmomkit::monte_carlo::evolve_quadrature_nodes(quadrature, velocity, 0.0, time_step);
    CHECK(nodes.ok());
    check_close(nodes.quadrature.abscissas[0], 1.25 - 0.4 * time_step, 1.0e-15);
    CHECK(nodes.diagnostics.step_count == 2);
}

void test_single_sample_linear_growth()
{
    const Ensemble ensemble{{0.8}};
    const double t0 = 0.0;
    const double t1 = 0.6;
    const double rate = 0.7;
    const qmomkit::PhaseSpaceVelocity callback = [rate](double, const std::vector<double>& coordinate)
    { return std::vector<double>{rate * coordinate[0]}; };

    const qmomkit::monte_carlo::EnsembleEvolutionResult result =
        qmomkit::monte_carlo::evolve_ensemble(ensemble, callback, t0, t1, accurate_options());

    CHECK(result.ok());
    check_close(result.ensemble[0][0], ensemble[0][0] * std::exp(rate * (t1 - t0)), 1.0e-9);
}

void test_deterministic_ensemble_rk23_max_step_convergence()
{
    const Ensemble initial{{0.5}, {1.0}, {2.0}};
    const std::vector<double> step_sizes{0.25, 0.125, 0.0625, 0.03125};
    const qmomkit::PhaseSpaceVelocity callback = [](double, const std::vector<double>& coordinate)
    {
        CHECK(coordinate.size() == 1);
        return std::vector<double>{coordinate[0]};
    };

    std::vector<double> errors;
    errors.reserve(step_sizes.size());
    for (const double step_size : step_sizes)
    {
        qmomkit::AdaptiveRkOptions options;
        options.relative_tolerance = 1.0e-2;
        options.absolute_tolerance = 1.0e-4;
        options.initial_step = step_size;
        options.max_step = step_size;
        options.min_step = 0.01 * step_size;
        options.max_steps = 1000;

        const qmomkit::monte_carlo::EnsembleEvolutionResult result = qmomkit::monte_carlo::evolve_ensemble(initial, callback, 0.0, 1.0, options);
        CHECK(result.ok());
        CHECK(result.ensemble.size() == initial.size());
        CHECK(result.total_accepted_steps == initial.size() * static_cast<std::size_t>(std::lround(1.0 / step_size)));
        CHECK(result.total_rejected_steps == 0);
        check_close(result.final_time, 1.0, 0.0);

        double maximum_error = 0.0;
        for (std::size_t sample = 0; sample < initial.size(); ++sample)
        {
            CHECK(result.ensemble[sample].size() == initial[sample].size());
            maximum_error = std::max(maximum_error, std::abs(result.ensemble[sample][0] - initial[sample][0] * std::exp(1.0)));
        }
        errors.push_back(maximum_error);
    }

    CHECK(initial.size() == 3);
    CHECK(initial[0][0] == 0.5);
    CHECK(initial[1][0] == 1.0);
    CHECK(initial[2][0] == 2.0);
    CHECK(errors.front() > 0.0);
    for (std::size_t refinement = 1; refinement < errors.size(); ++refinement)
    {
        CHECK(errors[refinement] < errors[refinement - 1]);
        const double observed_order = std::log(errors[refinement - 1] / errors[refinement]) / std::log(2.0);
        CHECK(observed_order > 2.0);
    }
}

void test_single_sample_affine_evolution()
{
    const Ensemble ensemble{{1.4}};
    const double t0 = 0.0;
    const double t1 = 0.5;
    const double rate = -0.35;
    const double shift = 0.9;
    const qmomkit::PhaseSpaceVelocity callback = [rate, shift](double, const std::vector<double>& coordinate)
    { return std::vector<double>{rate * coordinate[0] + shift}; };

    const qmomkit::monte_carlo::EnsembleEvolutionResult result =
        qmomkit::monte_carlo::evolve_ensemble(ensemble, callback, t0, t1, accurate_options());

    const double elapsed = t1 - t0;
    const double expected = (ensemble[0][0] + shift / rate) * std::exp(rate * elapsed) - shift / rate;
    CHECK(result.ok());
    check_close(result.ensemble[0][0], expected, 1.0e-9);
}

void test_ensemble_evolution_and_moments()
{
    const qmomkit::IndependentProductDistribution distribution{{qmomkit::NormalDistribution{0.25, 0.3}, qmomkit::GammaDistribution{2.0, 0.4}}, 2.5};
    const qmomkit::monte_carlo::SamplingResult sampling =
        qmomkit::monte_carlo::sample_independent_product(distribution, SamplingOptions{2500, 11001});
    CHECK(sampling.ok());

    const double t0 = 0.1;
    const double t1 = 0.4;
    const qmomkit::PhaseSpaceVelocity callback = [](double, const std::vector<double>& coordinate)
    {
        CHECK(coordinate.size() == 2);
        return std::vector<double>{0.2 + 0.1 * coordinate[0], -0.15 * coordinate[1]};
    };

    const qmomkit::monte_carlo::EnsembleEvolutionResult result =
        qmomkit::monte_carlo::evolve_ensemble(sampling.ensemble, callback, t0, t1, accurate_options());
    CHECK(result.ok());
    CHECK(result.ensemble.size() == sampling.ensemble.size());
    CHECK(result.total_accepted_steps > 0);
    CHECK(result.total_function_evaluations > 0);
    CHECK(result.maximum_rejected_steps <= result.total_rejected_steps);
    check_close(result.final_time, t1, 0.0);

    for (const Sample& sample : result.ensemble)
    {
        CHECK(sample.size() == 2);
        CHECK_FINITE(sample[0]);
        CHECK_FINITE(sample[1]);
    }

    const qmomkit::monte_carlo::TensorMomentEstimate moments =
        qmomkit::monte_carlo::compute_ensemble_moments(result.ensemble, sampling.number_density, std::vector<std::size_t>{2, 2});
    CHECK(moments.ok());
    check_close(moments.moments.values.front(), sampling.number_density, 0.0);
    for (const double value : moments.moments.values)
    {
        CHECK_FINITE(value);
    }
    for (const double error : moments.standard_errors)
    {
        CHECK_FINITE(error);
    }
}

void test_failure_inputs_preserve_original()
{
    const Ensemble original{{1.0}, {2.0}};
    const qmomkit::PhaseSpaceVelocity identity = [](double, const std::vector<double>& coordinate) { return coordinate; };
    const qmomkit::AdaptiveRkOptions options = accurate_options();

    const qmomkit::monte_carlo::EnsembleEvolutionResult empty = qmomkit::monte_carlo::evolve_ensemble(Ensemble{}, identity, 0.0, 1.0, options);
    CHECK(!empty.ok());
    CHECK(empty.status == EnsembleEvolutionStatus::invalid_input);
    CHECK(empty.ensemble.empty());

    const qmomkit::monte_carlo::EnsembleEvolutionResult inconsistent =
        qmomkit::monte_carlo::evolve_ensemble(Ensemble{{1.0}, {1.0, 2.0}}, identity, 0.0, 1.0, options);
    CHECK(!inconsistent.ok());
    CHECK(inconsistent.failed_sample_index == 1);
    CHECK(inconsistent.ensemble.empty());

    const qmomkit::monte_carlo::EnsembleEvolutionResult nonfinite =
        qmomkit::monte_carlo::evolve_ensemble(Ensemble{{std::numeric_limits<double>::quiet_NaN()}}, identity, 0.0, 1.0, options);
    CHECK(!nonfinite.ok());
    CHECK(nonfinite.failed_sample_index == 0);

    const qmomkit::monte_carlo::EnsembleEvolutionResult invalid_time = qmomkit::monte_carlo::evolve_ensemble(original, identity, 0.0, 0.0, options);
    CHECK(!invalid_time.ok());
    CHECK(invalid_time.status == EnsembleEvolutionStatus::invalid_input);

    const qmomkit::PhaseSpaceVelocity missing;
    const qmomkit::monte_carlo::EnsembleEvolutionResult missing_callback =
        qmomkit::monte_carlo::evolve_ensemble(original, missing, 0.0, 1.0, options);
    CHECK(!missing_callback.ok());
    CHECK(missing_callback.status == EnsembleEvolutionStatus::invalid_input);

    CHECK(original[0][0] == 1.0);
    CHECK(original[1][0] == 2.0);
}

void test_failure_during_integration()
{
    const Ensemble original{{1.0}, {2.0}};
    const qmomkit::AdaptiveRkOptions options = accurate_options();

    const qmomkit::PhaseSpaceVelocity nonfinite_callback = [](double, const std::vector<double>&)
    { return std::vector<double>{std::numeric_limits<double>::infinity()}; };
    const qmomkit::monte_carlo::EnsembleEvolutionResult nonfinite =
        qmomkit::monte_carlo::evolve_ensemble(original, nonfinite_callback, 0.0, 1.0, options);
    CHECK(!nonfinite.ok());
    CHECK(nonfinite.status == EnsembleEvolutionStatus::velocity_failed);
    CHECK(nonfinite.failed_sample_index == 0);
    CHECK(nonfinite.ensemble.empty());
    check_contains(nonfinite.message, "Phase-space velocity callback returned a non-finite value:");
    check_contains(nonfinite.message, "time: 0");
    check_contains(nonfinite.message, "dimension: 0");
    check_contains(nonfinite.message, "value: inf");
    check_contains(nonfinite.message, "state: [1]");
    check_contains(nonfinite.message, "velocity: [inf]");

    const qmomkit::PhaseSpaceVelocity wrong_dimension = [](double, const std::vector<double>&) { return std::vector<double>{1.0, 2.0}; };
    const qmomkit::monte_carlo::EnsembleEvolutionResult dimension =
        qmomkit::monte_carlo::evolve_ensemble(original, wrong_dimension, 0.0, 1.0, options);
    CHECK(!dimension.ok());
    CHECK(dimension.status == EnsembleEvolutionStatus::velocity_failed);
    CHECK(dimension.failed_sample_index == 0);
    check_contains(dimension.message, "Phase-space velocity callback returned the wrong dimension:");
    check_contains(dimension.message, "expected size: 1");
    check_contains(dimension.message, "returned size: 2");

    const qmomkit::PhaseSpaceVelocity throwing_callback = [](double, const std::vector<double>&) -> std::vector<double>
    { throw std::runtime_error("manufactured callback failure"); };
    const qmomkit::monte_carlo::EnsembleEvolutionResult thrown =
        qmomkit::monte_carlo::evolve_ensemble(original, throwing_callback, 0.0, 1.0, options);
    CHECK(!thrown.ok());
    CHECK(thrown.status == EnsembleEvolutionStatus::velocity_failed);
    CHECK(thrown.failed_sample_index == 0);

    qmomkit::AdaptiveRkOptions failing_options = accurate_options();
    failing_options.max_steps = 1;
    failing_options.max_step = 0.01;
    const qmomkit::PhaseSpaceVelocity linear = [](double, const std::vector<double>& coordinate) { return std::vector<double>{coordinate[0]}; };
    const qmomkit::monte_carlo::EnsembleEvolutionResult rk_failure =
        qmomkit::monte_carlo::evolve_ensemble(original, linear, 0.0, 1.0, failing_options);
    CHECK(!rk_failure.ok());
    CHECK(rk_failure.status == EnsembleEvolutionStatus::rk23_failed);
    CHECK(rk_failure.rk23_status == qmomkit::AdaptiveRkStatus::max_steps_exceeded);
    CHECK(rk_failure.failed_sample_index == 0);
    CHECK(rk_failure.ensemble.empty());

    const qmomkit::PhaseSpaceVelocity second_sample_failure = [](double, const std::vector<double>& coordinate)
    {
        if (coordinate[0] > 1.5)
        {
            throw std::runtime_error("manufactured second-sample failure");
        }
        return std::vector<double>{0.1};
    };
    const qmomkit::monte_carlo::EnsembleEvolutionResult later_failure =
        qmomkit::monte_carlo::evolve_ensemble(original, second_sample_failure, 0.0, 1.0, options);
    CHECK(!later_failure.ok());
    CHECK(later_failure.status == EnsembleEvolutionStatus::velocity_failed);
    CHECK(later_failure.failed_sample_index == 1);
    CHECK(later_failure.ensemble.empty());
    CHECK(later_failure.total_accepted_steps > 0);
    CHECK(later_failure.total_function_evaluations > 0);

    CHECK(original[0][0] == 1.0);
    CHECK(original[1][0] == 2.0);
}

void test_reproducibility()
{
    const qmomkit::DistributionSpecification distribution = qmomkit::NormalDistribution{0.2, 0.5};
    const SamplingOptions sampling_options{500, 12001};
    const qmomkit::monte_carlo::SamplingResult first_sampling = qmomkit::monte_carlo::sample_distribution(distribution, 1.8, sampling_options);
    const qmomkit::monte_carlo::SamplingResult second_sampling = qmomkit::monte_carlo::sample_distribution(distribution, 1.8, sampling_options);
    CHECK(first_sampling.ok());
    CHECK(second_sampling.ok());

    qmomkit::AdaptiveRkOptions options = accurate_options();
    options.max_step = 0.02;
    const qmomkit::PhaseSpaceVelocity callback = [](double time, const std::vector<double>& coordinate)
    { return std::vector<double>{0.1 + 0.2 * coordinate[0] + 0.03 * time}; };

    const qmomkit::monte_carlo::EnsembleEvolutionResult first =
        qmomkit::monte_carlo::evolve_ensemble(first_sampling.ensemble, callback, -0.2, 0.3, options);
    const qmomkit::monte_carlo::EnsembleEvolutionResult second =
        qmomkit::monte_carlo::evolve_ensemble(second_sampling.ensemble, callback, -0.2, 0.3, options);
    CHECK(first.ok());
    CHECK(second.ok());
    check_ensemble_equal(first.ensemble, second.ensemble);
    CHECK(first.total_accepted_steps == second.total_accepted_steps);
    CHECK(first.total_rejected_steps == second.total_rejected_steps);
    CHECK(first.total_function_evaluations == second.total_function_evaluations);
}

} // namespace

int main()
{
    test_single_sample_constant_velocity();
    test_convenience_evolution_uses_default_adaptive_initial_step();
    test_single_sample_linear_growth();
    test_deterministic_ensemble_rk23_max_step_convergence();
    test_single_sample_affine_evolution();
    test_ensemble_evolution_and_moments();
    test_failure_inputs_preserve_original();
    test_failure_during_integration();
    test_reproducibility();
    return 0;
}
