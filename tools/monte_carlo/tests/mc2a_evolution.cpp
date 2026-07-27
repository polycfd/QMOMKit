#include "qmomkit/distributions.hpp"
#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "monte_carlo/sampling.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{

using qmomkit::monte_carlo::Ensemble;
using qmomkit::monte_carlo::EnsembleEvolutionStatus;
using qmomkit::monte_carlo::Sample;
using qmomkit::monte_carlo::SamplingOptions;

qmomkit::AdaptiveRk23Options accurate_options()
{
    qmomkit::AdaptiveRk23Options options;
    options.relative_tolerance = 1.0e-10;
    options.absolute_tolerance = 1.0e-12;
    options.max_step = 0.01;
    return options;
}

void assert_close(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

void assert_ensemble_equal(const Ensemble& lhs, const Ensemble& rhs)
{
    assert(lhs.size() == rhs.size());
    for (std::size_t sample = 0; sample < lhs.size(); ++sample)
    {
        assert(lhs[sample] == rhs[sample]);
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
        assert(coordinate.size() == 1);
        return std::vector<double>{velocity};
    };

    const qmomkit::monte_carlo::EnsembleEvolutionResult result =
        qmomkit::monte_carlo::evolve_ensemble(ensemble, callback, t0, t1, accurate_options());

    assert(result.ok());
    assert(result.ensemble.size() == 1);
    assert(ensemble[0][0] == 1.25);
    assert_close(result.ensemble[0][0], 1.25 + velocity * (t1 - t0), 1.0e-12);
    assert_close(result.final_time, t1, 0.0);
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

    assert(result.ok());
    assert_close(result.ensemble[0][0], ensemble[0][0] * std::exp(rate * (t1 - t0)), 1.0e-9);
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
    assert(result.ok());
    assert_close(result.ensemble[0][0], expected, 1.0e-9);
}

void test_ensemble_evolution_and_moments()
{
    const qmomkit::IndependentProductDistribution distribution{{qmomkit::NormalDistribution{0.25, 0.3}, qmomkit::GammaDistribution{2.0, 0.4}},
                                                                 2.5};
    const qmomkit::monte_carlo::SamplingResult sampling =
        qmomkit::monte_carlo::sample_independent_product(distribution, SamplingOptions{2500, 11001});
    assert(sampling.ok());

    const double t0 = 0.1;
    const double t1 = 0.4;
    const qmomkit::PhaseSpaceVelocity callback = [](double, const std::vector<double>& coordinate)
    {
        assert(coordinate.size() == 2);
        return std::vector<double>{0.2 + 0.1 * coordinate[0], -0.15 * coordinate[1]};
    };

    const qmomkit::monte_carlo::EnsembleEvolutionResult result =
        qmomkit::monte_carlo::evolve_ensemble(sampling.ensemble, callback, t0, t1, accurate_options());
    assert(result.ok());
    assert(result.ensemble.size() == sampling.ensemble.size());
    assert(result.total_accepted_steps > 0);
    assert(result.total_function_evaluations > 0);
    assert(result.maximum_rejected_steps <= result.total_rejected_steps);
    assert_close(result.final_time, t1, 0.0);

    for (const Sample& sample : result.ensemble)
    {
        assert(sample.size() == 2);
        assert(std::isfinite(sample[0]));
        assert(std::isfinite(sample[1]));
    }

    const qmomkit::monte_carlo::TensorMomentEstimate moments =
        qmomkit::monte_carlo::compute_ensemble_moments(result.ensemble, sampling.number_density, std::vector<std::size_t>{2, 2});
    assert(moments.ok());
    assert_close(moments.moments.values.front(), sampling.number_density, 0.0);
    for (const double value : moments.moments.values)
    {
        assert(std::isfinite(value));
    }
    for (const double error : moments.standard_errors)
    {
        assert(std::isfinite(error));
    }
}

void test_failure_inputs_preserve_original()
{
    const Ensemble original{{1.0}, {2.0}};
    const qmomkit::PhaseSpaceVelocity identity = [](double, const std::vector<double>& coordinate) { return coordinate; };
    const qmomkit::AdaptiveRk23Options options = accurate_options();

    const qmomkit::monte_carlo::EnsembleEvolutionResult empty = qmomkit::monte_carlo::evolve_ensemble(Ensemble{}, identity, 0.0, 1.0, options);
    assert(!empty.ok());
    assert(empty.status == EnsembleEvolutionStatus::invalid_input);
    assert(empty.ensemble.empty());

    const qmomkit::monte_carlo::EnsembleEvolutionResult inconsistent =
        qmomkit::monte_carlo::evolve_ensemble(Ensemble{{1.0}, {1.0, 2.0}}, identity, 0.0, 1.0, options);
    assert(!inconsistent.ok());
    assert(inconsistent.failed_sample_index == 1);
    assert(inconsistent.ensemble.empty());

    const qmomkit::monte_carlo::EnsembleEvolutionResult nonfinite =
        qmomkit::monte_carlo::evolve_ensemble(Ensemble{{std::numeric_limits<double>::quiet_NaN()}}, identity, 0.0, 1.0, options);
    assert(!nonfinite.ok());
    assert(nonfinite.failed_sample_index == 0);

    const qmomkit::monte_carlo::EnsembleEvolutionResult invalid_time =
        qmomkit::monte_carlo::evolve_ensemble(original, identity, 0.0, 0.0, options);
    assert(!invalid_time.ok());
    assert(invalid_time.status == EnsembleEvolutionStatus::invalid_input);

    const qmomkit::PhaseSpaceVelocity missing;
    const qmomkit::monte_carlo::EnsembleEvolutionResult missing_callback =
        qmomkit::monte_carlo::evolve_ensemble(original, missing, 0.0, 1.0, options);
    assert(!missing_callback.ok());
    assert(missing_callback.status == EnsembleEvolutionStatus::invalid_input);

    assert(original[0][0] == 1.0);
    assert(original[1][0] == 2.0);
}

void test_failure_during_integration()
{
    const Ensemble original{{1.0}, {2.0}};
    const qmomkit::AdaptiveRk23Options options = accurate_options();

    const qmomkit::PhaseSpaceVelocity nonfinite_callback = [](double, const std::vector<double>&)
    { return std::vector<double>{std::numeric_limits<double>::infinity()}; };
    const qmomkit::monte_carlo::EnsembleEvolutionResult nonfinite =
        qmomkit::monte_carlo::evolve_ensemble(original, nonfinite_callback, 0.0, 1.0, options);
    assert(!nonfinite.ok());
    assert(nonfinite.status == EnsembleEvolutionStatus::velocity_failed);
    assert(nonfinite.failed_sample_index == 0);
    assert(nonfinite.ensemble.empty());

    const qmomkit::PhaseSpaceVelocity wrong_dimension = [](double, const std::vector<double>&) { return std::vector<double>{1.0, 2.0}; };
    const qmomkit::monte_carlo::EnsembleEvolutionResult dimension =
        qmomkit::monte_carlo::evolve_ensemble(original, wrong_dimension, 0.0, 1.0, options);
    assert(!dimension.ok());
    assert(dimension.status == EnsembleEvolutionStatus::velocity_failed);
    assert(dimension.failed_sample_index == 0);

    const qmomkit::PhaseSpaceVelocity throwing_callback = [](double, const std::vector<double>&) -> std::vector<double>
    { throw std::runtime_error("manufactured callback failure"); };
    const qmomkit::monte_carlo::EnsembleEvolutionResult thrown =
        qmomkit::monte_carlo::evolve_ensemble(original, throwing_callback, 0.0, 1.0, options);
    assert(!thrown.ok());
    assert(thrown.status == EnsembleEvolutionStatus::velocity_failed);
    assert(thrown.failed_sample_index == 0);

    qmomkit::AdaptiveRk23Options failing_options = accurate_options();
    failing_options.max_steps = 1;
    failing_options.max_step = 0.01;
    const qmomkit::PhaseSpaceVelocity linear = [](double, const std::vector<double>& coordinate) { return std::vector<double>{coordinate[0]}; };
    const qmomkit::monte_carlo::EnsembleEvolutionResult rk_failure =
        qmomkit::monte_carlo::evolve_ensemble(original, linear, 0.0, 1.0, failing_options);
    assert(!rk_failure.ok());
    assert(rk_failure.status == EnsembleEvolutionStatus::rk23_failed);
    assert(rk_failure.rk23_status == qmomkit::AdaptiveRk23Status::max_steps_exceeded);
    assert(rk_failure.failed_sample_index == 0);
    assert(rk_failure.ensemble.empty());

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
    assert(!later_failure.ok());
    assert(later_failure.status == EnsembleEvolutionStatus::velocity_failed);
    assert(later_failure.failed_sample_index == 1);
    assert(later_failure.ensemble.empty());
    assert(later_failure.total_accepted_steps > 0);
    assert(later_failure.total_function_evaluations > 0);

    assert(original[0][0] == 1.0);
    assert(original[1][0] == 2.0);
}

void test_reproducibility()
{
    const qmomkit::DistributionSpecification distribution = qmomkit::NormalDistribution{0.2, 0.5};
    const SamplingOptions sampling_options{500, 12001};
    const qmomkit::monte_carlo::SamplingResult first_sampling = qmomkit::monte_carlo::sample_distribution(distribution, 1.8, sampling_options);
    const qmomkit::monte_carlo::SamplingResult second_sampling = qmomkit::monte_carlo::sample_distribution(distribution, 1.8, sampling_options);
    assert(first_sampling.ok());
    assert(second_sampling.ok());

    qmomkit::AdaptiveRk23Options options = accurate_options();
    options.max_step = 0.02;
    const qmomkit::PhaseSpaceVelocity callback = [](double time, const std::vector<double>& coordinate)
    { return std::vector<double>{0.1 + 0.2 * coordinate[0] + 0.03 * time}; };

    const qmomkit::monte_carlo::EnsembleEvolutionResult first =
        qmomkit::monte_carlo::evolve_ensemble(first_sampling.ensemble, callback, -0.2, 0.3, options);
    const qmomkit::monte_carlo::EnsembleEvolutionResult second =
        qmomkit::monte_carlo::evolve_ensemble(second_sampling.ensemble, callback, -0.2, 0.3, options);
    assert(first.ok());
    assert(second.ok());
    assert_ensemble_equal(first.ensemble, second.ensemble);
    assert(first.total_accepted_steps == second.total_accepted_steps);
    assert(first.total_rejected_steps == second.total_rejected_steps);
    assert(first.total_function_evaluations == second.total_function_evaluations);
}

} // namespace

int main()
{
    test_single_sample_constant_velocity();
    test_single_sample_linear_growth();
    test_single_sample_affine_evolution();
    test_ensemble_evolution_and_moments();
    test_failure_inputs_preserve_original();
    test_failure_during_integration();
    test_reproducibility();
    return 0;
}
