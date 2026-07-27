#include "monte_carlo/sampled_initialization.hpp"

#include "qmomkit/distributions.hpp"
#include "qmomkit/multi_index.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

namespace
{

double global_maximum_error = 0.0;

void update_maximum_error(const double error)
{
    global_maximum_error = std::max(global_maximum_error, error);
}

void assert_close(const double actual, const double expected, const double tolerance)
{
    const double error = std::abs(actual - expected);
    update_maximum_error(error);
    assert(error <= tolerance);
}

void test_one_dimensional_distribution_density()
{
    const double number_density = 2.75;
    const qmomkit::monte_carlo::OneDimensionalSampledInitialization initialization = qmomkit::monte_carlo::sample_initialization(
        qmomkit::NormalDistribution{0.5, 1.2}, number_density, 3, qmomkit::monte_carlo::SamplingOptions{2000, 1201});

    assert(initialization.ok());
    assert(initialization.sampling.ensemble.size() == 2000);
    assert(initialization.sampling.number_density == number_density);
    assert(initialization.sampling.seed == 1201);
    assert(initialization.moments.moments.size() == 4);
    assert_close(initialization.moments.moments[0], number_density, 1.0e-14);
}

void test_mixture_absolute_density()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::NormalDistribution{-1.0, 0.4}, 1.25}, {qmomkit::GammaDistribution{2.0, 0.3}, 2.5}}};
    const qmomkit::monte_carlo::OneDimensionalSampledInitialization initialization =
        qmomkit::monte_carlo::sample_initialization(mixture, 2, qmomkit::monte_carlo::SamplingOptions{3000, 2201});

    assert(initialization.ok());
    assert_close(initialization.sampling.number_density, 3.75, 0.0);
    assert_close(initialization.moments.moments[0], 3.75, 1.0e-14);
}

void test_multidimensional_tensor_layout()
{
    const qmomkit::IndependentProductDistribution product{
        {qmomkit::GammaDistribution{2.0, 0.5}, qmomkit::NormalDistribution{0.25, 0.75}, qmomkit::BetaDistribution{2.0, 5.0}}, 4.0};
    const std::vector<std::size_t> maximum_orders{1, 2, 3};
    const qmomkit::monte_carlo::TensorSampledInitialization initialization =
        qmomkit::monte_carlo::sample_initialization(product, maximum_orders, qmomkit::monte_carlo::SamplingOptions{2500, 3201});

    assert(initialization.ok());
    assert((initialization.moments.moments.extents == std::vector<std::size_t>{2, 3, 4}));
    assert(initialization.moments.moments.values.size() == 24);
    assert_close(initialization.moments.moments.values.front(), product.number_density, 1.0e-14);

    const qmomkit::MultiIndex order{1, 2, 3};
    const std::size_t flat = qmomkit::flatten_multi_index(order, initialization.moments.moments.extents);
    assert(flat == 23);
    assert(qmomkit::unflatten_multi_index(flat, initialization.moments.moments.extents) == order);
}

void test_reproducibility()
{
    const qmomkit::DistributionSpecification distribution = qmomkit::LogNormalDistribution{-0.1, 0.35};
    const qmomkit::monte_carlo::SamplingOptions options{512, 4201};

    const qmomkit::monte_carlo::OneDimensionalSampledInitialization first =
        qmomkit::monte_carlo::sample_initialization(distribution, 1.5, 3, options);
    const qmomkit::monte_carlo::OneDimensionalSampledInitialization second =
        qmomkit::monte_carlo::sample_initialization(distribution, 1.5, 3, options);

    assert(first.ok());
    assert(second.ok());
    assert(first.sampling.ensemble == second.sampling.ensemble);
    assert(first.moments.moments.values() == second.moments.moments.values());
    assert(first.moments.standard_errors == second.moments.standard_errors);
}

void test_statistical_sanity_improves_with_sample_count()
{
    const qmomkit::DistributionSpecification distribution = qmomkit::NormalDistribution{0.3, 1.1};
    const double number_density = 2.0;
    const std::size_t maximum_order = 2;

    const qmomkit::RawMomentResult analytical = qmomkit::compute_raw_moments(distribution, maximum_order, number_density);
    assert(analytical.ok());

    const qmomkit::monte_carlo::OneDimensionalSampledInitialization small = qmomkit::monte_carlo::sample_initialization(
        distribution, number_density, maximum_order, qmomkit::monte_carlo::SamplingOptions{4000, 5201});
    const qmomkit::monte_carlo::OneDimensionalSampledInitialization large = qmomkit::monte_carlo::sample_initialization(
        distribution, number_density, maximum_order, qmomkit::monte_carlo::SamplingOptions{64000, 5201});

    assert(small.ok());
    assert(large.ok());

    const double small_error = std::abs(small.moments.moments[1] - analytical.moments[1]);
    const double large_error = std::abs(large.moments.moments[1] - analytical.moments[1]);
    update_maximum_error(small_error);
    update_maximum_error(large_error);
    assert(large_error < small_error);

    const double standard_error_ratio = large.moments.standard_errors[1] / small.moments.standard_errors[1];
    const double expected_ratio =
        std::sqrt(static_cast<double>(small.sampling.ensemble.size()) / static_cast<double>(large.sampling.ensemble.size()));
    assert(standard_error_ratio < 1.2 * expected_ratio);
    assert(standard_error_ratio > 0.8 * expected_ratio);
}

void test_error_propagation()
{
    const qmomkit::monte_carlo::OneDimensionalSampledInitialization invalid = qmomkit::monte_carlo::sample_initialization(
        qmomkit::NormalDistribution{0.0, -1.0}, 1.0, 2, qmomkit::monte_carlo::SamplingOptions{100, 6201});

    assert(!invalid.ok());
    assert(!invalid.sampling.ok());
    assert(!invalid.moments.ok());
}

} // namespace

int main()
{
    test_one_dimensional_distribution_density();
    test_mixture_absolute_density();
    test_multidimensional_tensor_layout();
    test_reproducibility();
    test_statistical_sanity_improves_with_sample_count();
    test_error_propagation();

    std::cout << "Monte Carlo sampled initialization maximum absolute error: " << global_maximum_error << '\n';
    return 0;
}
