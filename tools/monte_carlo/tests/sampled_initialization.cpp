#include "monte_carlo/sampled_initialization.hpp"

#include "qmomkit/distributions.hpp"
#include "qmomkit/multi_index.hpp"
#include "test_case_dispatch.hpp"
#include "test_checks.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

namespace
{

double global_maximum_error = 0.0;

void update_maximum_error(const double error)
{
    global_maximum_error = std::max(global_maximum_error, error);
}

void check_close(const double actual, const double expected, const double tolerance)
{
    const double error = std::abs(actual - expected);
    update_maximum_error(error);
    CHECK_CLOSE(actual, expected, tolerance);
}

void test_one_dimensional_distribution_density()
{
    const double number_density = 2.75;
    const qmomkit::monte_carlo::OneDimensionalSampledInitialization initialization = qmomkit::monte_carlo::sample_initialization(
        qmomkit::NormalDistribution{0.5, 1.2}, number_density, 3, qmomkit::monte_carlo::SamplingOptions{2000, 1201});

    CHECK(initialization.ok());
    CHECK(initialization.sampling.ensemble.size() == 2000);
    CHECK(initialization.sampling.number_density == number_density);
    CHECK(initialization.sampling.seed == 1201);
    CHECK(initialization.moments.moments.size() == 4);
    check_close(initialization.moments.moments[0], number_density, 1.0e-14);
}

void test_mixture_absolute_density()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::NormalDistribution{-1.0, 0.4}, 1.25}, {qmomkit::GammaDistribution{2.0, 0.3}, 2.5}}};
    const qmomkit::monte_carlo::OneDimensionalSampledInitialization initialization =
        qmomkit::monte_carlo::sample_initialization(mixture, 2, qmomkit::monte_carlo::SamplingOptions{3000, 2201});

    CHECK(initialization.ok());
    check_close(initialization.sampling.number_density, 3.75, 0.0);
    check_close(initialization.moments.moments[0], 3.75, 1.0e-14);
}

void test_multidimensional_tensor_layout()
{
    const qmomkit::IndependentProductDistribution product{
        {qmomkit::GammaDistribution{2.0, 0.5}, qmomkit::NormalDistribution{0.25, 0.75}, qmomkit::BetaDistribution{2.0, 5.0}}, 4.0};
    const std::vector<std::size_t> maximum_orders{1, 2, 3};
    const qmomkit::monte_carlo::TensorSampledInitialization initialization =
        qmomkit::monte_carlo::sample_initialization(product, maximum_orders, qmomkit::monte_carlo::SamplingOptions{2500, 3201});

    CHECK(initialization.ok());
    CHECK((initialization.moments.moments.extents == std::vector<std::size_t>{2, 3, 4}));
    CHECK(initialization.moments.moments.values.size() == 24);
    check_close(initialization.moments.moments.values.front(), product.number_density, 1.0e-14);

    const qmomkit::MultiIndex order{1, 2, 3};
    const std::size_t flat = qmomkit::flatten_multi_index(order, initialization.moments.moments.extents);
    CHECK(flat == 23);
    CHECK(qmomkit::unflatten_multi_index(flat, initialization.moments.moments.extents) == order);
}

void check_close_relative(const double actual, const double expected, const double relative_tolerance, const std::string& context)
{
    update_maximum_error(std::abs(actual - expected));
    CHECK_CLOSE_CONTEXT(actual, expected, 0.0, relative_tolerance, context);
}

void test_selected_initialization_matches_dense_tensor()
{
    const qmomkit::IndependentProductDistribution product{
        {qmomkit::GammaDistribution{2.0, 0.5}, qmomkit::NormalDistribution{0.25, 0.75}, qmomkit::BetaDistribution{2.0, 5.0}}, 4.0};
    const qmomkit::monte_carlo::SamplingOptions options{2500, 3301};
    const std::vector<std::size_t> maximum_orders{1, 2, 3};
    const std::vector<qmomkit::MultiIndex> requested{{0, 0, 0}, {1, 0, 0}, {0, 2, 0}, {1, 2, 3}, {0, 1, 1}};

    const qmomkit::monte_carlo::TensorSampledInitialization dense = qmomkit::monte_carlo::sample_initialization(product, maximum_orders, options);
    const qmomkit::monte_carlo::SelectedSampledInitialization selected =
        qmomkit::monte_carlo::sample_initialization(product, qmomkit::RawMomentRequest{requested}, options);

    CHECK(dense.ok());
    CHECK(selected.ok());

    // The two overloads differ only in which moments they evaluate. Identical
    // sampling options draw the identical ensemble, so every shared index must
    // carry the same estimate; that equality is what lets a caller choose the
    // requested form purely to avoid evaluating the rest of the rectangle.
    CHECK(selected.sampling.ensemble.size() == dense.sampling.ensemble.size());
    CHECK(selected.sampling.seed == dense.sampling.seed);
    CHECK(selected.sampling.number_density == dense.sampling.number_density);
    CHECK(selected.moments.indices == requested);
    CHECK(selected.moments.values.size() == requested.size());
    CHECK(selected.moments.standard_errors.size() == requested.size());

    for (std::size_t entry = 0; entry < requested.size(); ++entry)
    {
        const std::size_t flat = qmomkit::flatten_multi_index(requested[entry], dense.moments.moments.extents);
        check_close_relative(selected.moments.values[entry], dense.moments.moments.values[flat], 1.0e-14, "requested entry " + std::to_string(entry));
        check_close_relative(selected.moments.standard_errors[entry], dense.moments.standard_errors[flat], 1.0e-14,
                             "requested standard error " + std::to_string(entry));
    }
}

void test_selected_initialization_error_propagation()
{
    const qmomkit::IndependentProductDistribution invalid_product{{qmomkit::NormalDistribution{0.0, -1.0}}, 1.0};
    const qmomkit::monte_carlo::SelectedSampledInitialization invalid = qmomkit::monte_carlo::sample_initialization(
        invalid_product, qmomkit::RawMomentRequest{{{0}, {1}}}, qmomkit::monte_carlo::SamplingOptions{100, 6301});

    CHECK(!invalid.ok());
    CHECK(!invalid.sampling.ok());
    CHECK(!invalid.moments.ok());
    // A failed draw must not leave a partially populated estimate that a caller
    // could mistake for a valid zero-moment realization.
    CHECK(invalid.moments.indices.empty());
    CHECK(invalid.moments.values.empty());
    CHECK(invalid.moments.standard_errors.empty());
}

void test_reproducibility()
{
    const qmomkit::DistributionSpecification distribution = qmomkit::LogNormalDistribution{-0.1, 0.35};
    const qmomkit::monte_carlo::SamplingOptions options{512, 4201};

    const qmomkit::monte_carlo::OneDimensionalSampledInitialization first =
        qmomkit::monte_carlo::sample_initialization(distribution, 1.5, 3, options);
    const qmomkit::monte_carlo::OneDimensionalSampledInitialization second =
        qmomkit::monte_carlo::sample_initialization(distribution, 1.5, 3, options);

    CHECK(first.ok());
    CHECK(second.ok());
    CHECK(first.sampling.ensemble == second.sampling.ensemble);
    CHECK(first.moments.moments.values() == second.moments.moments.values());
    CHECK(first.moments.standard_errors == second.moments.standard_errors);
}

void test_statistical_sanity_improves_with_sample_count()
{
    const qmomkit::DistributionSpecification distribution = qmomkit::NormalDistribution{0.3, 1.1};
    const double number_density = 2.0;
    const std::size_t maximum_order = 2;

    const qmomkit::RawMomentResult analytical = qmomkit::compute_raw_moments(distribution, maximum_order, number_density);
    CHECK(analytical.ok());

    const qmomkit::monte_carlo::OneDimensionalSampledInitialization small =
        qmomkit::monte_carlo::sample_initialization(distribution, number_density, maximum_order, qmomkit::monte_carlo::SamplingOptions{4000, 5201});
    const qmomkit::monte_carlo::OneDimensionalSampledInitialization large =
        qmomkit::monte_carlo::sample_initialization(distribution, number_density, maximum_order, qmomkit::monte_carlo::SamplingOptions{64000, 5201});

    CHECK(small.ok());
    CHECK(large.ok());

    const double small_error = std::abs(small.moments.moments[1] - analytical.moments[1]);
    const double large_error = std::abs(large.moments.moments[1] - analytical.moments[1]);
    update_maximum_error(small_error);
    update_maximum_error(large_error);
    CHECK(large_error < small_error);

    const double standard_error_ratio = large.moments.standard_errors[1] / small.moments.standard_errors[1];
    const double expected_ratio =
        std::sqrt(static_cast<double>(small.sampling.ensemble.size()) / static_cast<double>(large.sampling.ensemble.size()));
    CHECK(standard_error_ratio < 1.2 * expected_ratio);
    CHECK(standard_error_ratio > 0.8 * expected_ratio);
}

void test_error_propagation()
{
    const qmomkit::monte_carlo::OneDimensionalSampledInitialization invalid =
        qmomkit::monte_carlo::sample_initialization(qmomkit::NormalDistribution{0.0, -1.0}, 1.0, 2, qmomkit::monte_carlo::SamplingOptions{100, 6201});

    CHECK(!invalid.ok());
    CHECK(!invalid.sampling.ok());
    CHECK(!invalid.moments.ok());
}

} // namespace

int main(const int argc, char** argv)
{
    const int status = qmomkit::monte_carlo::test::dispatch_test_case(
        argc, argv,
        {
            {"one_dimensional_distribution_density", test_one_dimensional_distribution_density},
            {"mixture_absolute_density", test_mixture_absolute_density},
            {"multidimensional_tensor_layout", test_multidimensional_tensor_layout},
            {"selected_initialization_matches_dense_tensor", test_selected_initialization_matches_dense_tensor},
            {"selected_initialization_error_propagation", test_selected_initialization_error_propagation},
            {"reproducibility", test_reproducibility},
            {"statistical_sanity_improves_with_sample_count", test_statistical_sanity_improves_with_sample_count},
            {"error_propagation", test_error_propagation},
        });
    std::cout << "Monte Carlo sampled initialization maximum absolute error: " << global_maximum_error << '\n';
    return status;
}
