#include "monte_carlo/moments.hpp"
#include "monte_carlo/output.hpp"
#include "monte_carlo/sampling.hpp"
#include "qmomkit/distributions.hpp"
#include "qmomkit/multi_index.hpp"
#include "test_checks.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace
{

using qmomkit::monte_carlo::Ensemble;
using qmomkit::monte_carlo::SamplingOptions;

double global_maximum_error = 0.0;

void update_maximum_error(const double error)
{
    global_maximum_error = std::max(global_maximum_error, error);
}

void check_close_with_standard_error(const double estimate, const double expected, const double standard_error, const std::size_t order,
                                     const SamplingOptions& options)
{
    const double absolute_error = std::abs(estimate - expected);
    update_maximum_error(absolute_error);
    const double tolerance = std::max(2.0e-2 * std::max(1.0, std::abs(expected)), 8.0 * standard_error + 1.0e-12);
    std::ostringstream context;
    context << "order=" << order << " sample_count=" << options.sample_count << " seed=" << options.seed << " standard_error=" << standard_error;
    CHECK_CLOSE_CONTEXT(estimate, expected, tolerance, 0.0, context.str());
}

void compare_one_dimensional_to_analytical(const qmomkit::DistributionSpecification& distribution, const double number_density,
                                           const std::size_t maximum_order, const SamplingOptions& options)
{
    const qmomkit::monte_carlo::SamplingResult sampling = qmomkit::monte_carlo::sample_distribution(distribution, number_density, options);
    CHECK(sampling.ok());
    CHECK(sampling.number_density == number_density);

    const qmomkit::monte_carlo::OneDimensionalMomentEstimate estimate =
        qmomkit::monte_carlo::compute_ensemble_moments(sampling.ensemble, sampling.number_density, maximum_order);
    CHECK(estimate.ok());

    const qmomkit::RawMomentResult analytical = qmomkit::compute_raw_moments(distribution, maximum_order, number_density);
    CHECK(analytical.ok());
    CHECK(estimate.moments.size() == analytical.moments.size());
    for (std::size_t order = 0; order <= maximum_order; ++order)
    {
        check_close_with_standard_error(estimate.moments[order], analytical.moments[order], estimate.standard_errors[order], order, options);
    }
}

void test_single_distribution_validation()
{
    const std::vector<std::size_t> sample_counts{40000, 120000};
    for (const std::size_t sample_count : sample_counts)
    {
        compare_one_dimensional_to_analytical(qmomkit::NormalDistribution{0.75, 1.2}, 2.5, 4, SamplingOptions{sample_count, 4101});
        compare_one_dimensional_to_analytical(qmomkit::GammaDistribution{2.0, 0.7}, 1.75, 3, SamplingOptions{sample_count, 4102});
        compare_one_dimensional_to_analytical(qmomkit::LogNormalDistribution{-0.2, 0.35}, 3.0, 3, SamplingOptions{sample_count, 4103});
    }
}

void test_all_supported_single_family_sampling_smoke()
{
    const std::vector<qmomkit::DistributionSpecification> distributions{
        qmomkit::NormalDistribution{0.0, 1.0},  qmomkit::LogNormalDistribution{0.0, 0.3},
        qmomkit::BetaDistribution{2.0, 5.0},    qmomkit::GammaDistribution{3.0, 0.4},
        qmomkit::ChiDistribution{4.0},          qmomkit::ChiSquaredDistribution{5.0},
        qmomkit::WeibullDistribution{1.5, 0.8},
    };

    for (std::size_t index = 0; index < distributions.size(); ++index)
    {
        const qmomkit::monte_carlo::SamplingResult sampling =
            qmomkit::monte_carlo::sample_distribution(distributions[index], 1.25, SamplingOptions{2000, 5000 + index});
        CHECK(sampling.ok());
        CHECK(sampling.ensemble.size() == 2000);
        const qmomkit::monte_carlo::OneDimensionalMomentEstimate estimate =
            qmomkit::monte_carlo::compute_ensemble_moments(sampling.ensemble, sampling.number_density, 1);
        CHECK(estimate.ok());
        CHECK(std::abs(estimate.moments[0] - 1.25) < 1.0e-14);
    }
}

void compare_mixture_to_analytical(const qmomkit::DistributionMixture& mixture, const std::size_t maximum_order, const SamplingOptions& options)
{
    const qmomkit::monte_carlo::SamplingResult sampling = qmomkit::monte_carlo::sample_mixture(mixture, options);
    CHECK(sampling.ok());

    const qmomkit::monte_carlo::OneDimensionalMomentEstimate estimate =
        qmomkit::monte_carlo::compute_ensemble_moments(sampling.ensemble, sampling.number_density, maximum_order);
    CHECK(estimate.ok());

    const qmomkit::RawMomentResult analytical = qmomkit::compute_raw_moments(mixture, maximum_order);
    CHECK(analytical.ok());
    CHECK(std::abs(estimate.moments[0] - analytical.moments[0]) < 1.0e-14);
    for (std::size_t order = 0; order <= maximum_order; ++order)
    {
        check_close_with_standard_error(estimate.moments[order], analytical.moments[order], estimate.standard_errors[order], order, options);
    }
}

void test_mixture_validation()
{
    compare_mixture_to_analytical(
        qmomkit::DistributionMixture{{{qmomkit::NormalDistribution{-1.0, 0.5}, 1.25}, {qmomkit::NormalDistribution{2.0, 0.8}, 2.75}}}, 4,
        SamplingOptions{150000, 6101});

    compare_mixture_to_analytical(
        qmomkit::DistributionMixture{{{qmomkit::GammaDistribution{2.5, 0.4}, 1.2}, {qmomkit::LogNormalDistribution{0.1, 0.25}, 2.3}}}, 3,
        SamplingOptions{150000, 6102});
}

void compare_product_to_analytical(const qmomkit::IndependentProductDistribution& product, const std::vector<std::size_t>& maximum_orders,
                                   const SamplingOptions& options)
{
    const qmomkit::monte_carlo::SamplingResult sampling = qmomkit::monte_carlo::sample_independent_product(product, options);
    CHECK(sampling.ok());

    const qmomkit::monte_carlo::TensorMomentEstimate estimate =
        qmomkit::monte_carlo::compute_ensemble_moments(sampling.ensemble, sampling.number_density, maximum_orders);
    CHECK(estimate.ok());

    const qmomkit::TensorMomentResult analytical = qmomkit::compute_raw_moments(product, maximum_orders);
    CHECK(analytical.ok());
    CHECK(estimate.moments.extents == analytical.moments.extents);
    CHECK(estimate.moments.values.size() == analytical.moments.values.size());

    for (std::size_t flat = 0; flat < estimate.moments.values.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, estimate.moments.extents);
        CHECK(qmomkit::flatten_multi_index(index, estimate.moments.extents) == flat);
        check_close_with_standard_error(estimate.moments.values[flat], analytical.moments.values[flat], estimate.standard_errors[flat], flat,
                                        options);
    }
    CHECK(std::abs(estimate.moments.values.front() - product.number_density) < 1.0e-14);
}

void test_independent_product_validation()
{
    compare_product_to_analytical(
        qmomkit::IndependentProductDistribution{{qmomkit::GammaDistribution{2.0, 0.5}, qmomkit::NormalDistribution{0.5, 0.7}}, 2.2},
        std::vector<std::size_t>{2, 2}, SamplingOptions{160000, 7101});

    compare_product_to_analytical(
        qmomkit::IndependentProductDistribution{
            {qmomkit::LogNormalDistribution{-0.1, 0.25}, qmomkit::GammaDistribution{3.0, 0.2}, qmomkit::BetaDistribution{2.0, 4.0}}, 1.8},
        std::vector<std::size_t>{1, 2, 2}, SamplingOptions{180000, 7102});
}

void test_standard_error_convergence()
{
    const qmomkit::DistributionSpecification distribution = qmomkit::NormalDistribution{0.2, 1.1};
    const qmomkit::monte_carlo::SamplingResult small = qmomkit::monte_carlo::sample_distribution(distribution, 2.0, SamplingOptions{8000, 8101});
    const qmomkit::monte_carlo::SamplingResult large = qmomkit::monte_carlo::sample_distribution(distribution, 2.0, SamplingOptions{128000, 8101});
    CHECK(small.ok());
    CHECK(large.ok());

    const qmomkit::monte_carlo::OneDimensionalMomentEstimate small_moments =
        qmomkit::monte_carlo::compute_ensemble_moments(small.ensemble, small.number_density, 1);
    const qmomkit::monte_carlo::OneDimensionalMomentEstimate large_moments =
        qmomkit::monte_carlo::compute_ensemble_moments(large.ensemble, large.number_density, 1);
    CHECK(small_moments.ok());
    CHECK(large_moments.ok());

    const double expected_ratio = std::sqrt(static_cast<double>(small.ensemble.size()) / static_cast<double>(large.ensemble.size()));
    const double observed_ratio = large_moments.standard_errors[1] / small_moments.standard_errors[1];
    CHECK(observed_ratio < 1.2 * expected_ratio);
    CHECK(observed_ratio > 0.8 * expected_ratio);
}

void test_asymmetric_tensor_moment_layout()
{
    const Ensemble ensemble{{2.0, 3.0, 5.0}, {7.0, 11.0, 13.0}};
    const double number_density = 4.0;
    const std::vector<std::size_t> maximum_orders{1, 2, 3};
    const qmomkit::monte_carlo::TensorMomentEstimate estimate =
        qmomkit::monte_carlo::compute_ensemble_moments(ensemble, number_density, maximum_orders);
    CHECK(estimate.ok());
    CHECK((estimate.moments.extents == std::vector<std::size_t>{2, 3, 4}));

    const qmomkit::MultiIndex asymmetric_order{1, 2, 3};
    const std::size_t flat = qmomkit::flatten_multi_index(asymmetric_order, estimate.moments.extents);
    CHECK(flat == 23);

    const double expected = number_density * 0.5 * (2.0 * 3.0 * 3.0 * 5.0 * 5.0 * 5.0 + 7.0 * 11.0 * 11.0 * 13.0 * 13.0 * 13.0);
    CHECK(std::abs(estimate.moments.values[flat] - expected) < 1.0e-10 * std::max(1.0, std::abs(expected)));
    CHECK(std::abs(estimate.moments.values.front() - number_density) < 1.0e-14);
}

void test_invalid_inputs_and_output()
{
    const qmomkit::monte_carlo::SamplingResult invalid =
        qmomkit::monte_carlo::sample_distribution(qmomkit::NormalDistribution{0.0, -1.0}, 1.0, SamplingOptions{10, 1});
    CHECK(!invalid.ok());

    const qmomkit::monte_carlo::OneDimensionalMomentEstimate empty = qmomkit::monte_carlo::compute_ensemble_moments(Ensemble{}, 1.0, 2);
    CHECK(!empty.ok());

    const Ensemble one_dimensional{{1.0}};
    const qmomkit::monte_carlo::OneDimensionalMomentEstimate huge_order =
        qmomkit::monte_carlo::compute_ensemble_moments(one_dimensional, 1.0, std::numeric_limits<std::size_t>::max());
    CHECK(!huge_order.ok());

    const Ensemble two_dimensional{{1.0, 2.0}};
    const qmomkit::monte_carlo::TensorMomentEstimate huge_tensor_order =
        qmomkit::monte_carlo::compute_ensemble_moments(two_dimensional, 1.0, std::vector<std::size_t>{1, std::numeric_limits<std::size_t>::max()});
    CHECK(!huge_tensor_order.ok());

    const qmomkit::monte_carlo::TensorMomentEstimate overflowing_tensor = qmomkit::monte_carlo::compute_ensemble_moments(
        two_dimensional, 1.0, std::vector<std::size_t>{std::numeric_limits<std::size_t>::max() / 2, std::numeric_limits<std::size_t>::max() / 2});
    CHECK(!overflowing_tensor.ok());

    const qmomkit::monte_carlo::SamplingResult sampling =
        qmomkit::monte_carlo::sample_distribution(qmomkit::NormalDistribution{0.0, 1.0}, 1.5, SamplingOptions{500, 9101});
    CHECK(sampling.ok());
    const qmomkit::monte_carlo::OneDimensionalMomentEstimate estimate =
        qmomkit::monte_carlo::compute_ensemble_moments(sampling.ensemble, sampling.number_density, 2);
    CHECK(estimate.ok());

    std::ostringstream output;
    qmomkit::monte_carlo::write_moment_text(output, estimate, sampling.ensemble.size(), sampling.seed, sampling.number_density);
    const std::string text = output.str();
    CHECK(text.find("# samples=500") != std::string::npos);
    CHECK(text.find("# seed=9101") != std::string::npos);
    CHECK(text.find("# alpha_0 value standard_error") != std::string::npos);
    CHECK(text.find(',') == std::string::npos);
}

} // namespace

int main()
{
    test_all_supported_single_family_sampling_smoke();
    test_single_distribution_validation();
    test_mixture_validation();
    test_independent_product_validation();
    test_standard_error_convergence();
    test_asymmetric_tensor_moment_layout();
    test_invalid_inputs_and_output();

    std::cout << "Monte Carlo sampling moment validation maximum analytical absolute error: " << global_maximum_error << '\n';
    return 0;
}
