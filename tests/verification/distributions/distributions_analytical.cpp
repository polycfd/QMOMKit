// Raw moments computed from named distributions, against their textbook
// closed-form expressions.
//
// What is verified
//   That initializing moments from a named distribution - normal, log-normal,
//   beta, gamma, exponential, Weibull, chi, chi-squared, and mixtures and
//   independent products of them - produces the raw moments that distribution
//   is known to have, and that the moments of narrow gamma distributions and
//   of one configured independent product invert to quadratures that
//   reproduce them.
//
// How the reference is constructed
//   Not invented for this test. Each expected moment is evaluated from the
//   published expression for that family, written out in this file:
//
//     log-normal   E[X^k] = exp(k*mu + k^2*sigma^2/2)
//     beta         E[X^k] = product over j < k of (alpha + j)/(alpha + beta + j)
//     Weibull      E[X^k] = scale^k * Gamma(1 + k/shape)
//     chi-squared  E[X^k] = product over j < k of (nu + 2j)
//     chi          E[X^k] = 2^(k/2) Gamma((nu + k)/2) / Gamma(nu/2), with the
//                  gamma values at integers and half-integers built from
//                  Gamma(1) = 1 and Gamma(1/2) = sqrt(pi) by Gamma(x+1) = x Gamma(x)
//     gamma, normal, exponential   the corresponding standard expressions
//
//   Moments are scaled by the configured total number density, since QMOMKit
//   carries raw moments of the number density function rather than of a
//   normalized distribution, and m0 is the total number density rather than 1.
//
// Error metric and pass criteria
//   Each moment is compared relative to max(|expected|, 1), so a single bound
//   works across orders whose magnitudes differ by many decades - a high-order
//   log-normal moment is enormous while a low-order one is near unity, and an
//   absolute bound alone would be meaningless for one of them. Tolerances are
//   1.0e-14 to 1.0e-13 for single families, 1.0e-13 to 1.0e-12 for mixtures,
//   configured moments and product tensors, with 1.0e-11 absolute where a
//   tensor holds log-normal and Weibull factors, 1.0e-12 (1.0e-14 for the
//   first) for the odd moments of two normals that cancel to zero, and
//   5.0e-10 for moments rebuilt from a CQMOM reconstruction. Measured errors
//   stay at a few rounding errors, so the looser bounds are safety bounds.
//
//   The chi, chi-squared, gamma and beta cases with closed-form products, and
//   the selected-order cases, compare purely relatively, since their moments
//   span from below 1e-4 to above 1e15. The library forms these moments as
//   products of one simple factor per order, so they are accurate to a few
//   rounding errors whatever the parameters. Against the closed-form products,
//   themselves accurate to a few rounding errors, measured errors stay below
//   5e-16; the tolerance of 4e-15 leaves a margin of about eight.
//
//   For parameters up to 1e8 the references are the same products multiplied
//   out in error-free double-double arithmetic, accurate to about 1e-30, and a
//   moment of order k must lie within 4 (k + 1) epsilon; measured errors stay
//   below 2.2 epsilon. Odd chi moments are compared with values of E[X]
//   computed with mpmath at 40 digits, as recorded with that case.
//
//   Narrow gamma distributions must also invert: their moments are inverted
//   with Wheeler's algorithm and must be reproduced by the result to 1e-13,
//   where measured errors stay below 2e-15.
//
// Reproduce
//   ctest -R verification_distributions
//
// Build modes
//   Both. Failures are reported through an exit code, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/distributions.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/wheeler.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace
{

double maximum_relative_error = 0.0;
double maximum_absolute_error = 0.0;

void assert_near_relative(const double actual, const double expected, const double tolerance)
{
    maximum_absolute_error = std::max(maximum_absolute_error, std::abs(actual - expected));
    const double scale = std::max(1.0, std::abs(expected));
    const double error = std::abs(actual - expected) / scale;
    maximum_relative_error = std::max(maximum_relative_error, error);
    QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(actual, expected, 0.0, tolerance);
}

void assert_near_absolute_relative(const double actual, const double expected, const double absolute_tolerance, const double relative_tolerance)
{
    const double absolute_error = std::abs(actual - expected);
    maximum_absolute_error = std::max(maximum_absolute_error, absolute_error);
    const double scale = std::max(1.0, std::abs(expected));
    const double relative_error = absolute_error / scale;
    maximum_relative_error = std::max(maximum_relative_error, relative_error);
    QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(actual, expected, absolute_tolerance, relative_tolerance);
}

void assert_moments(const qmomkit::RawMomentResult& result, const std::vector<double>& expected, const double tolerance)
{
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.moments.size() == expected.size());
    for (std::size_t order = 0; order < expected.size(); ++order)
    {
        assert_near_relative(result.moments[order], expected[order], tolerance);
    }
}

void assert_moment_sets_match(const qmomkit::MomentSet& actual, const qmomkit::MomentSet& expected, const double tolerance)
{
    QMOMKIT_VERIFY_CHECK(actual.size() == expected.size());
    for (std::size_t order = 0; order < actual.size(); ++order)
    {
        assert_near_relative(actual[order], expected[order], tolerance);
    }
}

void assert_tensor(const qmomkit::TensorMomentResult& result, const std::vector<std::size_t>& expected_extents,
                   const std::vector<double>& expected_values, const double absolute_tolerance, const double relative_tolerance)
{
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.moments.extents == expected_extents);
    QMOMKIT_VERIFY_CHECK(result.moments.values.size() == expected_values.size());
    for (std::size_t flat_index = 0; flat_index < expected_values.size(); ++flat_index)
    {
        assert_near_absolute_relative(result.moments.values[flat_index], expected_values[flat_index], absolute_tolerance, relative_tolerance);
    }
}

void assert_tensor_values(const qmomkit::RawMomentTensor& tensor, const std::vector<std::size_t>& expected_extents,
                          const std::vector<double>& expected_values, const double absolute_tolerance, const double relative_tolerance)
{
    QMOMKIT_VERIFY_CHECK(tensor.extents == expected_extents);
    QMOMKIT_VERIFY_CHECK(tensor.values.size() == expected_values.size());
    for (std::size_t flat_index = 0; flat_index < expected_values.size(); ++flat_index)
    {
        assert_near_absolute_relative(tensor.values[flat_index], expected_values[flat_index], absolute_tolerance, relative_tolerance);
    }
}

double factorial(const std::size_t n)
{
    double value = 1.0;
    for (std::size_t factor = 2; factor <= n; ++factor)
    {
        value *= static_cast<double>(factor);
    }
    return value;
}

double double_factorial_odd(const std::size_t n)
{
    double value = 1.0;
    for (std::size_t factor = 1; factor <= n; factor += 2)
    {
        value *= static_cast<double>(factor);
    }
    return value;
}

double normal_raw_moment_sum(const double mean, const double standard_deviation, const std::size_t order)
{
    double value = 0.0;
    for (std::size_t j = 0; j <= order / 2; ++j)
    {
        const std::size_t remaining = order - 2 * j;
        const double coefficient = factorial(order) / (std::pow(2.0, static_cast<double>(j)) * factorial(j) * factorial(remaining));
        value += coefficient * std::pow(standard_deviation, static_cast<double>(2 * j)) * std::pow(mean, static_cast<double>(remaining));
    }
    return value;
}

double log_normal_raw_moment(const double normal_mean, const double normal_standard_deviation, const std::size_t order)
{
    const double k = static_cast<double>(order);
    return std::exp(k * normal_mean + 0.5 * k * k * normal_standard_deviation * normal_standard_deviation);
}

double beta_raw_moment(const double alpha, const double beta, const std::size_t order)
{
    double value = 1.0;
    for (std::size_t j = 0; j < order; ++j)
    {
        value *= (alpha + static_cast<double>(j)) / (alpha + beta + static_cast<double>(j));
    }
    return value;
}

double gamma_raw_moment(const double shape, const double scale, const std::size_t order)
{
    double value = 1.0;
    for (std::size_t j = 0; j < order; ++j)
    {
        value *= scale * (shape + static_cast<double>(j));
    }
    return value;
}

// Gamma(m / 2) for a positive integer m, from Gamma(1) = 1 and
// Gamma(1/2) = sqrt(pi) by the recurrence Gamma(x + 1) = x Gamma(x).
double gamma_of_half_integer(const std::size_t m)
{
    double value = m % 2 == 0 ? 1.0 : std::sqrt(std::acos(-1.0));
    for (std::size_t twice_x = m % 2 == 0 ? 2 : 1; twice_x < m; twice_x += 2)
    {
        value *= 0.5 * static_cast<double>(twice_x);
    }
    return value;
}

double chi_raw_moment(const std::size_t degrees_of_freedom, const std::size_t order)
{
    return std::pow(2.0, 0.5 * static_cast<double>(order)) * gamma_of_half_integer(degrees_of_freedom + order) /
           gamma_of_half_integer(degrees_of_freedom);
}

double chi_squared_raw_moment(const double degrees_of_freedom, const std::size_t order)
{
    double value = 1.0;
    for (std::size_t j = 0; j < order; ++j)
    {
        value *= degrees_of_freedom + 2.0 * static_cast<double>(j);
    }
    return value;
}

// A product kept as the unevaluated sum hi + lo of two doubles. Each factor is
// multiplied in error-free arithmetic, since the exact product of two doubles
// is p + fma(hi, x, -p), so a product of a few dozen factors stays accurate to
// about 1e-30, far below the rounding of a double result.
struct ExactProduct
{
    double hi = 1.0;
    double lo = 0.0;

    void multiply(const double factor)
    {
        const double product = hi * factor;
        const double error = std::fma(hi, factor, -product);
        const double low = lo * factor + error;
        hi = product + low;
        lo = low - (hi - product);
    }

    [[nodiscard]] double value() const
    {
        return hi + lo;
    }
};

// numerator / denominator for two exact products, with one correction step.
double exact_quotient(const ExactProduct& numerator, const ExactProduct& denominator)
{
    const double quotient = numerator.hi / denominator.hi;
    const double product = quotient * denominator.hi;
    const double product_error = std::fma(quotient, denominator.hi, -product);
    const double remainder = ((numerator.hi - product) - product_error) + numerator.lo - quotient * denominator.lo;
    return quotient + remainder / denominator.hi;
}

// The library forms a moment of order k as a product of k simple factors, each
// adding a few rounding errors.
double product_tolerance(const std::size_t order)
{
    return 4.0 * static_cast<double>(order + 1) * std::numeric_limits<double>::epsilon();
}

// Pure relative comparison, for moments that span many decades and fall below
// one.
void assert_relative(const double actual, const double expected, const double tolerance)
{
    QMOMKIT_VERIFY_NEAR(actual, expected, tolerance * std::abs(expected));
}

template <class Distribution, class Reference>
void assert_moments_relative(const Distribution& distribution, const double total_weight, const std::size_t maximum_order, const Reference& reference,
                             const double tolerance)
{
    const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(distribution, maximum_order, total_weight);
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.moments.size() == maximum_order + 1);
    for (std::size_t order = 0; order <= maximum_order; ++order)
    {
        assert_relative(result.moments[order], total_weight * reference(order), tolerance);
    }
}

// The selected-order path, with the orders requested out of sequence.
template <class Distribution, class Reference>
void assert_selected_moments_relative(const Distribution& distribution, const double total_weight, const Reference& reference, const double tolerance)
{
    const std::vector<std::size_t> orders{10, 0, 3, 7};
    qmomkit::RawMomentRequest request;
    for (const std::size_t order : orders)
    {
        request.indices.push_back(qmomkit::MultiIndex{order});
    }
    const qmomkit::SelectedRawMomentResult result = qmomkit::compute_raw_moments(distribution, request, total_weight);
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.values.size() == orders.size());
    for (std::size_t index = 0; index < orders.size(); ++index)
    {
        assert_relative(result.values[index], total_weight * reference(orders[index]), tolerance);
    }
}

double weibull_raw_moment(const double shape, const double scale, const std::size_t order)
{
    const double k = static_cast<double>(order);
    return std::pow(scale, k) * std::tgamma(1.0 + k / shape);
}

void verify_standard_normal()
{
    std::vector<double> expected(9, 0.0);
    expected[0] = 1.0;
    for (std::size_t order = 2; order < expected.size(); order += 2)
    {
        expected[order] = double_factorial_odd(order - 1);
    }

    assert_moments(qmomkit::compute_raw_moments(qmomkit::NormalDistribution{0.0, 1.0}, 8), expected, 1.0e-14);
}

void verify_shifted_scaled_normal()
{
    constexpr double total_weight = 2.5;
    std::vector<double> expected(7, 0.0);
    for (std::size_t order = 0; order < expected.size(); ++order)
    {
        expected[order] = total_weight * normal_raw_moment_sum(1.25, 0.75, order);
    }

    assert_moments(qmomkit::compute_raw_moments(qmomkit::NormalDistribution{1.25, 0.75}, 6, total_weight), expected, 1.0e-13);
}

void verify_log_normal()
{
    constexpr double total_weight = 3.0;
    std::vector<double> expected(6, 0.0);
    for (std::size_t order = 0; order < expected.size(); ++order)
    {
        const double k = static_cast<double>(order);
        expected[order] = total_weight * std::exp(k * 0.25 + 0.5 * k * k * 0.16);
    }

    assert_moments(qmomkit::compute_raw_moments(qmomkit::LogNormalDistribution{0.25, 0.4}, 5, total_weight), expected, 1.0e-13);
}

void verify_beta()
{
    constexpr double total_weight = 4.0;
    std::vector<double> uniform(7, 0.0);
    for (std::size_t order = 0; order < uniform.size(); ++order)
    {
        uniform[order] = total_weight / static_cast<double>(order + 1);
    }
    assert_moments(qmomkit::compute_raw_moments(qmomkit::BetaDistribution{1.0, 1.0}, 6, total_weight), uniform, 1.0e-14);

    std::vector<double> nonsymmetric(6, 0.0);
    double expectation = 1.0;
    for (std::size_t order = 0; order < nonsymmetric.size(); ++order)
    {
        if (order > 0)
        {
            expectation *= (2.0 + static_cast<double>(order - 1)) / (7.0 + static_cast<double>(order - 1));
        }
        nonsymmetric[order] = total_weight * expectation;
    }
    assert_moments(qmomkit::compute_raw_moments(qmomkit::BetaDistribution{2.0, 5.0}, 5, total_weight), nonsymmetric, 1.0e-14);
}

void verify_gamma()
{
    constexpr double total_weight = 2.0;
    std::vector<double> exponential(7, 0.0);
    for (std::size_t order = 0; order < exponential.size(); ++order)
    {
        exponential[order] = total_weight * std::pow(3.0, static_cast<double>(order)) * factorial(order);
    }
    assert_moments(qmomkit::compute_raw_moments(qmomkit::GammaDistribution{1.0, 3.0}, 6, total_weight), exponential, 1.0e-13);

    std::vector<double> half_shape{1.0, 1.0, 3.0, 15.0, 105.0};
    assert_moments(qmomkit::compute_raw_moments(qmomkit::GammaDistribution{0.5, 2.0}, 4), half_shape, 1.0e-13);
}

void verify_chi_squared()
{
    constexpr double total_weight = 1.7;
    const qmomkit::RawMomentResult chi_squared = qmomkit::compute_raw_moments(qmomkit::ChiSquaredDistribution{5.0}, 6, total_weight);
    const qmomkit::RawMomentResult gamma = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{2.5, 2.0}, 6, total_weight);

    QMOMKIT_VERIFY_CHECK(chi_squared.ok());
    QMOMKIT_VERIFY_CHECK(gamma.ok());
    QMOMKIT_VERIFY_CHECK(chi_squared.moments.size() == gamma.moments.size());
    for (std::size_t order = 0; order < chi_squared.moments.size(); ++order)
    {
        assert_near_relative(chi_squared.moments[order], gamma.moments[order], 1.0e-13);
    }
}

void verify_weibull_gamma_identity()
{
    constexpr double total_weight = 2.75;
    const qmomkit::RawMomentResult weibull = qmomkit::compute_raw_moments(qmomkit::WeibullDistribution{1.0, 1.5}, 6, total_weight);
    const qmomkit::RawMomentResult gamma = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{1.0, 1.5}, 6, total_weight);

    QMOMKIT_VERIFY_CHECK(weibull.ok());
    QMOMKIT_VERIFY_CHECK(gamma.ok());
    QMOMKIT_VERIFY_CHECK(weibull.moments.size() == gamma.moments.size());
    for (std::size_t order = 0; order < weibull.moments.size(); ++order)
    {
        assert_near_relative(weibull.moments[order], gamma.moments[order], 1.0e-13);
    }
}

void verify_chi()
{
    constexpr double total_weight = 2.25;
    const double expected_mean = total_weight * 2.0 * std::sqrt(2.0 / std::acos(-1.0));
    const std::vector<double> expected{total_weight, expected_mean, total_weight * 3.0, total_weight * 8.0 * std::sqrt(2.0 / std::acos(-1.0)),
                                       total_weight * 15.0};

    assert_moments(qmomkit::compute_raw_moments(qmomkit::ChiDistribution{3.0}, 4, total_weight), expected, 1.0e-13);
}

void verify_chi_squared_closed_form()
{
    for (const double degrees_of_freedom : {1.0, 4.0, 5.0, 12.0})
    {
        assert_moments_relative(
            qmomkit::ChiSquaredDistribution{degrees_of_freedom}, 1.7, 10,
            [degrees_of_freedom](const std::size_t order) { return chi_squared_raw_moment(degrees_of_freedom, order); }, 4.0e-15);
    }
}

// Even degrees of freedom give integer gamma arguments at even orders and
// half-integer ones at odd orders; odd degrees of freedom the reverse.
void verify_chi_closed_form()
{
    for (const std::size_t degrees_of_freedom : {1, 2, 4, 5, 9})
    {
        assert_moments_relative(
            qmomkit::ChiDistribution{static_cast<double>(degrees_of_freedom)}, 2.25, 10,
            [degrees_of_freedom](const std::size_t order) { return chi_raw_moment(degrees_of_freedom, order); }, 4.0e-15);
    }
}

void verify_gamma_large_shape()
{
    for (const double shape : {10.0, 25.5, 50.0})
    {
        assert_moments_relative(
            qmomkit::GammaDistribution{shape, 0.5}, 3.0, 10, [shape](const std::size_t order) { return gamma_raw_moment(shape, 0.5, order); },
            4.0e-15);
    }
}

void verify_beta_large_parameters()
{
    for (const std::pair<double, double>& parameters : {std::pair<double, double>{20.0, 30.0}, {12.5, 37.5}, {2.0, 48.0}})
    {
        assert_moments_relative(
            qmomkit::BetaDistribution{parameters.first, parameters.second}, 0.8, 10,
            [parameters](const std::size_t order) { return beta_raw_moment(parameters.first, parameters.second, order); }, 4.0e-15);
    }
}

void verify_selected_orders_closed_form()
{
    assert_selected_moments_relative(
        qmomkit::ChiSquaredDistribution{5.0}, 1.7, [](const std::size_t order) { return chi_squared_raw_moment(5.0, order); }, 4.0e-15);
    assert_selected_moments_relative(qmomkit::ChiDistribution{4.0}, 2.25, [](const std::size_t order) { return chi_raw_moment(4, order); }, 4.0e-15);
    assert_selected_moments_relative(qmomkit::ChiDistribution{5.0}, 2.25, [](const std::size_t order) { return chi_raw_moment(5, order); }, 4.0e-15);
    assert_selected_moments_relative(
        qmomkit::GammaDistribution{25.5, 0.5}, 3.0, [](const std::size_t order) { return gamma_raw_moment(25.5, 0.5, order); }, 4.0e-15);
    assert_selected_moments_relative(
        qmomkit::BetaDistribution{20.0, 30.0}, 0.8, [](const std::size_t order) { return beta_raw_moment(20.0, 30.0, order); }, 4.0e-15);
}

// Large parameters, where differences of log-gamma values would lose most of
// their digits. The references are the moments of exactly the double
// parameters the library receives, multiplied out in error-free arithmetic.
void verify_large_parameter_moments_against_exact_products()
{
    constexpr double total_weight = 2.5;
    for (const double shape : {1.0e3, 1.0e5, 1.0e8})
    {
        const double scale = 1.0 / shape;
        const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{shape, scale}, 10, total_weight);
        QMOMKIT_VERIFY_CHECK(result.ok());
        ExactProduct reference;
        reference.multiply(total_weight);
        for (std::size_t order = 1; order <= 10; ++order)
        {
            reference.multiply(shape + static_cast<double>(order - 1));
            reference.multiply(scale);
            assert_relative(result.moments[order], reference.value(), product_tolerance(order));
        }
    }

    for (const double degrees_of_freedom : {1.0e4, 1.0e6})
    {
        const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(qmomkit::ChiSquaredDistribution{degrees_of_freedom}, 10, total_weight);
        QMOMKIT_VERIFY_CHECK(result.ok());
        ExactProduct reference;
        reference.multiply(total_weight);
        for (std::size_t order = 1; order <= 10; ++order)
        {
            reference.multiply(degrees_of_freedom + 2.0 * static_cast<double>(order - 1));
            assert_relative(result.moments[order], reference.value(), product_tolerance(order));
        }
    }

    for (const std::pair<double, double>& parameters : {std::pair<double, double>{3.0e3, 7.0e3}, {3.0e5, 7.0e5}})
    {
        const qmomkit::RawMomentResult result =
            qmomkit::compute_raw_moments(qmomkit::BetaDistribution{parameters.first, parameters.second}, 10, total_weight);
        QMOMKIT_VERIFY_CHECK(result.ok());
        ExactProduct numerator;
        numerator.multiply(total_weight);
        ExactProduct denominator;
        for (std::size_t order = 1; order <= 10; ++order)
        {
            numerator.multiply(parameters.first + static_cast<double>(order - 1));
            denominator.multiply(parameters.first + parameters.second + static_cast<double>(order - 1));
            assert_relative(result.moments[order], exact_quotient(numerator, denominator), product_tolerance(order));
        }
    }

    // Even chi moments: E[X^(2j)] = prod over i < j of (nu + 2i).
    for (const double degrees_of_freedom : {1.0e4, 1000001.0})
    {
        const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(qmomkit::ChiDistribution{degrees_of_freedom}, 10, total_weight);
        QMOMKIT_VERIFY_CHECK(result.ok());
        ExactProduct reference;
        reference.multiply(total_weight);
        for (std::size_t order = 2; order <= 10; order += 2)
        {
            reference.multiply(degrees_of_freedom + static_cast<double>(order - 2));
            assert_relative(result.moments[order], reference.value(), product_tolerance(order));
        }
    }
}

// Odd chi moments need E[X] = sqrt(2) Gamma((nu + 1)/2) / Gamma(nu/2), which
// has no finite product form for general nu. The reference values are the
// doubles nearest to it, computed with mpmath at 40 digits:
//   python3 -c "from mpmath import mp, mpf, gamma, sqrt; mp.dps = 40;
//               nu = mpf('3.3'); print(repr(float(sqrt(2) * gamma((nu + 1) / 2) / gamma(nu / 2))))"
// The higher odd orders follow exactly from E[X^k] = E[X^(k-2)] (nu + k - 2).
void verify_chi_odd_orders_against_high_precision_values()
{
    const std::vector<std::pair<double, double>> first_moments{{0.7, 0.6179099646614185},
                                                               {3.3, 1.6858334189146345},
                                                               {19.5, 4.359652366400788},
                                                               {1234.5, 35.12833758230752},
                                                               {1000001.0, 1000.0002500000312}};
    const std::vector<std::size_t> selected_orders{9, 1, 5};
    for (const auto& [degrees_of_freedom, first_moment] : first_moments)
    {
        const qmomkit::ChiDistribution distribution{degrees_of_freedom};
        const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(distribution, 9);
        QMOMKIT_VERIFY_CHECK(result.ok());
        std::vector<double> references(10, 0.0);
        ExactProduct reference;
        reference.multiply(first_moment);
        for (std::size_t order = 1; order <= 9; order += 2)
        {
            if (order >= 3)
            {
                reference.multiply(degrees_of_freedom + static_cast<double>(order - 2));
            }
            references[order] = reference.value();
            assert_relative(result.moments[order], references[order], product_tolerance(order));
        }

        qmomkit::RawMomentRequest request;
        for (const std::size_t order : selected_orders)
        {
            request.indices.push_back(qmomkit::MultiIndex{order});
        }
        const qmomkit::SelectedRawMomentResult selected = qmomkit::compute_raw_moments(distribution, request);
        QMOMKIT_VERIFY_CHECK(selected.ok());
        for (std::size_t index = 0; index < selected_orders.size(); ++index)
        {
            assert_relative(selected.values[index], references[selected_orders[index]], product_tolerance(selected_orders[index]));
        }
    }
}

// The narrower a distribution, the more digits of its raw moments the
// inversion needs. With moments accurate to a few units in the last place, a
// gamma distribution of shape up to 1e7 inverts at two nodes, up to 1e6 at
// three, and 1e4 at four. At shape 1e7 with three nodes and 1e5 with four, the
// last Wheeler sigma diagonal is zero within its roundoff: the moments no
// longer determine the rule, which then differs from the exact Gauss rule by
// tens of standard deviations, and the fixed-count inversion refuses it.
void verify_narrow_gamma_moments_invert()
{
    const std::vector<std::pair<double, std::size_t>> cases{{1.0e5, 2}, {1.0e5, 3}, {1.0e6, 2}, {1.0e6, 3}, {1.0e7, 2}, {1.0e4, 4}};
    qmomkit::WheelerOptions fixed_options;
    fixed_options.adaptive = false;
    for (const auto& [shape, node_count] : {std::pair{1.0e7, std::size_t{3}}, std::pair{1.0e5, std::size_t{4}}})
    {
        const qmomkit::RawMomentResult moments = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{shape, 1.0 / shape}, 2 * node_count - 1);
        QMOMKIT_VERIFY_CHECK(moments.ok());
        const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments.moments, node_count, fixed_options);
        QMOMKIT_VERIFY_CHECK(quadrature.status == qmomkit::QuadratureStatus::not_realizable);
    }
    for (const auto& [shape, node_count] : cases)
    {
        const qmomkit::RawMomentResult moments = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{shape, 1.0 / shape}, 2 * node_count - 1);
        QMOMKIT_VERIFY_CHECK(moments.ok());
        const qmomkit::Quadrature1D quadrature = qmomkit::compute_quadrature(moments.moments, node_count, fixed_options);
        QMOMKIT_VERIFY_CHECK(quadrature.ok());
        QMOMKIT_VERIFY_CHECK(quadrature.weights.size() == node_count);
        for (std::size_t order = 0; order < moments.moments.size(); ++order)
        {
            double reconstructed = 0.0;
            for (std::size_t node = 0; node < node_count; ++node)
            {
                reconstructed += quadrature.weights[node] * std::pow(quadrature.abscissas[node], static_cast<double>(order));
            }
            assert_relative(reconstructed, moments.moments[order], 1.0e-13);
        }
    }
}

void verify_weibull()
{
    constexpr double total_weight = 5.0;
    std::vector<double> exponential(6, 0.0);
    for (std::size_t order = 0; order < exponential.size(); ++order)
    {
        exponential[order] = total_weight * std::pow(2.0, static_cast<double>(order)) * factorial(order);
    }
    assert_moments(qmomkit::compute_raw_moments(qmomkit::WeibullDistribution{1.0, 2.0}, 5, total_weight), exponential, 1.0e-13);

    const std::vector<double> shape_two{1.0, 3.0 * 0.5 * std::sqrt(std::acos(-1.0)), 9.0, 27.0 * 0.75 * std::sqrt(std::acos(-1.0)), 162.0};
    assert_moments(qmomkit::compute_raw_moments(qmomkit::WeibullDistribution{2.0, 3.0}, 4), shape_two, 1.0e-13);
}

void verify_maximum_order_zero_and_invalid_parameters()
{
    assert_moments(qmomkit::compute_raw_moments(qmomkit::ChiDistribution{2.0}, 0, 6.0), {6.0}, 0.0);
    QMOMKIT_VERIFY_CHECK(qmomkit::compute_raw_moments(qmomkit::NormalDistribution{0.0, -1.0}, 3).status == qmomkit::RawMomentStatus::invalid_input);
    QMOMKIT_VERIFY_CHECK(qmomkit::compute_raw_moments(qmomkit::LogNormalDistribution{0.0, 0.0}, 3).status == qmomkit::RawMomentStatus::invalid_input);
    QMOMKIT_VERIFY_CHECK(qmomkit::compute_raw_moments(qmomkit::BetaDistribution{1.0, 0.0}, 3).status == qmomkit::RawMomentStatus::invalid_input);
    QMOMKIT_VERIFY_CHECK(qmomkit::compute_raw_moments(qmomkit::GammaDistribution{0.0, 1.0}, 3).status == qmomkit::RawMomentStatus::invalid_input);
    QMOMKIT_VERIFY_CHECK(qmomkit::compute_raw_moments(qmomkit::ChiDistribution{0.0}, 3).status == qmomkit::RawMomentStatus::invalid_input);
    QMOMKIT_VERIFY_CHECK(qmomkit::compute_raw_moments(qmomkit::ChiSquaredDistribution{0.0}, 3).status == qmomkit::RawMomentStatus::invalid_input);
    QMOMKIT_VERIFY_CHECK(qmomkit::compute_raw_moments(qmomkit::WeibullDistribution{0.0, 1.0}, 3).status == qmomkit::RawMomentStatus::invalid_input);
}

void verify_configuration_distribution_initialization()
{
    struct Case
    {
        std::string configuration_fields;
        qmomkit::DistributionSpecification distribution;
        std::size_t maximum_order = 0;
        double number_density = 1.0;
    };

    const std::vector<Case> cases{
        {"family = normal\nnumber_density = 2.5\nmaximum_order = 5\nmean = 1.25\nstandard_deviation = 0.75\n",
         qmomkit::NormalDistribution{1.25, 0.75}, 5, 2.5},
        {"family = lognormal\nnumber_density = 3.5\nmaximum_order = 6\nmu_log = -0.2\nsigma_log = 0.4\n", qmomkit::LogNormalDistribution{-0.2, 0.4},
         6, 3.5},
        {"family = gamma\nnumber_density = 4.25\nmaximum_order = 5\nshape = 2.0\nscale = 1.5\n", qmomkit::GammaDistribution{2.0, 1.5}, 5, 4.25},
        {"family = beta\nnumber_density = 1.75\nmaximum_order = 6\nalpha = 2.0\nbeta = 5.0\n", qmomkit::BetaDistribution{2.0, 5.0}, 6, 1.75}};

    for (const Case& test_case : cases)
    {
        const std::string text = "[phase_space]\nnode_counts = 3\n[initial]\ntype = distribution\n" + test_case.configuration_fields;
        const auto configured = qmomkit::config::parse_qmom(text);
        const qmomkit::RawMomentResult direct =
            qmomkit::compute_raw_moments(test_case.distribution, test_case.maximum_order, test_case.number_density);

        QMOMKIT_VERIFY_CHECK(configured.ok());
        QMOMKIT_VERIFY_CHECK(direct.ok());
        QMOMKIT_VERIFY_CHECK(configured.value().moments.size() == test_case.maximum_order + 1);
        assert_near_relative(configured.value().moments[0], test_case.number_density, 0.0);
        assert_moment_sets_match(configured.value().moments, direct.moments, 1.0e-13);
    }
}

void verify_configuration_mixture_initialization()
{
    struct Case
    {
        std::string configuration_components;
        qmomkit::DistributionMixture mixture;
        std::size_t maximum_order = 0;
        double expected_density = 0.0;
    };

    const std::vector<Case> cases{{R"(
[initial.mixture.component.0]
family = normal
number_density = 1.2
mean = -1.0
standard_deviation = 0.5

[initial.mixture.component.1]
family = normal
number_density = 3.4
mean = 2.0
standard_deviation = 1.5
)",
                                   {{{qmomkit::NormalDistribution{-1.0, 0.5}, 1.2}, {qmomkit::NormalDistribution{2.0, 1.5}, 3.4}}},
                                   6,
                                   4.6},
                                  {R"(
[initial.mixture.component.0]
family = gamma
number_density = 1.5
shape = 2.0
scale = 3.0

[initial.mixture.component.1]
family = lognormal
number_density = 0.75
mu_log = 0.25
sigma_log = 0.4
)",
                                   {{{qmomkit::GammaDistribution{2.0, 3.0}, 1.5}, {qmomkit::LogNormalDistribution{0.25, 0.4}, 0.75}}},
                                   5,
                                   2.25},
                                  {R"(
[initial.mixture.component.0]
family = normal
number_density = 2.0
mean = -0.5
standard_deviation = 0.75

[initial.mixture.component.1]
family = beta
number_density = 1.25
alpha = 2.0
beta = 5.0

[initial.mixture.component.2]
family = weibull
number_density = 3.5
shape = 2.0
scale = 1.5
)",
                                   {{{qmomkit::NormalDistribution{-0.5, 0.75}, 2.0},
                                     {qmomkit::BetaDistribution{2.0, 5.0}, 1.25},
                                     {qmomkit::WeibullDistribution{2.0, 1.5}, 3.5}}},
                                   4,
                                   6.75}};

    for (const Case& test_case : cases)
    {
        const std::string text =
            "[phase_space]\nnode_counts = 2\n[initial]\ntype = mixture\nmaximum_order = " + std::to_string(test_case.maximum_order) + "\n" +
            test_case.configuration_components;
        const auto configured = qmomkit::config::parse_qmom(text);
        const qmomkit::RawMomentResult direct = qmomkit::compute_raw_moments(test_case.mixture, test_case.maximum_order);

        QMOMKIT_VERIFY_CHECK(configured.ok());
        QMOMKIT_VERIFY_CHECK(direct.ok());
        QMOMKIT_VERIFY_CHECK(configured.value().moments.size() == test_case.maximum_order + 1);
        assert_near_relative(configured.value().moments[0], test_case.expected_density, 0.0);
        assert_moment_sets_match(configured.value().moments, direct.moments, 1.0e-12);
    }
}

// The subject of these checks is the moment tensor a configuration produces, so
// they parse the initial condition and initialize the tensor directly rather
// than building a CQMOM case, which holds only the moments it reconstructs from.
qmomkit::config::ConfigurationResult<qmomkit::RawMomentTensor> parse_configured_moment_tensor(const std::string_view text)
{
    const auto initial = qmomkit::config::parse_cqmom_initial_condition(text);
    if (!initial.ok())
    {
        return qmomkit::config::ConfigurationResult<qmomkit::RawMomentTensor>::failure(initial.error());
    }
    return qmomkit::config::initialize_cqmom_moment_tensor(initial.value());
}

void verify_configuration_independent_product_two_dimensional()
{
    constexpr double number_density = 2.75;
    const std::vector<std::size_t> extents{7, 3};
    std::vector<double> expected(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_index = 0; flat_index < expected.size(); ++flat_index)
    {
        const qmomkit::MultiIndex order = qmomkit::unflatten_multi_index(flat_index, extents);
        expected[flat_index] = number_density * gamma_raw_moment(1.5, 2.0, order[0]) * gamma_raw_moment(2.5, 0.5, order[1]);
    }

    const auto configured = parse_configured_moment_tensor(R"(
[phase_space]
node_counts = 3, 1

[initial]
type = independent_product
number_density = 2.75
maximum_orders = 6, 2

[initial.dimension.0]
family = gamma
shape = 1.5
scale = 2.0

[initial.dimension.1]
family = gamma
shape = 2.5
scale = 0.5
)");

    QMOMKIT_VERIFY_CHECK(configured.ok());
    assert_tensor_values(configured.value(), extents, expected, 1.0e-12, 1.0e-13);
    assert_near_absolute_relative(configured.value().values[0], number_density, 0.0, 0.0);
}

void verify_configuration_independent_product_normal_gamma()
{
    constexpr double number_density = 1.75;
    constexpr double normal_mean = -0.5;
    constexpr double normal_standard_deviation = 1.25;
    constexpr double gamma_shape = 2.5;
    constexpr double gamma_scale = 0.75;
    const std::vector<std::size_t> extents{4, 4};
    std::vector<double> expected(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_index = 0; flat_index < expected.size(); ++flat_index)
    {
        const qmomkit::MultiIndex order = qmomkit::unflatten_multi_index(flat_index, extents);
        expected[flat_index] = number_density * normal_raw_moment_sum(normal_mean, normal_standard_deviation, order[0]) *
                               gamma_raw_moment(gamma_shape, gamma_scale, order[1]);
    }

    const std::string configuration_text = R"(
[phase_space]
names = x, y
node_counts = 2, 2
conditioning_order = 1, 0

[initial]
type = independent_product
number_density = 1.75
maximum_orders = 3, 3

[initial.dimension.0]
family = normal
mean = -0.5
standard_deviation = 1.25

[initial.dimension.1]
family = gamma
shape = 2.5
scale = 0.75
)";

    const auto configured = parse_configured_moment_tensor(configuration_text);
    QMOMKIT_VERIFY_CHECK(configured.ok());

    // The conditioning order is phase-space configuration rather than moment
    // data, so it is checked on the parsed case rather than on the tensor.
    const auto configured_cqmom_case_result = qmomkit::config::parse_cqmom(configuration_text);
    QMOMKIT_VERIFY_CHECK(configured_cqmom_case_result.ok());
    QMOMKIT_VERIFY_CHECK(configured_cqmom_case_result.value().moments.layout.hierarchy.conditioning_order() == std::vector<std::size_t>({1, 0}));

    assert_tensor_values(configured.value(), extents, expected, 1.0e-12, 1.0e-12);
    assert_near_absolute_relative(configured.value().values[qmomkit::flatten_multi_index({0, 0}, extents)], number_density, 0.0, 0.0);
    assert_near_absolute_relative(configured.value().values[qmomkit::flatten_multi_index({1, 2}, extents)],
                                  number_density * normal_raw_moment_sum(normal_mean, normal_standard_deviation, 1) *
                                      gamma_raw_moment(gamma_shape, gamma_scale, 2),
                                  1.0e-12, 1.0e-12);
}

void verify_configuration_independent_product_beta_log_normal_weibull()
{
    constexpr double number_density = 3.25;
    const std::vector<std::size_t> extents{3, 3, 3};
    std::vector<double> expected(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_index = 0; flat_index < expected.size(); ++flat_index)
    {
        const qmomkit::MultiIndex order = qmomkit::unflatten_multi_index(flat_index, extents);
        expected[flat_index] = number_density * beta_raw_moment(2.0, 5.0, order[0]) * log_normal_raw_moment(0.2, 0.35, order[1]) *
                               weibull_raw_moment(1.5, 2.0, order[2]);
    }

    const auto configured = parse_configured_moment_tensor(R"(
[phase_space]
node_counts = 1, 1, 1

[initial]
type = independent_product
number_density = 3.25
maximum_orders = 2, 2, 2

[initial.dimension.0]
family = beta
alpha = 2.0
beta = 5.0

[initial.dimension.1]
family = lognormal
mu_log = 0.2
sigma_log = 0.35

[initial.dimension.2]
family = weibull
shape = 1.5
scale = 2.0
)");

    QMOMKIT_VERIFY_CHECK(configured.ok());
    assert_tensor_values(configured.value(), extents, expected, 1.0e-11, 1.0e-12);
    assert_near_absolute_relative(configured.value().values[0], number_density, 0.0, 0.0);
    QMOMKIT_VERIFY_CHECK(qmomkit::flatten_multi_index({0, 0, 1}, extents) == 1);
    QMOMKIT_VERIFY_CHECK(qmomkit::flatten_multi_index({0, 1, 0}, extents) == 3);
    QMOMKIT_VERIFY_CHECK(qmomkit::flatten_multi_index({1, 0, 0}, extents) == 9);
}

void verify_configuration_independent_product_dimension_sweep()
{
    for (std::size_t dimensions = 4; dimensions <= 6; ++dimensions)
    {
        std::string node_counts;
        std::string maximum_orders;
        std::string dimension_sections;
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
        {
            if (dimension != 0)
            {
                node_counts += ", ";
                maximum_orders += ", ";
            }
            node_counts += "1";
            maximum_orders += "1";
            dimension_sections += "\n[initial.dimension." + std::to_string(dimension) +
                                  "]\nfamily = gamma\nshape = " + std::to_string(static_cast<double>(dimension + 1)) + "\nscale = 0.5\n";
        }

        const double number_density = 1.0 + static_cast<double>(dimensions);
        const std::string text = "[phase_space]\nnode_counts = " + node_counts +
                                 "\n[initial]\ntype = independent_product\nnumber_density = " + std::to_string(number_density) +
                                 "\nmaximum_orders = " + maximum_orders + "\n" + dimension_sections;
        const auto configured = parse_configured_moment_tensor(text);
        QMOMKIT_VERIFY_CHECK(configured.ok());

        const std::vector<std::size_t> extents(dimensions, 2);
        std::vector<double> expected(qmomkit::product_extents(extents), 0.0);
        for (std::size_t flat_index = 0; flat_index < expected.size(); ++flat_index)
        {
            const qmomkit::MultiIndex order = qmomkit::unflatten_multi_index(flat_index, extents);
            double value = number_density;
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
            {
                value *= gamma_raw_moment(static_cast<double>(dimension + 1), 0.5, order[dimension]);
            }
            expected[flat_index] = value;
        }

        assert_tensor_values(configured.value(), extents, expected, 1.0e-12, 1.0e-12);
        QMOMKIT_VERIFY_CHECK(configured.value().extents.size() == dimensions);
        assert_near_absolute_relative(configured.value().values[0], number_density, 0.0, 0.0);
        qmomkit::MultiIndex last_order(dimensions, 1);
        assert_near_absolute_relative(configured.value().values[qmomkit::flatten_multi_index(last_order, extents)],
                                      expected[qmomkit::flatten_multi_index(last_order, extents)], 1.0e-12, 1.0e-12);
    }
}

void verify_configuration_independent_product_mixture_marginal()
{
    constexpr double number_density = 2.75;
    constexpr double mixture_density = 1.5 + 3.25;
    constexpr double normal_mean = -0.4;
    constexpr double normal_standard_deviation = 0.7;
    const std::string configuration = R"(
[phase_space]
node_counts = 2, 2

[initial]
type = independent_product
number_density = 2.75
maximum_orders = 3, 3

[initial.dimension.0]
family = normal
mean = -0.4
standard_deviation = 0.7

[initial.dimension.1]
type = mixture

[initial.dimension.1.mixture.component.0]
family = lognormal
number_density = 1.5
mu_log = -0.3
sigma_log = 0.25

[initial.dimension.1.mixture.component.1]
family = lognormal
number_density = 3.25
mu_log = 0.4
sigma_log = 0.45
)";

    const auto initial = qmomkit::config::parse_cqmom_initial_condition(configuration);
    QMOMKIT_VERIFY_CHECK(initial.ok());
    const auto* const independent_product = std::get_if<qmomkit::config::IndependentProductInitialization>(&initial.value());
    QMOMKIT_VERIFY_CHECK(independent_product != nullptr);
    QMOMKIT_VERIFY_CHECK(independent_product->distribution.marginals.size() == 2);

    const auto mixture_raw_moment = [](const std::size_t order)
    { return 1.5 * log_normal_raw_moment(-0.3, 0.25, order) + 3.25 * log_normal_raw_moment(0.4, 0.45, order); };
    const auto analytical_moment = [&](const qmomkit::MultiIndex& order)
    {
        return number_density * normal_raw_moment_sum(normal_mean, normal_standard_deviation, order[0]) * mixture_raw_moment(order[1]) /
               mixture_density;
    };

    const std::vector<std::size_t> extents{4, 4};
    std::vector<double> expected(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_index = 0; flat_index < expected.size(); ++flat_index)
    {
        expected[flat_index] = analytical_moment(qmomkit::unflatten_multi_index(flat_index, extents));
    }

    const qmomkit::TensorMomentResult tensor_result = qmomkit::compute_raw_moments(independent_product->distribution, {3, 3});
    assert_tensor(tensor_result, extents, expected, 1.0e-12, 1.0e-12);
    assert_near_absolute_relative(tensor_result.moments.values[qmomkit::flatten_multi_index({0, 0}, extents)], number_density, 0.0, 0.0);

    const qmomkit::RawMomentRequest request{{{0, 0}, {1, 0}, {0, 1}, {1, 1}, {0, 2}, {2, 3}}};
    const qmomkit::SelectedRawMomentResult selected = qmomkit::compute_raw_moments(independent_product->distribution, request);
    QMOMKIT_VERIFY_CHECK(selected.ok());
    QMOMKIT_VERIFY_CHECK(selected.indices == request.indices);
    for (std::size_t index = 0; index < request.indices.size(); ++index)
    {
        assert_near_absolute_relative(selected.values[index], analytical_moment(request.indices[index]), 1.0e-12, 1.0e-12);
    }

    const auto configured = qmomkit::config::parse_cqmom_prepared(configuration);
    QMOMKIT_VERIFY_CHECK(configured.ok());
    const qmomkit::CqmomState initial_state = configured.value().initial_state();
    QMOMKIT_VERIFY_CHECK(initial_state.moment_indices().size() == initial_state.moment_values().size());
    for (std::size_t index = 0; index < initial_state.moment_values().size(); ++index)
    {
        assert_near_absolute_relative(initial_state.moment_values()[index], analytical_moment(initial_state.moment_indices()[index]), 1.0e-12,
                                      1.0e-12);
    }

    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(configured.value());
    QMOMKIT_VERIFY_CHECK(cqmom.ok());
    const qmomkit::SelectedRawMomentResult reconstructed =
        qmomkit::reconstruct_raw_moments(cqmom.quadrature, qmomkit::RawMomentRequest{initial_state.moment_indices()});
    QMOMKIT_VERIFY_CHECK(reconstructed.ok());
    QMOMKIT_VERIFY_CHECK(reconstructed.indices == initial_state.moment_indices());
    for (std::size_t index = 0; index < reconstructed.values.size(); ++index)
    {
        assert_near_absolute_relative(reconstructed.values[index], initial_state.moment_values()[index], 5.0e-10, 5.0e-10);
    }
}

void verify_one_component_mixture_equivalence()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::GammaDistribution{2.0, 3.0}, 2.5}}};
    const qmomkit::RawMomentResult mixture_result = qmomkit::compute_raw_moments(mixture, 6);
    const qmomkit::RawMomentResult distribution_result = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{2.0, 3.0}, 6, 2.5);

    QMOMKIT_VERIFY_CHECK(mixture_result.ok());
    QMOMKIT_VERIFY_CHECK(distribution_result.ok());
    QMOMKIT_VERIFY_CHECK(mixture_result.moments.size() == distribution_result.moments.size());
    for (std::size_t order = 0; order < mixture_result.moments.size(); ++order)
    {
        assert_near_relative(mixture_result.moments[order], distribution_result.moments[order], 1.0e-13);
    }
}

void verify_bimodal_normal_mixture()
{
    constexpr double density_a = 1.2;
    constexpr double density_b = 3.4;
    std::vector<double> expected(7, 0.0);
    for (std::size_t order = 0; order < expected.size(); ++order)
    {
        expected[order] = density_a * normal_raw_moment_sum(-1.0, 0.5, order) + density_b * normal_raw_moment_sum(2.0, 1.5, order);
    }

    const qmomkit::DistributionMixture mixture{
        {{qmomkit::NormalDistribution{-1.0, 0.5}, density_a}, {qmomkit::NormalDistribution{2.0, 1.5}, density_b}}};
    assert_moments(qmomkit::compute_raw_moments(mixture, 6), expected, 1.0e-12);
}

void verify_heterogeneous_gamma_log_normal_mixture()
{
    constexpr double gamma_density = 1.5;
    constexpr double log_normal_density = 0.75;
    std::vector<double> expected(6, 0.0);
    for (std::size_t order = 0; order < expected.size(); ++order)
    {
        expected[order] = gamma_density * gamma_raw_moment(2.0, 3.0, order) + log_normal_density * log_normal_raw_moment(0.25, 0.4, order);
    }

    const qmomkit::DistributionMixture mixture{
        {{qmomkit::GammaDistribution{2.0, 3.0}, gamma_density}, {qmomkit::LogNormalDistribution{0.25, 0.4}, log_normal_density}}};
    assert_moments(qmomkit::compute_raw_moments(mixture, 5), expected, 1.0e-12);
}

void verify_three_component_multimodal_mixture()
{
    std::vector<double> expected(5, 0.0);
    for (std::size_t order = 0; order < expected.size(); ++order)
    {
        expected[order] =
            2.0 * normal_raw_moment_sum(-0.5, 0.75, order) + 1.25 * beta_raw_moment(2.0, 5.0, order) + 3.5 * weibull_raw_moment(2.0, 1.5, order);
    }

    const qmomkit::DistributionMixture mixture{
        {{qmomkit::NormalDistribution{-0.5, 0.75}, 2.0}, {qmomkit::BetaDistribution{2.0, 5.0}, 1.25}, {qmomkit::WeibullDistribution{2.0, 1.5}, 3.5}}};
    assert_moments(qmomkit::compute_raw_moments(mixture, 4), expected, 1.0e-12);
}

void verify_zero_density_mixture_component_policy()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::NormalDistribution{0.0, 0.0}, 0.0}, {qmomkit::GammaDistribution{1.0, 2.0}, 3.0}}};
    const std::vector<double> expected{3.0, 6.0, 24.0, 144.0};

    assert_moments(qmomkit::compute_raw_moments(mixture, 3), expected, 1.0e-13);
}

void verify_strongly_unequal_component_densities()
{
    constexpr double large_density = 1.0e12;
    constexpr double small_density = 1.0;
    std::vector<double> expected(5, 0.0);
    for (std::size_t order = 0; order < expected.size(); ++order)
    {
        expected[order] = large_density * beta_raw_moment(2.0, 3.0, order) + small_density * gamma_raw_moment(3.0, 0.5, order);
    }

    const qmomkit::DistributionMixture mixture{
        {{qmomkit::BetaDistribution{2.0, 3.0}, large_density}, {qmomkit::GammaDistribution{3.0, 0.5}, small_density}}};
    assert_moments(qmomkit::compute_raw_moments(mixture, 4), expected, 1.0e-13);
}

void verify_normal_odd_moment_cancellation()
{
    constexpr double density = 2.0;
    std::vector<double> expected(6, 0.0);
    for (std::size_t order = 0; order < expected.size(); ++order)
    {
        expected[order] = density * normal_raw_moment_sum(-10.0, 0.5, order) + density * normal_raw_moment_sum(10.0, 0.5, order);
    }

    const qmomkit::DistributionMixture mixture{
        {{qmomkit::NormalDistribution{-10.0, 0.5}, density}, {qmomkit::NormalDistribution{10.0, 0.5}, density}}};
    const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(mixture, 5);

    assert_moments(result, expected, 1.0e-12);
    assert_near_relative(result.moments[1], 0.0, 1.0e-14);
    assert_near_relative(result.moments[3], 0.0, 1.0e-12);
    assert_near_relative(result.moments[5], 0.0, 1.0e-10);
}

void verify_component_order_permutation_invariance()
{
    const qmomkit::DistributionMixture original{{{qmomkit::GammaDistribution{2.0, 1.5}, 2.0},
                                                 {qmomkit::NormalDistribution{-0.25, 0.75}, 3.0},
                                                 {qmomkit::LogNormalDistribution{0.1, 0.2}, 4.0}}};
    const qmomkit::DistributionMixture permuted{{{qmomkit::LogNormalDistribution{0.1, 0.2}, 4.0},
                                                 {qmomkit::GammaDistribution{2.0, 1.5}, 2.0},
                                                 {qmomkit::NormalDistribution{-0.25, 0.75}, 3.0}}};

    const qmomkit::RawMomentResult original_result = qmomkit::compute_raw_moments(original, 6);
    const qmomkit::RawMomentResult permuted_result = qmomkit::compute_raw_moments(permuted, 6);

    QMOMKIT_VERIFY_CHECK(original_result.ok());
    QMOMKIT_VERIFY_CHECK(permuted_result.ok());
    QMOMKIT_VERIFY_CHECK(original_result.moments.size() == permuted_result.moments.size());
    for (std::size_t order = 0; order < original_result.moments.size(); ++order)
    {
        assert_near_relative(permuted_result.moments[order], original_result.moments[order], 1.0e-14);
    }
}

void verify_large_total_density_is_not_normalized()
{
    constexpr double beta_density = 250000.0;
    constexpr double gamma_density = 750000.0;
    std::vector<double> expected(4, 0.0);
    for (std::size_t order = 0; order < expected.size(); ++order)
    {
        expected[order] = beta_density * beta_raw_moment(1.0, 1.0, order) + gamma_density * gamma_raw_moment(1.0, 1.0, order);
    }

    const qmomkit::DistributionMixture mixture{
        {{qmomkit::BetaDistribution{1.0, 1.0}, beta_density}, {qmomkit::GammaDistribution{1.0, 1.0}, gamma_density}}};
    assert_moments(qmomkit::compute_raw_moments(mixture, 3), expected, 1.0e-13);
}

void verify_mixture_accumulation_overflow()
{
    const double large_mean = 0.75 * std::numeric_limits<double>::max();
    const qmomkit::DistributionMixture mixture{
        {{qmomkit::NormalDistribution{large_mean, 1.0}, 1.0}, {qmomkit::NormalDistribution{large_mean, 1.0}, 1.0}}};
    const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(mixture, 1);

    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::RawMomentStatus::overflow);
    QMOMKIT_VERIFY_CHECK(result.failed_order == 1);
    QMOMKIT_VERIFY_CHECK(result.failed_component == 1);
}

void verify_invalid_mixture_density()
{
    const qmomkit::DistributionMixture mixture{{{qmomkit::GammaDistribution{1.0, 1.0}, -1.0}}};
    const qmomkit::RawMomentResult result = qmomkit::compute_raw_moments(mixture, 1);

    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::RawMomentStatus::invalid_input);
    QMOMKIT_VERIFY_CHECK(result.failed_component == 0);
}

void verify_independent_product_one_dimensional_equivalence()
{
    constexpr double number_density = 2.5;
    const qmomkit::RawMomentResult one_dimensional = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{1.5, 2.0}, 6, number_density);
    const qmomkit::TensorMomentResult tensor =
        qmomkit::compute_raw_moments(qmomkit::IndependentProductDistribution{{qmomkit::GammaDistribution{1.5, 2.0}}, number_density}, {6});

    QMOMKIT_VERIFY_CHECK(one_dimensional.ok());
    QMOMKIT_VERIFY_CHECK(tensor.ok());
    QMOMKIT_VERIFY_CHECK(tensor.moments.extents == std::vector<std::size_t>{7});
    QMOMKIT_VERIFY_CHECK(tensor.moments.values.size() == one_dimensional.moments.size());
    for (std::size_t order = 0; order < one_dimensional.moments.size(); ++order)
    {
        assert_near_absolute_relative(tensor.moments.values[order], one_dimensional.moments[order], 1.0e-12, 1.0e-13);
    }
}

void verify_independent_product_normal_gamma()
{
    constexpr double number_density = 1.75;
    constexpr double normal_mean = -0.5;
    constexpr double normal_standard_deviation = 1.25;
    constexpr double gamma_shape = 2.5;
    constexpr double gamma_scale = 0.75;
    const std::vector<std::size_t> extents{5, 4};
    std::vector<double> expected(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_index = 0; flat_index < expected.size(); ++flat_index)
    {
        const qmomkit::MultiIndex order = qmomkit::unflatten_multi_index(flat_index, extents);
        expected[flat_index] = number_density * normal_raw_moment_sum(normal_mean, normal_standard_deviation, order[0]) *
                               gamma_raw_moment(gamma_shape, gamma_scale, order[1]);
    }

    const qmomkit::TensorMomentResult result =
        qmomkit::compute_raw_moments(qmomkit::IndependentProductDistribution{{qmomkit::NormalDistribution{normal_mean, normal_standard_deviation},
                                                                              qmomkit::GammaDistribution{gamma_shape, gamma_scale}},
                                                                             number_density},
                                     {4, 3});

    assert_tensor(result, extents, expected, 1.0e-12, 1.0e-12);
    assert_near_absolute_relative(result.moments.values[qmomkit::flatten_multi_index({0, 0}, extents)], number_density, 0.0, 0.0);
    assert_near_absolute_relative(result.moments.values[qmomkit::flatten_multi_index({1, 2}, extents)],
                                  expected[qmomkit::flatten_multi_index({1, 2}, extents)], 1.0e-12, 1.0e-12);
    assert_near_absolute_relative(result.moments.values[qmomkit::flatten_multi_index({4, 3}, extents)],
                                  expected[qmomkit::flatten_multi_index({4, 3}, extents)], 1.0e-12, 1.0e-12);
}

void verify_independent_product_beta_log_normal_weibull()
{
    constexpr double number_density = 3.25;
    const std::vector<std::size_t> extents{4, 3, 5};
    std::vector<double> expected(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_index = 0; flat_index < expected.size(); ++flat_index)
    {
        const qmomkit::MultiIndex order = qmomkit::unflatten_multi_index(flat_index, extents);
        expected[flat_index] = number_density * beta_raw_moment(2.0, 5.0, order[0]) * log_normal_raw_moment(0.2, 0.35, order[1]) *
                               weibull_raw_moment(1.5, 2.0, order[2]);
    }

    const qmomkit::TensorMomentResult result = qmomkit::compute_raw_moments(
        qmomkit::IndependentProductDistribution{
            {qmomkit::BetaDistribution{2.0, 5.0}, qmomkit::LogNormalDistribution{0.2, 0.35}, qmomkit::WeibullDistribution{1.5, 2.0}}, number_density},
        {3, 2, 4});

    assert_tensor(result, extents, expected, 1.0e-11, 1.0e-12);
    assert_near_absolute_relative(result.moments.values[0], number_density, 0.0, 0.0);
    assert_near_absolute_relative(result.moments.values[qmomkit::flatten_multi_index({0, 2, 0}, extents)],
                                  expected[qmomkit::flatten_multi_index({0, 2, 0}, extents)], 1.0e-12, 1.0e-12);
    QMOMKIT_VERIFY_CHECK(qmomkit::flatten_multi_index({0, 0, 1}, extents) == 1);
    QMOMKIT_VERIFY_CHECK(qmomkit::flatten_multi_index({0, 1, 0}, extents) == 5);
    QMOMKIT_VERIFY_CHECK(qmomkit::flatten_multi_index({1, 0, 0}, extents) == 15);
}

void verify_independent_product_signed_odd_moments()
{
    constexpr double number_density = 2.0;
    const std::vector<std::size_t> extents{4, 3};
    std::vector<double> expected(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_index = 0; flat_index < expected.size(); ++flat_index)
    {
        const qmomkit::MultiIndex order = qmomkit::unflatten_multi_index(flat_index, extents);
        expected[flat_index] = number_density * normal_raw_moment_sum(-2.0, 0.5, order[0]) * normal_raw_moment_sum(0.0, 1.0, order[1]);
    }

    const qmomkit::TensorMomentResult result = qmomkit::compute_raw_moments(
        qmomkit::IndependentProductDistribution{{qmomkit::NormalDistribution{-2.0, 0.5}, qmomkit::NormalDistribution{0.0, 1.0}}, number_density},
        {3, 2});

    assert_tensor(result, extents, expected, 1.0e-12, 1.0e-12);
    QMOMKIT_VERIFY_CHECK(result.moments.values[qmomkit::flatten_multi_index({1, 0}, extents)] < 0.0);
    assert_near_absolute_relative(result.moments.values[qmomkit::flatten_multi_index({0, 1}, extents)], 0.0, 0.0, 0.0);
    assert_near_absolute_relative(result.moments.values[qmomkit::flatten_multi_index({1, 1}, extents)], 0.0, 0.0, 0.0);
}

void verify_independent_product_zero_density()
{
    const qmomkit::TensorMomentResult result = qmomkit::compute_raw_moments(
        qmomkit::IndependentProductDistribution{{qmomkit::NormalDistribution{0.0, 0.0}, qmomkit::GammaDistribution{0.0, 1.0}}, 0.0}, {2, 3});

    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.moments.extents == std::vector<std::size_t>({3, 4}));
    QMOMKIT_VERIFY_CHECK(result.moments.values.size() == 12);
    for (double value : result.moments.values)
    {
        assert_near_absolute_relative(value, 0.0, 0.0, 0.0);
    }
}

void verify_independent_product_tensor_overflow()
{
    const double large_mean = 0.75 * std::numeric_limits<double>::max();
    const qmomkit::TensorMomentResult result = qmomkit::compute_raw_moments(
        qmomkit::IndependentProductDistribution{{qmomkit::NormalDistribution{large_mean, 1.0}, qmomkit::NormalDistribution{large_mean, 1.0}}, 1.0},
        {1, 1});

    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::RawMomentStatus::overflow);
    QMOMKIT_VERIFY_CHECK(result.failed_dimension == qmomkit::invalid_dimension_index);
    QMOMKIT_VERIFY_CHECK(result.failed_order == std::vector<std::size_t>({1, 1}));
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"standard_normal", verify_standard_normal},
            {"shifted_scaled_normal", verify_shifted_scaled_normal},
            {"log_normal", verify_log_normal},
            {"beta", verify_beta},
            {"gamma", verify_gamma},
            {"chi_squared", verify_chi_squared},
            {"weibull_gamma_identity", verify_weibull_gamma_identity},
            {"chi", verify_chi},
            {"chi_closed_form", verify_chi_closed_form},
            {"chi_squared_closed_form", verify_chi_squared_closed_form},
            {"gamma_large_shape", verify_gamma_large_shape},
            {"beta_large_parameters", verify_beta_large_parameters},
            {"selected_orders_closed_form", verify_selected_orders_closed_form},
            {"large_parameter_moments_against_exact_products", verify_large_parameter_moments_against_exact_products},
            {"chi_odd_orders_against_high_precision_values", verify_chi_odd_orders_against_high_precision_values},
            {"narrow_gamma_moments_invert", verify_narrow_gamma_moments_invert},
            {"weibull", verify_weibull},
            {"maximum_order_zero_and_invalid_parameters", verify_maximum_order_zero_and_invalid_parameters},
            {"configuration_distribution_initialization", verify_configuration_distribution_initialization},
            {"configuration_mixture_initialization", verify_configuration_mixture_initialization},
            {"configuration_independent_product_two_dimensional", verify_configuration_independent_product_two_dimensional},
            {"configuration_independent_product_normal_gamma", verify_configuration_independent_product_normal_gamma},
            {"configuration_independent_product_beta_log_normal_weibull", verify_configuration_independent_product_beta_log_normal_weibull},
            {"configuration_independent_product_dimension_sweep", verify_configuration_independent_product_dimension_sweep},
            {"configuration_independent_product_mixture_marginal", verify_configuration_independent_product_mixture_marginal},
            {"one_component_mixture_equivalence", verify_one_component_mixture_equivalence},
            {"bimodal_normal_mixture", verify_bimodal_normal_mixture},
            {"heterogeneous_gamma_log_normal_mixture", verify_heterogeneous_gamma_log_normal_mixture},
            {"three_component_multimodal_mixture", verify_three_component_multimodal_mixture},
            {"zero_density_mixture_component_policy", verify_zero_density_mixture_component_policy},
            {"strongly_unequal_component_densities", verify_strongly_unequal_component_densities},
            {"normal_odd_moment_cancellation", verify_normal_odd_moment_cancellation},
            {"component_order_permutation_invariance", verify_component_order_permutation_invariance},
            {"large_total_density_is_not_normalized", verify_large_total_density_is_not_normalized},
            {"mixture_accumulation_overflow", verify_mixture_accumulation_overflow},
            {"invalid_mixture_density", verify_invalid_mixture_density},
            {"independent_product_one_dimensional_equivalence", verify_independent_product_one_dimensional_equivalence},
            {"independent_product_normal_gamma", verify_independent_product_normal_gamma},
            {"independent_product_beta_log_normal_weibull", verify_independent_product_beta_log_normal_weibull},
            {"independent_product_signed_odd_moments", verify_independent_product_signed_odd_moments},
            {"independent_product_zero_density", verify_independent_product_zero_density},
            {"independent_product_tensor_overflow", verify_independent_product_tensor_overflow},
        });
}
