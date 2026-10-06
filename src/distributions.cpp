#include "qmomkit/distributions.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace qmomkit
{
namespace
{

[[nodiscard]] bool is_positive_finite(const double value)
{
    return std::isfinite(value) && value > 0.0;
}

[[nodiscard]] RawMomentResult make_failure(const RawMomentStatus status, std::string message, const std::size_t failed_order,
                                           const std::size_t failed_component = invalid_component_index)
{
    RawMomentResult result;
    result.status = status;
    result.message = std::move(message);
    result.failed_order = failed_order;
    result.failed_component = failed_component;
    return result;
}

[[nodiscard]] RawMomentResult invalid_input(std::string message)
{
    return make_failure(RawMomentStatus::invalid_input, std::move(message), 0);
}

[[nodiscard]] TensorMomentResult make_tensor_failure(const RawMomentStatus status, std::string message,
                                                     const std::size_t failed_dimension = invalid_dimension_index, MultiIndex failed_order = {})
{
    TensorMomentResult result;
    result.status = status;
    result.message = std::move(message);
    result.failed_dimension = failed_dimension;
    result.failed_order = std::move(failed_order);
    return result;
}

[[nodiscard]] SelectedRawMomentResult make_selected_failure(const RawMomentStatus status, std::string message,
                                                            const std::size_t failed_dimension = invalid_dimension_index,
                                                            MultiIndex failed_order = {})
{
    SelectedRawMomentResult result;
    result.status = status;
    result.message = std::move(message);
    result.failed_dimension = failed_dimension;
    result.failed_order = std::move(failed_order);
    return result;
}

[[nodiscard]] SelectedRawMomentResult make_selected_component_failure(const RawMomentStatus status, std::string message,
                                                                      const std::size_t failed_component,
                                                                      const std::size_t failed_dimension = invalid_dimension_index,
                                                                      MultiIndex failed_order = {})
{
    SelectedRawMomentResult result = make_selected_failure(status, std::move(message), failed_dimension, std::move(failed_order));
    result.failed_component = failed_component;
    return result;
}

[[nodiscard]] TensorMomentResult invalid_tensor_input(std::string message)
{
    return make_tensor_failure(RawMomentStatus::invalid_input, std::move(message));
}

[[nodiscard]] RawMomentResult failure_at_order(const RawMomentStatus status, const char* reason, const std::size_t order)
{
    return make_failure(status, std::string(reason), order);
}

[[nodiscard]] SelectedRawMomentResult selected_failure_at_order(const RawMomentStatus status, const char* reason, const std::size_t order)
{
    return make_selected_failure(status, std::string(reason), 0, {order});
}

[[nodiscard]] RawMomentStatus finite_failure_status(const double value)
{
    if (std::isinf(value))
    {
        return RawMomentStatus::overflow;
    }
    return RawMomentStatus::non_finite;
}

[[nodiscard]] RawMomentResult check_common_inputs(const std::size_t maximum_order, const double total_weight)
{
    if (maximum_order == std::numeric_limits<std::size_t>::max())
    {
        return invalid_input("Maximum moment order is too large.");
    }
    if (!std::isfinite(total_weight) || total_weight < 0.0)
    {
        return invalid_input("Total weight must be finite and non-negative.");
    }
    return {};
}

[[nodiscard]] RawMomentResult check_maximum_order(const std::size_t maximum_order)
{
    if (maximum_order == std::numeric_limits<std::size_t>::max())
    {
        return invalid_input("Maximum moment order is too large.");
    }
    return {};
}

[[nodiscard]] RawMomentResult finish_moments(std::vector<double> moments)
{
    return RawMomentResult{MomentSet(std::move(moments)), RawMomentStatus::success, "", 0};
}

struct MomentValueResult
{
    double value = 0.0;
    RawMomentStatus status = RawMomentStatus::success;
    const char* message = "";
    std::size_t failed_order = 0;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == RawMomentStatus::success;
    }
};

[[nodiscard]] MomentValueResult moment_value_failure(const RawMomentStatus status, const char* message, const std::size_t failed_order)
{
    return MomentValueResult{0.0, status, message, failed_order};
}

[[nodiscard]] MomentValueResult normalize_marginal_moment(const double raw_moment, const double zeroth_moment, const std::size_t order)
{
    if (!std::isfinite(zeroth_moment) || zeroth_moment <= 0.0)
    {
        return moment_value_failure(RawMomentStatus::invalid_input, "Independent product marginal zeroth raw moment must be finite and positive.",
                                    order);
    }

    const double normalized_moment = raw_moment / zeroth_moment;
    if (!std::isfinite(normalized_moment))
    {
        return moment_value_failure(finite_failure_status(normalized_moment), "Normalized marginal raw moment is not finite.", order);
    }
    if (raw_moment != 0.0 && normalized_moment == 0.0)
    {
        return moment_value_failure(RawMomentStatus::underflow, "Normalized marginal raw moment underflows double precision.", order);
    }
    return MomentValueResult{normalized_moment, RawMomentStatus::success, "", 0};
}

[[nodiscard]] MomentValueResult evaluate_log_expectation(const double log_expectation, const double total_weight, const std::size_t order)
{
    if (order == 0)
    {
        return MomentValueResult{total_weight, RawMomentStatus::success, "", 0};
    }

    const double log_total_weight = std::log(total_weight);
    const double log_max = std::log(std::numeric_limits<double>::max());
    const double log_min = std::log(std::numeric_limits<double>::denorm_min());
    const double log_moment = log_total_weight + log_expectation;
    if (!std::isfinite(log_moment))
    {
        return moment_value_failure(RawMomentStatus::non_finite, "Raw moment logarithm is not finite.", order);
    }
    if (log_moment > log_max)
    {
        return moment_value_failure(RawMomentStatus::overflow, "Raw moment overflows double precision.", order);
    }
    if (log_moment < log_min)
    {
        return moment_value_failure(RawMomentStatus::underflow, "Positive raw moment underflows double precision.", order);
    }

    const double moment = std::exp(log_moment);
    if (!std::isfinite(moment))
    {
        return moment_value_failure(finite_failure_status(moment), "Raw moment is not finite.", order);
    }
    return MomentValueResult{moment, RawMomentStatus::success, "", 0};
}

[[nodiscard]] MomentValueResult evaluate_log_normal_moment(const LogNormalDistribution& distribution, const double variance,
                                                           const double total_weight, const std::size_t order)
{
    // The variance enters every positive order, so an overflowing variance is
    // reported at the first such order that is actually requested.
    if (order > 0 && !std::isfinite(variance))
    {
        return moment_value_failure(RawMomentStatus::overflow, "Underlying log-normal variance overflows double precision.", order);
    }
    const double k = static_cast<double>(order);
    return evaluate_log_expectation(k * distribution.normal_mean + 0.5 * k * k * variance, total_weight, order);
}

// A positive product held as mantissa * 2^exponent and renormalized after
// every factor, so that a long product of large or small factors cannot
// overflow or underflow before its final value is known. Each factor adds one
// rounding to its own evaluation, so a moment that is a product of k simple
// factors is accurate to a few times k units in the last place, however large
// the distribution parameters are.
class ScaledProduct
{
public:
    explicit ScaledProduct(const double initial)
    {
        int exponent = 0;
        mantissa_ = std::frexp(initial, &exponent);
        exponent_ = exponent;
    }

    void multiply(const double factor)
    {
        int factor_exponent = 0;
        const double factor_mantissa = std::frexp(factor, &factor_exponent);
        int product_exponent = 0;
        mantissa_ = std::frexp(mantissa_ * factor_mantissa, &product_exponent);
        exponent_ += static_cast<long long>(factor_exponent) + product_exponent;
    }

    // Subnormal values are returned as they are; only a value outside the
    // range of a double fails.
    [[nodiscard]] MomentValueResult value(const std::size_t order) const
    {
        // The mantissa is below one, so the value is finite exactly when the
        // exponent does not exceed max_exponent.
        if (mantissa_ != 0.0 && exponent_ > std::numeric_limits<double>::max_exponent)
        {
            return moment_value_failure(RawMomentStatus::overflow, "Raw moment overflows double precision.", order);
        }
        constexpr long long lowest_exponent = std::numeric_limits<double>::min_exponent - std::numeric_limits<double>::digits - 1;
        const double value = std::ldexp(mantissa_, static_cast<int>(std::max(exponent_, lowest_exponent)));
        if (value == 0.0)
        {
            return moment_value_failure(RawMomentStatus::underflow, "Positive raw moment underflows double precision.", order);
        }
        return MomentValueResult{value, RawMomentStatus::success, "", 0};
    }

private:
    double mantissa_ = 0.0;
    long long exponent_ = 0;
};

// The factor that takes each family from E[X^(k-1)] to E[X^k]:
//   beta         E[X^k] = prod over i < k of (alpha + i) / (alpha + beta + i)
//   gamma        E[X^k] = prod over i < k of scale (shape + i)
//   chi-squared  E[X^k] = prod over i < k of (nu + 2i)
void multiply_beta_step(ScaledProduct& product, const BetaDistribution& distribution, const std::size_t order)
{
    const double i = static_cast<double>(order - 1);
    product.multiply((distribution.alpha + i) / (distribution.alpha + distribution.beta + i));
}

void multiply_gamma_step(ScaledProduct& product, const GammaDistribution& distribution, const std::size_t order)
{
    product.multiply(distribution.shape + static_cast<double>(order - 1));
    product.multiply(distribution.scale);
}

void multiply_chi_squared_step(ScaledProduct& product, const ChiSquaredDistribution& distribution, const std::size_t order)
{
    product.multiply(distribution.degrees_of_freedom + 2.0 * static_cast<double>(order - 1));
}

// total_weight * E[X^order] for a family whose moments are products of one
// factor per order.
template <typename Step>
[[nodiscard]] MomentValueResult evaluate_product_moment(const double total_weight, const std::size_t order, const Step& multiply_step)
{
    ScaledProduct product(total_weight);
    for (std::size_t current = 1; current <= order; ++current)
    {
        multiply_step(product, current);
    }
    return product.value(order);
}

// Gamma(x + 1/2) / Gamma(x) for x > 0. Subtracting the expansions of
// ln Gamma(x + h) at h = 1/2 and h = 0 (DLMF 5.11.8), with
// B_m(1/2) = -(1 - 2^(1-m)) B_m (DLMF 24.4.27), gives
//   ln(r(x) / sqrt(x)) ~ sum over odd n of (2^-n - 2) B_(n+1) / (n (n+1) x^n),
// whose coefficients are -1/8, 1/192, -1/640, ... From x = 10 on, the eight
// terms below leave a truncation error under 1e-17. Smaller x are first raised
// to that range by r(x) = r(x + 1) x / (x + 1/2).
[[nodiscard]] double half_step_gamma_ratio(double x)
{
    double factor = 1.0;
    while (x < 10.0)
    {
        factor *= x / (x + 0.5);
        x += 1.0;
    }
    constexpr std::array<double, 8> coefficients{-1.0 / 8.0,      1.0 / 192.0,      -1.0 / 640.0,       17.0 / 14336.0,
                                                 -31.0 / 18432.0, 691.0 / 180224.0, -5461.0 / 425984.0, 929569.0 / 15728640.0};
    const double inverse_square = 1.0 / (x * x);
    double series = 0.0;
    for (auto coefficient = coefficients.rbegin(); coefficient != coefficients.rend(); ++coefficient)
    {
        series = series * inverse_square + *coefficient;
    }
    return factor * std::sqrt(x) * std::exp(series / x);
}

// E[X^k] = 2^(k/2) Gamma((nu + k)/2) / Gamma(nu/2) satisfies
// E[X^k] = E[X^(k-2)] (nu + k - 2), so even orders are products from
// E[X^0] = 1 and odd orders from E[X^1] = sqrt(2) Gamma(nu/2 + 1/2) / Gamma(nu/2).
[[nodiscard]] double chi_first_moment(const ChiDistribution& distribution)
{
    constexpr double square_root_two = 1.41421356237309504880;
    return square_root_two * half_step_gamma_ratio(0.5 * distribution.degrees_of_freedom);
}

void multiply_chi_step(ScaledProduct& product, const ChiDistribution& distribution, const std::size_t order)
{
    product.multiply(distribution.degrees_of_freedom + static_cast<double>(order - 2));
}

[[nodiscard]] MomentValueResult evaluate_chi_moment(const ChiDistribution& distribution, const double first_moment, const double total_weight,
                                                    const std::size_t order)
{
    ScaledProduct product(total_weight);
    if (order % 2 == 1)
    {
        product.multiply(first_moment);
    }
    for (std::size_t current = order % 2 + 2; current <= order; current += 2)
    {
        multiply_chi_step(product, distribution, current);
    }
    return product.value(order);
}

[[nodiscard]] MomentValueResult evaluate_weibull_moment(const WeibullDistribution& distribution, const double log_scale, const double total_weight,
                                                        const std::size_t order)
{
    const double k = static_cast<double>(order);
    return evaluate_log_expectation(k * log_scale + std::lgamma(1.0 + k / distribution.shape), total_weight, order);
}

// Each moment is the previous one times the family's step factor, which is
// the same arithmetic evaluate_product_moment performs for a single order.
template <typename Step>
[[nodiscard]] RawMomentResult compute_from_products(const std::size_t maximum_order, const double total_weight, const Step& multiply_step)
{
    if (total_weight == 0.0)
    {
        return finish_moments(std::vector<double>(maximum_order + 1, 0.0));
    }

    std::vector<double> moments(maximum_order + 1, 0.0);
    moments[0] = total_weight;
    ScaledProduct product(total_weight);
    for (std::size_t order = 1; order <= maximum_order; ++order)
    {
        multiply_step(product, order);
        const MomentValueResult moment = product.value(order);
        if (!moment.ok())
        {
            return failure_at_order(moment.status, moment.message, moment.failed_order);
        }
        moments[order] = moment.value;
    }
    return finish_moments(std::move(moments));
}

template <typename MomentEvaluator>
[[nodiscard]] RawMomentResult compute_from_log_moments(const std::size_t maximum_order, const double total_weight, MomentEvaluator&& evaluate_moment)
{
    if (total_weight == 0.0)
    {
        return finish_moments(std::vector<double>(maximum_order + 1, 0.0));
    }

    std::vector<double> moments(maximum_order + 1, 0.0);
    moments[0] = total_weight;
    for (std::size_t order = 1; order <= maximum_order; ++order)
    {
        const MomentValueResult moment = evaluate_moment(order);
        if (!moment.ok())
        {
            return failure_at_order(moment.status, moment.message, moment.failed_order);
        }
        moments[order] = moment.value;
    }
    return finish_moments(std::move(moments));
}

[[nodiscard]] SelectedRawMomentResult check_selected_common_inputs(const RawMomentRequest& request, const double total_weight,
                                                                   std::vector<std::size_t>& orders, std::size_t& maximum_order)
{
    if (!std::isfinite(total_weight) || total_weight < 0.0)
    {
        return make_selected_failure(RawMomentStatus::invalid_input, "Total weight must be finite and non-negative.");
    }

    orders.clear();
    orders.reserve(request.indices.size());
    maximum_order = 0;
    for (const MultiIndex& index : request.indices)
    {
        if (index.size() != 1)
        {
            return make_selected_failure(RawMomentStatus::invalid_input, "Selected one-dimensional raw moment indices must have rank one.",
                                         invalid_dimension_index, index);
        }
        const std::size_t order = index[0];
        if (order == std::numeric_limits<std::size_t>::max())
        {
            return make_selected_failure(RawMomentStatus::invalid_input, "Selected raw moment order is too large.", 0, index);
        }
        orders.push_back(order);
        if (order > maximum_order)
        {
            maximum_order = order;
        }
    }
    return {};
}

[[nodiscard]] SelectedRawMomentResult finish_selected_orders(const RawMomentRequest& request, std::vector<double> values)
{
    return make_selected_raw_moment_result(request, std::move(values));
}

[[nodiscard]] SelectedRawMomentResult zero_selected_moments(const RawMomentRequest& request)
{
    return finish_selected_orders(request, std::vector<double>(request.indices.size(), 0.0));
}

[[nodiscard]] SelectedRawMomentResult selected_failure_from_moment(const MomentValueResult& result)
{
    return selected_failure_at_order(result.status, result.message, result.failed_order);
}

struct NormalExpectationsResult
{
    std::vector<double> expectations;
    RawMomentStatus status = RawMomentStatus::success;
    const char* message = "";
    std::size_t failed_order = 0;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == RawMomentStatus::success;
    }
};

[[nodiscard]] NormalExpectationsResult normal_expectations(const NormalDistribution& distribution, const std::size_t maximum_order)
{
    std::vector<double> expectations(maximum_order + 1, 0.0);
    expectations[0] = 1.0;
    if (maximum_order >= 1)
    {
        expectations[1] = distribution.mean;
    }
    // The variance enters from order 2 on, so it can only make a requested
    // moment unrepresentable when such an order is requested.
    if (maximum_order < 2)
    {
        return NormalExpectationsResult{std::move(expectations), RawMomentStatus::success, "", 0};
    }
    const double variance = distribution.standard_deviation * distribution.standard_deviation;
    if (!std::isfinite(variance))
    {
        return NormalExpectationsResult{{}, RawMomentStatus::overflow, "Normal variance overflows double precision.", 2};
    }

    for (std::size_t order = 2; order <= maximum_order; ++order)
    {
        // Multiply the (possibly exactly zero, for odd orders when mean == 0) prior
        // expectation by the small integer coefficient first, before multiplying by
        // variance. Reversing this order (variance first) can overflow the
        // intermediate product for a large but finite variance even though the
        // multiplication by an exact zero would otherwise make the true term zero,
        // spuriously turning `0 * inf` into NaN instead of the correct 0.0.
        expectations[order] = distribution.mean * expectations[order - 1] + variance * (expectations[order - 2] * static_cast<double>(order - 1));
        if (!std::isfinite(expectations[order]))
        {
            return NormalExpectationsResult{{}, finite_failure_status(expectations[order]), "Normal raw moment is not finite.", order};
        }
        // Even-order moments are strictly positive, and odd-order moments are
        // nonzero whenever the mean is, so an exact zero here is an underflow.
        if (expectations[order] == 0.0 && (order % 2 == 0 || distribution.mean != 0.0))
        {
            return NormalExpectationsResult{{}, RawMomentStatus::underflow, "Normal raw moment underflows double precision.", order};
        }
    }

    return NormalExpectationsResult{std::move(expectations), RawMomentStatus::success, "", 0};
}

[[nodiscard]] SelectedRawMomentResult compute_normal_selected(const NormalDistribution& distribution, const RawMomentRequest& request,
                                                              const std::vector<std::size_t>& orders, const std::size_t maximum_order,
                                                              const double total_weight)
{
    if (!std::isfinite(distribution.mean) || !is_positive_finite(distribution.standard_deviation))
    {
        return make_selected_failure(RawMomentStatus::invalid_input,
                                     "Normal distribution requires finite mean and positive finite standard deviation.");
    }
    if (total_weight == 0.0)
    {
        return zero_selected_moments(request);
    }

    NormalExpectationsResult expectations = normal_expectations(distribution, maximum_order);
    if (!expectations.ok())
    {
        return selected_failure_at_order(expectations.status, expectations.message, expectations.failed_order);
    }

    std::vector<double> values;
    values.reserve(orders.size());
    for (const std::size_t order : orders)
    {
        const double moment = total_weight * expectations.expectations[order];
        if (!std::isfinite(moment))
        {
            return selected_failure_at_order(finite_failure_status(moment), "Weighted normal raw moment is not finite.", order);
        }
        if (moment == 0.0 && expectations.expectations[order] != 0.0)
        {
            return selected_failure_at_order(RawMomentStatus::underflow, "Weighted normal raw moment underflows double precision.", order);
        }
        values.push_back(moment);
    }
    return finish_selected_orders(request, std::move(values));
}

[[nodiscard]] SelectedRawMomentResult compute_log_normal_selected(const LogNormalDistribution& distribution, const RawMomentRequest& request,
                                                                  const std::vector<std::size_t>& orders, const double total_weight)
{
    if (!std::isfinite(distribution.normal_mean) || !is_positive_finite(distribution.normal_standard_deviation))
    {
        return make_selected_failure(RawMomentStatus::invalid_input,
                                     "Log-normal distribution requires finite underlying mean and positive finite underlying standard deviation.");
    }

    const double variance = distribution.normal_standard_deviation * distribution.normal_standard_deviation;
    if (total_weight == 0.0)
    {
        return zero_selected_moments(request);
    }

    std::vector<double> values;
    values.reserve(orders.size());
    for (const std::size_t order : orders)
    {
        const MomentValueResult moment = evaluate_log_normal_moment(distribution, variance, total_weight, order);
        if (!moment.ok())
        {
            return selected_failure_from_moment(moment);
        }
        values.push_back(moment.value);
    }
    return finish_selected_orders(request, std::move(values));
}

[[nodiscard]] SelectedRawMomentResult compute_beta_selected(const BetaDistribution& distribution, const RawMomentRequest& request,
                                                            const std::vector<std::size_t>& orders, const double total_weight)
{
    if (!is_positive_finite(distribution.alpha) || !is_positive_finite(distribution.beta))
    {
        return make_selected_failure(RawMomentStatus::invalid_input, "Beta distribution requires positive finite alpha and beta.");
    }
    if (total_weight == 0.0)
    {
        return zero_selected_moments(request);
    }

    std::vector<double> values;
    values.reserve(orders.size());
    for (const std::size_t order : orders)
    {
        const MomentValueResult moment = evaluate_product_moment(total_weight, order, [&distribution](ScaledProduct& product, const std::size_t step)
                                                                 { multiply_beta_step(product, distribution, step); });
        if (!moment.ok())
        {
            return selected_failure_from_moment(moment);
        }
        values.push_back(moment.value);
    }
    return finish_selected_orders(request, std::move(values));
}

[[nodiscard]] SelectedRawMomentResult compute_gamma_selected(const GammaDistribution& distribution, const RawMomentRequest& request,
                                                             const std::vector<std::size_t>& orders, const double total_weight)
{
    if (!is_positive_finite(distribution.shape) || !is_positive_finite(distribution.scale))
    {
        return make_selected_failure(RawMomentStatus::invalid_input, "Gamma distribution requires positive finite shape and scale.");
    }
    if (total_weight == 0.0)
    {
        return zero_selected_moments(request);
    }

    std::vector<double> values;
    values.reserve(orders.size());
    for (const std::size_t order : orders)
    {
        const MomentValueResult moment = evaluate_product_moment(total_weight, order, [&distribution](ScaledProduct& product, const std::size_t step)
                                                                 { multiply_gamma_step(product, distribution, step); });
        if (!moment.ok())
        {
            return selected_failure_from_moment(moment);
        }
        values.push_back(moment.value);
    }
    return finish_selected_orders(request, std::move(values));
}

[[nodiscard]] SelectedRawMomentResult compute_chi_selected(const ChiDistribution& distribution, const RawMomentRequest& request,
                                                           const std::vector<std::size_t>& orders, const double total_weight)
{
    if (!is_positive_finite(distribution.degrees_of_freedom))
    {
        return make_selected_failure(RawMomentStatus::invalid_input, "Chi distribution requires positive finite degrees of freedom.");
    }
    if (total_weight == 0.0)
    {
        return zero_selected_moments(request);
    }

    const double first_moment = chi_first_moment(distribution);
    std::vector<double> values;
    values.reserve(orders.size());
    for (const std::size_t order : orders)
    {
        const MomentValueResult moment = evaluate_chi_moment(distribution, first_moment, total_weight, order);
        if (!moment.ok())
        {
            return selected_failure_from_moment(moment);
        }
        values.push_back(moment.value);
    }
    return finish_selected_orders(request, std::move(values));
}

[[nodiscard]] SelectedRawMomentResult compute_chi_squared_selected(const ChiSquaredDistribution& distribution, const RawMomentRequest& request,
                                                                   const std::vector<std::size_t>& orders, const double total_weight)
{
    if (!is_positive_finite(distribution.degrees_of_freedom))
    {
        return make_selected_failure(RawMomentStatus::invalid_input, "Chi-squared distribution requires positive finite degrees of freedom.");
    }
    if (total_weight == 0.0)
    {
        return zero_selected_moments(request);
    }

    std::vector<double> values;
    values.reserve(orders.size());
    for (const std::size_t order : orders)
    {
        const MomentValueResult moment = evaluate_product_moment(total_weight, order, [&distribution](ScaledProduct& product, const std::size_t step)
                                                                 { multiply_chi_squared_step(product, distribution, step); });
        if (!moment.ok())
        {
            return selected_failure_from_moment(moment);
        }
        values.push_back(moment.value);
    }
    return finish_selected_orders(request, std::move(values));
}

[[nodiscard]] SelectedRawMomentResult compute_weibull_selected(const WeibullDistribution& distribution, const RawMomentRequest& request,
                                                               const std::vector<std::size_t>& orders, const double total_weight)
{
    if (!is_positive_finite(distribution.shape) || !is_positive_finite(distribution.scale))
    {
        return make_selected_failure(RawMomentStatus::invalid_input, "Weibull distribution requires positive finite shape and scale.");
    }
    if (total_weight == 0.0)
    {
        return zero_selected_moments(request);
    }

    const double log_scale = std::log(distribution.scale);
    std::vector<double> values;
    values.reserve(orders.size());
    for (const std::size_t order : orders)
    {
        const MomentValueResult moment = evaluate_weibull_moment(distribution, log_scale, total_weight, order);
        if (!moment.ok())
        {
            return selected_failure_from_moment(moment);
        }
        values.push_back(moment.value);
    }
    return finish_selected_orders(request, std::move(values));
}

[[nodiscard]] RawMomentResult compute_normal(const NormalDistribution& distribution, const std::size_t maximum_order, const double total_weight)
{
    if (!std::isfinite(distribution.mean) || !is_positive_finite(distribution.standard_deviation))
    {
        return invalid_input("Normal distribution requires finite mean and positive finite standard deviation.");
    }
    if (total_weight == 0.0)
    {
        return finish_moments(std::vector<double>(maximum_order + 1, 0.0));
    }

    NormalExpectationsResult expectations = normal_expectations(distribution, maximum_order);
    if (!expectations.ok())
    {
        return failure_at_order(expectations.status, expectations.message, expectations.failed_order);
    }

    std::vector<double> moments(maximum_order + 1, 0.0);
    for (std::size_t order = 0; order <= maximum_order; ++order)
    {
        moments[order] = total_weight * expectations.expectations[order];
        if (!std::isfinite(moments[order]))
        {
            return failure_at_order(finite_failure_status(moments[order]), "Weighted normal raw moment is not finite.", order);
        }
        if (moments[order] == 0.0 && expectations.expectations[order] != 0.0)
        {
            return failure_at_order(RawMomentStatus::underflow, "Weighted normal raw moment underflows double precision.", order);
        }
    }
    return finish_moments(std::move(moments));
}

[[nodiscard]] RawMomentResult compute_log_normal(const LogNormalDistribution& distribution, const std::size_t maximum_order,
                                                 const double total_weight)
{
    if (!std::isfinite(distribution.normal_mean) || !is_positive_finite(distribution.normal_standard_deviation))
    {
        return invalid_input("Log-normal distribution requires finite underlying mean and positive finite underlying standard deviation.");
    }

    const double variance = distribution.normal_standard_deviation * distribution.normal_standard_deviation;
    return compute_from_log_moments(maximum_order, total_weight,
                                    [&](const std::size_t order) { return evaluate_log_normal_moment(distribution, variance, total_weight, order); });
}

[[nodiscard]] RawMomentResult compute_beta(const BetaDistribution& distribution, const std::size_t maximum_order, const double total_weight)
{
    if (!is_positive_finite(distribution.alpha) || !is_positive_finite(distribution.beta))
    {
        return invalid_input("Beta distribution requires positive finite alpha and beta.");
    }

    return compute_from_products(maximum_order, total_weight, [&distribution](ScaledProduct& product, const std::size_t order)
                                 { multiply_beta_step(product, distribution, order); });
}

[[nodiscard]] RawMomentResult compute_gamma(const GammaDistribution& distribution, const std::size_t maximum_order, const double total_weight)
{
    if (!is_positive_finite(distribution.shape) || !is_positive_finite(distribution.scale))
    {
        return invalid_input("Gamma distribution requires positive finite shape and scale.");
    }

    return compute_from_products(maximum_order, total_weight, [&distribution](ScaledProduct& product, const std::size_t order)
                                 { multiply_gamma_step(product, distribution, order); });
}

[[nodiscard]] RawMomentResult compute_chi(const ChiDistribution& distribution, const std::size_t maximum_order, const double total_weight)
{
    if (!is_positive_finite(distribution.degrees_of_freedom))
    {
        return invalid_input("Chi distribution requires positive finite degrees of freedom.");
    }

    if (total_weight == 0.0)
    {
        return finish_moments(std::vector<double>(maximum_order + 1, 0.0));
    }

    // Even and odd orders form separate chains, E[X^k] = E[X^(k-2)] (nu + k - 2),
    // each the same arithmetic evaluate_chi_moment performs for a single order.
    std::vector<double> moments(maximum_order + 1, 0.0);
    moments[0] = total_weight;
    ScaledProduct even(total_weight);
    ScaledProduct odd(total_weight);
    odd.multiply(chi_first_moment(distribution));
    for (std::size_t order = 1; order <= maximum_order; ++order)
    {
        ScaledProduct& chain = order % 2 == 0 ? even : odd;
        if (order >= 2)
        {
            multiply_chi_step(chain, distribution, order);
        }
        const MomentValueResult moment = chain.value(order);
        if (!moment.ok())
        {
            return failure_at_order(moment.status, moment.message, moment.failed_order);
        }
        moments[order] = moment.value;
    }
    return finish_moments(std::move(moments));
}

[[nodiscard]] RawMomentResult compute_chi_squared(const ChiSquaredDistribution& distribution, const std::size_t maximum_order,
                                                  const double total_weight)
{
    if (!is_positive_finite(distribution.degrees_of_freedom))
    {
        return invalid_input("Chi-squared distribution requires positive finite degrees of freedom.");
    }

    return compute_from_products(maximum_order, total_weight, [&distribution](ScaledProduct& product, const std::size_t order)
                                 { multiply_chi_squared_step(product, distribution, order); });
}

[[nodiscard]] RawMomentResult compute_weibull(const WeibullDistribution& distribution, const std::size_t maximum_order, const double total_weight)
{
    if (!is_positive_finite(distribution.shape) || !is_positive_finite(distribution.scale))
    {
        return invalid_input("Weibull distribution requires positive finite shape and scale.");
    }

    const double log_scale = std::log(distribution.scale);
    return compute_from_log_moments(maximum_order, total_weight,
                                    [&](const std::size_t order) { return evaluate_weibull_moment(distribution, log_scale, total_weight, order); });
}

struct RawMomentVisitor
{
    std::size_t maximum_order;
    double total_weight;

    [[nodiscard]] RawMomentResult operator()(const NormalDistribution& distribution) const
    {
        return compute_normal(distribution, maximum_order, total_weight);
    }

    [[nodiscard]] RawMomentResult operator()(const LogNormalDistribution& distribution) const
    {
        return compute_log_normal(distribution, maximum_order, total_weight);
    }

    [[nodiscard]] RawMomentResult operator()(const BetaDistribution& distribution) const
    {
        return compute_beta(distribution, maximum_order, total_weight);
    }

    [[nodiscard]] RawMomentResult operator()(const GammaDistribution& distribution) const
    {
        return compute_gamma(distribution, maximum_order, total_weight);
    }

    [[nodiscard]] RawMomentResult operator()(const ChiDistribution& distribution) const
    {
        return compute_chi(distribution, maximum_order, total_weight);
    }

    [[nodiscard]] RawMomentResult operator()(const ChiSquaredDistribution& distribution) const
    {
        return compute_chi_squared(distribution, maximum_order, total_weight);
    }

    [[nodiscard]] RawMomentResult operator()(const WeibullDistribution& distribution) const
    {
        return compute_weibull(distribution, maximum_order, total_weight);
    }
};

struct SelectedRawMomentVisitor
{
    const RawMomentRequest& request;
    const std::vector<std::size_t>& orders;
    std::size_t maximum_order;
    double total_weight;

    [[nodiscard]] SelectedRawMomentResult operator()(const NormalDistribution& distribution) const
    {
        return compute_normal_selected(distribution, request, orders, maximum_order, total_weight);
    }

    [[nodiscard]] SelectedRawMomentResult operator()(const LogNormalDistribution& distribution) const
    {
        return compute_log_normal_selected(distribution, request, orders, total_weight);
    }

    [[nodiscard]] SelectedRawMomentResult operator()(const BetaDistribution& distribution) const
    {
        return compute_beta_selected(distribution, request, orders, total_weight);
    }

    [[nodiscard]] SelectedRawMomentResult operator()(const GammaDistribution& distribution) const
    {
        return compute_gamma_selected(distribution, request, orders, total_weight);
    }

    [[nodiscard]] SelectedRawMomentResult operator()(const ChiDistribution& distribution) const
    {
        return compute_chi_selected(distribution, request, orders, total_weight);
    }

    [[nodiscard]] SelectedRawMomentResult operator()(const ChiSquaredDistribution& distribution) const
    {
        return compute_chi_squared_selected(distribution, request, orders, total_weight);
    }

    [[nodiscard]] SelectedRawMomentResult operator()(const WeibullDistribution& distribution) const
    {
        return compute_weibull_selected(distribution, request, orders, total_weight);
    }
};

[[nodiscard]] RawMomentResult add_compensated(double& sum, double& correction, const double value, const std::size_t order,
                                              const std::size_t component_index)
{
    if (!std::isfinite(value))
    {
        return make_failure(finite_failure_status(value), "Mixture component raw moment is not finite.", order, component_index);
    }

    const double updated_sum = sum + value;
    if (!std::isfinite(updated_sum))
    {
        return make_failure(finite_failure_status(updated_sum), "Mixture raw moment accumulation is not finite.", order, component_index);
    }

    if (std::abs(sum) >= std::abs(value))
    {
        correction += (sum - updated_sum) + value;
    }
    else
    {
        correction += (value - updated_sum) + sum;
    }
    if (!std::isfinite(correction))
    {
        return make_failure(finite_failure_status(correction), "Mixture compensation accumulation is not finite.", order, component_index);
    }

    sum = updated_sum;
    return {};
}

[[nodiscard]] SelectedRawMomentResult add_selected_compensated(double& sum, double& correction, const double value, const std::size_t order,
                                                               const std::size_t component_index)
{
    if (!std::isfinite(value))
    {
        return make_selected_component_failure(finite_failure_status(value), "Mixture component raw moment is not finite.", component_index, 0,
                                               {order});
    }

    const double updated_sum = sum + value;
    if (!std::isfinite(updated_sum))
    {
        return make_selected_component_failure(finite_failure_status(updated_sum), "Mixture raw moment accumulation is not finite.", component_index,
                                               0, {order});
    }

    if (std::abs(sum) >= std::abs(value))
    {
        correction += (sum - updated_sum) + value;
    }
    else
    {
        correction += (value - updated_sum) + sum;
    }
    if (!std::isfinite(correction))
    {
        return make_selected_component_failure(finite_failure_status(correction), "Mixture compensation accumulation is not finite.", component_index,
                                               0, {order});
    }

    sum = updated_sum;
    return {};
}

[[nodiscard]] bool contains_order(const std::vector<std::size_t>& orders, const std::size_t order)
{
    for (const std::size_t existing : orders)
    {
        if (existing == order)
        {
            return true;
        }
    }
    return false;
}

[[nodiscard]] std::size_t order_position(const std::vector<std::size_t>& orders, const std::size_t order)
{
    for (std::size_t index = 0; index < orders.size(); ++index)
    {
        if (orders[index] == order)
        {
            return index;
        }
    }
    throw std::out_of_range("Requested marginal moment order was not evaluated.");
}

[[nodiscard]] RawMomentResult compute_independent_product_marginal_moments(const IndependentProductMarginal& marginal,
                                                                           const std::size_t maximum_order)
{
    if (const auto* distribution = std::get_if<DistributionSpecification>(&marginal))
    {
        return compute_raw_moments(*distribution, maximum_order, 1.0);
    }
    return compute_raw_moments(std::get<DistributionMixture>(marginal), maximum_order);
}

[[nodiscard]] SelectedRawMomentResult compute_independent_product_marginal_moments(const IndependentProductMarginal& marginal,
                                                                                   const RawMomentRequest& request)
{
    if (const auto* distribution = std::get_if<DistributionSpecification>(&marginal))
    {
        return compute_raw_moments(*distribution, request, 1.0);
    }
    return compute_raw_moments(std::get<DistributionMixture>(marginal), request);
}

[[nodiscard]] MultiIndex first_requested_index_with_order(const RawMomentRequest& request, const std::size_t dimension, const std::size_t order)
{
    for (const MultiIndex& index : request.indices)
    {
        if (dimension < index.size() && index[dimension] == order)
        {
            return index;
        }
    }
    return request.indices.empty() ? MultiIndex{} : request.indices.front();
}

} // namespace

SelectedRawMomentResult make_selected_raw_moment_result(const RawMomentRequest& request, std::vector<double> values)
{
    if (request.indices.size() != values.size())
    {
        return make_selected_failure(RawMomentStatus::invalid_input, "Selected raw moment value count must match the request index count.");
    }

    if (request.indices.empty())
    {
        return SelectedRawMomentResult{{}, std::move(values), RawMomentStatus::success, "", invalid_component_index, invalid_dimension_index, {}};
    }

    const std::size_t rank = request.indices.front().size();
    if (rank == 0)
    {
        return make_selected_failure(RawMomentStatus::invalid_input, "Selected raw moment indices must have positive rank.");
    }

    for (std::size_t request_index = 0; request_index < request.indices.size(); ++request_index)
    {
        const MultiIndex& moment_order = request.indices[request_index];
        if (moment_order.size() != rank)
        {
            return make_selected_failure(RawMomentStatus::invalid_input, "Selected raw moment indices must have one common rank.",
                                         invalid_dimension_index, moment_order);
        }
        for (std::size_t dimension = 0; dimension < moment_order.size(); ++dimension)
        {
            if (moment_order[dimension] == std::numeric_limits<std::size_t>::max())
            {
                return make_selected_failure(RawMomentStatus::invalid_input, "Selected raw moment order is too large.", dimension, moment_order);
            }
        }
        if (!std::isfinite(values[request_index]))
        {
            return make_selected_failure(finite_failure_status(values[request_index]), "Selected raw moment value is not finite.",
                                         invalid_dimension_index, moment_order);
        }
    }

    return SelectedRawMomentResult{
        request.indices, std::move(values), RawMomentStatus::success, "", invalid_component_index, invalid_dimension_index, {}};
}

SelectedRawMomentResult compute_raw_moments(const DistributionSpecification& distribution, const RawMomentRequest& request, const double total_weight)
{
    std::vector<std::size_t> orders;
    std::size_t maximum_order = 0;
    SelectedRawMomentResult common = check_selected_common_inputs(request, total_weight, orders, maximum_order);
    if (!common.ok())
    {
        return common;
    }
    if (request.indices.empty())
    {
        return make_selected_raw_moment_result(request, {});
    }
    return std::visit(SelectedRawMomentVisitor{request, orders, maximum_order, total_weight}, distribution);
}

RawMomentResult compute_raw_moments(const DistributionSpecification& distribution, const std::size_t maximum_order, const double total_weight)
{
    RawMomentResult common = check_common_inputs(maximum_order, total_weight);
    if (!common.ok())
    {
        return common;
    }
    return std::visit(RawMomentVisitor{maximum_order, total_weight}, distribution);
}

RawMomentResult compute_raw_moments(const DistributionMixture& mixture, const std::size_t maximum_order)
{
    RawMomentResult maximum_order_check = check_maximum_order(maximum_order);
    if (!maximum_order_check.ok())
    {
        return maximum_order_check;
    }
    if (mixture.components.empty())
    {
        return invalid_input("Distribution mixture requires at least one component.");
    }

    bool has_positive_density = false;
    for (std::size_t component_index = 0; component_index < mixture.components.size(); ++component_index)
    {
        const double number_density = mixture.components[component_index].number_density;
        if (!std::isfinite(number_density) || number_density < 0.0)
        {
            return make_failure(RawMomentStatus::invalid_input, "Mixture component number density must be finite and non-negative.", 0,
                                component_index);
        }
        has_positive_density = has_positive_density || number_density > 0.0;
    }
    if (!has_positive_density)
    {
        return invalid_input("Distribution mixture requires at least one positive-density component.");
    }

    std::vector<double> sums(maximum_order + 1, 0.0);
    std::vector<double> corrections(maximum_order + 1, 0.0);

    for (std::size_t component_index = 0; component_index < mixture.components.size(); ++component_index)
    {
        const DistributionComponent& component = mixture.components[component_index];
        if (component.number_density == 0.0)
        {
            continue;
        }

        RawMomentResult component_result = compute_raw_moments(component.distribution, maximum_order, component.number_density);
        if (!component_result.ok())
        {
            return make_failure(component_result.status,
                                "Mixture component " + std::to_string(component_index) + " failed: " + component_result.message,
                                component_result.failed_order, component_index);
        }

        for (std::size_t order = 0; order <= maximum_order; ++order)
        {
            RawMomentResult accumulation_result =
                add_compensated(sums[order], corrections[order], component_result.moments[order], order, component_index);
            if (!accumulation_result.ok())
            {
                return accumulation_result;
            }
        }
    }

    std::vector<double> moments(maximum_order + 1, 0.0);
    for (std::size_t order = 0; order <= maximum_order; ++order)
    {
        moments[order] = sums[order] + corrections[order];
        if (!std::isfinite(moments[order]))
        {
            return make_failure(finite_failure_status(moments[order]), "Mixture corrected raw moment is not finite.", order);
        }
    }
    return finish_moments(std::move(moments));
}

SelectedRawMomentResult compute_raw_moments(const DistributionMixture& mixture, const RawMomentRequest& request)
{
    std::vector<std::size_t> orders;
    std::size_t maximum_order = 0;
    SelectedRawMomentResult request_check = check_selected_common_inputs(request, 1.0, orders, maximum_order);
    if (!request_check.ok())
    {
        return request_check;
    }
    if (mixture.components.empty())
    {
        return make_selected_failure(RawMomentStatus::invalid_input, "Distribution mixture requires at least one component.");
    }

    bool has_positive_density = false;
    for (std::size_t component_index = 0; component_index < mixture.components.size(); ++component_index)
    {
        const double number_density = mixture.components[component_index].number_density;
        if (!std::isfinite(number_density) || number_density < 0.0)
        {
            return make_selected_component_failure(RawMomentStatus::invalid_input,
                                                   "Mixture component number density must be finite and non-negative.", component_index);
        }
        has_positive_density = has_positive_density || number_density > 0.0;
    }
    if (!has_positive_density)
    {
        return make_selected_failure(RawMomentStatus::invalid_input, "Distribution mixture requires at least one positive-density component.");
    }

    if (request.indices.empty())
    {
        return make_selected_raw_moment_result(request, {});
    }

    std::vector<double> sums(request.indices.size(), 0.0);
    std::vector<double> corrections(request.indices.size(), 0.0);

    for (std::size_t component_index = 0; component_index < mixture.components.size(); ++component_index)
    {
        const DistributionComponent& component = mixture.components[component_index];
        if (component.number_density == 0.0)
        {
            continue;
        }

        SelectedRawMomentResult component_result = compute_raw_moments(component.distribution, request, component.number_density);
        if (!component_result.ok())
        {
            const std::size_t failed_dimension = component_result.failed_dimension == invalid_dimension_index ? 0 : component_result.failed_dimension;
            MultiIndex failed_order = component_result.failed_order.empty() ? MultiIndex{0} : component_result.failed_order;
            SelectedRawMomentResult result = make_selected_component_failure(
                component_result.status, "Mixture component " + std::to_string(component_index) + " failed: " + component_result.message,
                component_index, failed_dimension, std::move(failed_order));
            return result;
        }

        for (std::size_t request_index = 0; request_index < request.indices.size(); ++request_index)
        {
            SelectedRawMomentResult accumulation_result = add_selected_compensated(
                sums[request_index], corrections[request_index], component_result.values[request_index], orders[request_index], component_index);
            if (!accumulation_result.ok())
            {
                return accumulation_result;
            }
        }
    }

    std::vector<double> values(request.indices.size(), 0.0);
    for (std::size_t request_index = 0; request_index < request.indices.size(); ++request_index)
    {
        values[request_index] = sums[request_index] + corrections[request_index];
        if (!std::isfinite(values[request_index]))
        {
            return make_selected_failure(finite_failure_status(values[request_index]), "Mixture corrected raw moment is not finite.", 0,
                                         request.indices[request_index]);
        }
    }
    return make_selected_raw_moment_result(request, std::move(values));
}

TensorMomentResult compute_raw_moments(const IndependentProductDistribution& distribution, const std::vector<std::size_t>& maximum_orders)
{
    const std::size_t n_dimensions = distribution.marginals.size();
    if (n_dimensions == 0)
    {
        return invalid_tensor_input("Independent product distribution requires at least one marginal.");
    }
    if (maximum_orders.size() != n_dimensions)
    {
        return invalid_tensor_input("Maximum-order rank must match the number of marginals.");
    }
    if (!std::isfinite(distribution.number_density) || distribution.number_density < 0.0)
    {
        return invalid_tensor_input("Independent product number density must be finite and non-negative.");
    }

    std::vector<std::size_t> extents(n_dimensions, 0);
    for (std::size_t dimension = 0; dimension < n_dimensions; ++dimension)
    {
        if (maximum_orders[dimension] == std::numeric_limits<std::size_t>::max())
        {
            MultiIndex failed_order(n_dimensions, 0);
            failed_order[dimension] = maximum_orders[dimension];
            return make_tensor_failure(RawMomentStatus::invalid_input, "Maximum moment order is too large.", dimension, std::move(failed_order));
        }
        extents[dimension] = maximum_orders[dimension] + 1;
    }

    std::size_t n_tensor_values = 0;
    try
    {
        n_tensor_values = product_extents(extents);
    }
    catch (const std::exception& exception)
    {
        return invalid_tensor_input(std::string("The extents of the full tensor of raw moments are invalid: ") + exception.what());
    }

    if (distribution.number_density == 0.0)
    {
        return TensorMomentResult{RawMomentTensor{std::move(extents), std::vector<double>(n_tensor_values, 0.0)},
                                  RawMomentStatus::success,
                                  "",
                                  invalid_dimension_index,
                                  {}};
    }

    std::vector<MomentSet> marginal_moments;
    marginal_moments.reserve(n_dimensions);
    for (std::size_t dimension = 0; dimension < n_dimensions; ++dimension)
    {
        RawMomentResult marginal_result = compute_independent_product_marginal_moments(distribution.marginals[dimension], maximum_orders[dimension]);
        if (!marginal_result.ok())
        {
            MultiIndex failed_order(n_dimensions, 0);
            failed_order[dimension] = marginal_result.failed_order;
            return make_tensor_failure(marginal_result.status,
                                       "Independent product marginal " + std::to_string(dimension) + " failed: " + marginal_result.message, dimension,
                                       std::move(failed_order));
        }
        const double zeroth_moment = marginal_result.moments[0];
        if (!std::isfinite(zeroth_moment) || zeroth_moment <= 0.0)
        {
            MultiIndex failed_order(n_dimensions, 0);
            return make_tensor_failure(RawMomentStatus::invalid_input, "Independent product marginal zeroth raw moment must be finite and positive.",
                                       dimension, std::move(failed_order));
        }
        marginal_moments.push_back(std::move(marginal_result.moments));
    }

    std::vector<double> values(n_tensor_values, 0.0);
    for (std::size_t flat_index = 0; flat_index < n_tensor_values; ++flat_index)
    {
        const MultiIndex moment_order = unflatten_multi_index(flat_index, extents);
        double value = distribution.number_density;
        for (std::size_t dimension = 0; dimension < n_dimensions; ++dimension)
        {
            const MomentValueResult normalized = normalize_marginal_moment(marginal_moments[dimension][moment_order[dimension]],
                                                                           marginal_moments[dimension][0], moment_order[dimension]);
            if (!normalized.ok())
            {
                return make_tensor_failure(normalized.status, normalized.message, dimension, moment_order);
            }
            const double factor = normalized.value;
            const double product = value * factor;
            if (!std::isfinite(product))
            {
                return make_tensor_failure(finite_failure_status(product), "Independent product raw moment is not finite.", invalid_dimension_index,
                                           moment_order);
            }
            if (value != 0.0 && factor != 0.0 && product == 0.0)
            {
                return make_tensor_failure(RawMomentStatus::underflow, "Independent product raw moment underflows double precision.",
                                           invalid_dimension_index, moment_order);
            }
            value = product;
        }
        values[flat_index] = value;
    }

    return TensorMomentResult{RawMomentTensor{std::move(extents), std::move(values)}, RawMomentStatus::success, "", invalid_dimension_index, {}};
}

SelectedRawMomentResult compute_raw_moments(const IndependentProductDistribution& distribution, const RawMomentRequest& request)
{
    const std::size_t n_dimensions = distribution.marginals.size();
    if (n_dimensions == 0)
    {
        return make_selected_failure(RawMomentStatus::invalid_input, "Independent product distribution requires at least one marginal.");
    }
    if (!std::isfinite(distribution.number_density) || distribution.number_density < 0.0)
    {
        return make_selected_failure(RawMomentStatus::invalid_input, "Independent product number density must be finite and non-negative.");
    }

    std::vector<std::vector<std::size_t>> marginal_orders(n_dimensions);
    for (const MultiIndex& index : request.indices)
    {
        if (index.size() != n_dimensions)
        {
            return make_selected_failure(RawMomentStatus::invalid_input,
                                         "Selected independent product raw moment index rank must match the number of marginals.",
                                         invalid_dimension_index, index);
        }
        for (std::size_t dimension = 0; dimension < n_dimensions; ++dimension)
        {
            if (index[dimension] == std::numeric_limits<std::size_t>::max())
            {
                return make_selected_failure(RawMomentStatus::invalid_input, "Selected raw moment order is too large.", dimension, index);
            }
            if (!contains_order(marginal_orders[dimension], index[dimension]))
            {
                marginal_orders[dimension].push_back(index[dimension]);
            }
        }
    }

    if (request.indices.empty())
    {
        return make_selected_raw_moment_result(request, {});
    }
    if (distribution.number_density == 0.0)
    {
        return make_selected_raw_moment_result(request, std::vector<double>(request.indices.size(), 0.0));
    }

    for (std::vector<std::size_t>& orders : marginal_orders)
    {
        if (!contains_order(orders, 0))
        {
            orders.push_back(0);
        }
    }

    std::vector<std::vector<double>> marginal_values(n_dimensions);
    for (std::size_t dimension = 0; dimension < n_dimensions; ++dimension)
    {
        RawMomentRequest marginal_request;
        marginal_request.indices.reserve(marginal_orders[dimension].size());
        for (const std::size_t order : marginal_orders[dimension])
        {
            marginal_request.indices.push_back({order});
        }

        SelectedRawMomentResult marginal_result = compute_independent_product_marginal_moments(distribution.marginals[dimension], marginal_request);
        if (!marginal_result.ok())
        {
            const std::size_t marginal_order = marginal_result.failed_order.empty() ? 0 : marginal_result.failed_order.front();
            const MultiIndex failed_order = first_requested_index_with_order(request, dimension, marginal_order);
            SelectedRawMomentResult result = make_selected_failure(
                marginal_result.status, "Independent product marginal " + std::to_string(dimension) + " failed: " + marginal_result.message,
                dimension, failed_order);
            result.failed_component = marginal_result.failed_component;
            return result;
        }
        const std::size_t zeroth_position = order_position(marginal_orders[dimension], 0);
        const double zeroth_moment = marginal_result.values[zeroth_position];
        if (!std::isfinite(zeroth_moment) || zeroth_moment <= 0.0)
        {
            return make_selected_failure(RawMomentStatus::invalid_input,
                                         "Independent product marginal zeroth raw moment must be finite and positive.", dimension,
                                         first_requested_index_with_order(request, dimension, 0));
        }
        marginal_values[dimension] = std::move(marginal_result.values);
    }

    std::vector<double> values;
    values.reserve(request.indices.size());
    for (const MultiIndex& moment_order : request.indices)
    {
        double value = distribution.number_density;
        for (std::size_t dimension = 0; dimension < n_dimensions; ++dimension)
        {
            const std::size_t position = order_position(marginal_orders[dimension], moment_order[dimension]);
            const std::size_t zeroth_position = order_position(marginal_orders[dimension], 0);
            const MomentValueResult normalized =
                normalize_marginal_moment(marginal_values[dimension][position], marginal_values[dimension][zeroth_position], moment_order[dimension]);
            if (!normalized.ok())
            {
                return make_selected_failure(normalized.status, normalized.message, dimension, moment_order);
            }
            const double factor = normalized.value;
            const double product = value * factor;
            if (!std::isfinite(product))
            {
                return make_selected_failure(finite_failure_status(product), "Independent product raw moment is not finite.", invalid_dimension_index,
                                             moment_order);
            }
            if (value != 0.0 && factor != 0.0 && product == 0.0)
            {
                return make_selected_failure(RawMomentStatus::underflow, "Independent product raw moment underflows double precision.",
                                             invalid_dimension_index, moment_order);
            }
            value = product;
        }
        values.push_back(value);
    }

    return make_selected_raw_moment_result(request, std::move(values));
}

} // namespace qmomkit
