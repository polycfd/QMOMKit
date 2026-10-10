#pragma once

#include "qmomkit/moment_set.hpp"
#include "qmomkit/multi_index.hpp"
#include "qmomkit/raw_moment_tensor.hpp"

#include <cstddef>
#include <string>
#include <variant>
#include <vector>

namespace qmomkit
{

/**
 * Normal distribution parameterized by mean and standard deviation.
 */
struct NormalDistribution
{
    double mean = 0.0;
    double standard_deviation = 1.0;
};

/**
 * Log-normal distribution parameterized by the underlying normal variable.
 */
struct LogNormalDistribution
{
    double normal_mean = 0.0;
    double normal_standard_deviation = 1.0;
};

/**
 * Beta distribution on [0, 1] parameterized by shape parameters alpha and beta.
 */
struct BetaDistribution
{
    double alpha = 1.0;
    double beta = 1.0;
};

/**
 * Gamma distribution parameterized by shape and scale.
 */
struct GammaDistribution
{
    double shape = 1.0;
    double scale = 1.0;
};

/**
 * Chi distribution parameterized by degrees of freedom.
 */
struct ChiDistribution
{
    double degrees_of_freedom = 1.0;
};

/**
 * Chi-squared distribution parameterized by degrees of freedom.
 */
struct ChiSquaredDistribution
{
    double degrees_of_freedom = 1.0;
};

/**
 * Weibull distribution parameterized by shape and scale.
 */
struct WeibullDistribution
{
    double shape = 1.0;
    double scale = 1.0;
};

using DistributionSpecification = std::variant<NormalDistribution, LogNormalDistribution, BetaDistribution, GammaDistribution, ChiDistribution,
                                               ChiSquaredDistribution, WeibullDistribution>;

struct DistributionComponent
{
    DistributionSpecification distribution;
    double number_density = 1.0;
};

struct DistributionMixture
{
    std::vector<DistributionComponent> components;
};

using IndependentProductMarginal = std::variant<DistributionSpecification, DistributionMixture>;

struct IndependentProductDistribution
{
    std::vector<IndependentProductMarginal> marginals;
    double number_density = 1.0;
};

enum class RawMomentStatus
{
    success,
    invalid_input,
    overflow,
    underflow,
    non_finite
};

inline constexpr std::size_t invalid_component_index = static_cast<std::size_t>(-1);
inline constexpr std::size_t invalid_dimension_index = static_cast<std::size_t>(-1);

struct RawMomentResult
{
    MomentSet moments;
    RawMomentStatus status = RawMomentStatus::success;
    std::string message;
    std::size_t failed_order = 0;
    std::size_t failed_component = invalid_component_index;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == RawMomentStatus::success;
    }
};

/**
 * Request selected raw moments by phase-space-order moment multi-index.
 *
 * Request order is significant. A result for this request stores values in the
 * same order, one value per requested index. Duplicate indices are allowed and
 * produce duplicate result entries. An empty request is valid and produces an
 * empty successful result when evaluated.
 */
struct RawMomentRequest
{
    std::vector<MultiIndex> indices;
};

/**
 * Selected raw moment values aligned with a `RawMomentRequest`.
 *
 * `values[i]` stores the raw moment for `indices[i]`. The container is a
 * request/result pair, not a tensor representation: indices are not sorted,
 * deduplicated, or reshaped.
 */
struct SelectedRawMomentResult
{
    std::vector<MultiIndex> indices;
    std::vector<double> values;
    RawMomentStatus status = RawMomentStatus::success;
    std::string message;
    std::size_t failed_component = invalid_component_index;
    std::size_t failed_dimension = invalid_dimension_index;
    MultiIndex failed_order;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == RawMomentStatus::success;
    }
};

struct TensorMomentResult
{
    RawMomentTensor moments;
    RawMomentStatus status = RawMomentStatus::success;
    std::string message;
    std::size_t failed_dimension = invalid_dimension_index;
    MultiIndex failed_order;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == RawMomentStatus::success;
    }
};

/**
 * Build a selected raw moment result after validating request/result shape.
 *
 * A valid request is empty, or contains only non-empty multi-indices of one
 * common rank. The maximum `std::size_t` value is rejected as a moment order,
 * as in rectangular moment generation, where it would overflow an extent.
 * Values must be finite and have the same count as requested indices.
 */
[[nodiscard]] SelectedRawMomentResult make_selected_raw_moment_result(const RawMomentRequest& request, std::vector<double> values);

/**
 * Generate selected one-dimensional raw moments for an analytical distribution.
 *
 * Each requested index must have rank one, with `index[0]` storing the moment
 * order. Request order and duplicate entries are preserved in the result.
 */
[[nodiscard]] SelectedRawMomentResult compute_raw_moments(const DistributionSpecification& distribution, const RawMomentRequest& request,
                                                          double total_weight = 1.0);

/**
 * Generate one-dimensional raw moments m_0 through m_maximum_order.
 *
 * The returned sequence is unnormalized unless total_weight is one:
 *
 *     m_k = total_weight * E[xi^k]
 *
 * Therefore m_0 is always total_weight on success.
 */
[[nodiscard]] RawMomentResult compute_raw_moments(const DistributionSpecification& distribution, std::size_t maximum_order,
                                                  double total_weight = 1.0);

/**
 * Generate one-dimensional raw moments for a heterogeneous mixture.
 *
 * Each component uses an absolute number density n_r:
 *
 *     m_k = sum_r n_r * E_r[xi^k]
 *
 * Component densities are not normalized internally. Therefore m_0 is the sum
 * of component densities on success.
 */
[[nodiscard]] RawMomentResult compute_raw_moments(const DistributionMixture& mixture, std::size_t maximum_order);

/**
 * Generate selected one-dimensional raw moments for a heterogeneous mixture.
 *
 * Request order and duplicate entries are preserved. Component number
 * densities are absolute NDF contributions, matching the maximum-order mixture
 * API.
 */
[[nodiscard]] SelectedRawMomentResult compute_raw_moments(const DistributionMixture& mixture, const RawMomentRequest& request);

/**
 * Generate a raw moment tensor for an independent multidimensional product.
 *
 * Marginals are normalized through their zeroth raw moments. The global number
 * density is applied exactly once:
 *
 *     M[k_0, ..., k_D-1] = number_density * product_d E[xi_d^k_d]
 *
 * The requested maximum orders are inclusive, so
 * `extents[d] = maximum_orders[d] + 1`.
 */
[[nodiscard]] TensorMomentResult compute_raw_moments(const IndependentProductDistribution& distribution,
                                                     const std::vector<std::size_t>& maximum_orders);

/**
 * Generate selected raw moments for an independent multidimensional product.
 *
 * Requested multi-indices are interpreted in phase-space dimension order.
 * Request order and duplicate entries are preserved. The global number density
 * is applied exactly once.
 */
[[nodiscard]] SelectedRawMomentResult compute_raw_moments(const IndependentProductDistribution& distribution, const RawMomentRequest& request);

} // namespace qmomkit
