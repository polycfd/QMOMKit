#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/distributions.hpp"
#include "qmomkit/moment_set.hpp"
#include "qmomkit/raw_moment_tensor.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace qmomkit::monte_carlo
{

/**
 * Status for estimating raw moments from a finite sample ensemble.
 *
 * Success means all requested estimates and standard errors are finite. It
 * does not mean the estimates equal analytical distribution moments exactly.
 */
enum class EnsembleMomentStatus
{
    success,
    invalid_input,
    non_finite_sample,
    non_finite_moment
};

struct OneDimensionalMomentEstimate
{
    // Raw moments in QMOMKit order: moments[k] == m_k.
    qmomkit::MomentSet moments;
    // Estimated standard error for each moment estimate.
    std::vector<double> standard_errors;
    EnsembleMomentStatus status = EnsembleMomentStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == EnsembleMomentStatus::success;
    }
};

struct TensorMomentEstimate
{
    // Raw tensor moments in phase-space-order, row-major QMOMKit storage.
    qmomkit::RawMomentTensor moments;
    // Estimated standard error for each flat tensor-moment entry.
    std::vector<double> standard_errors;
    EnsembleMomentStatus status = EnsembleMomentStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == EnsembleMomentStatus::success;
    }
};

/**
 * Selected multidimensional raw moments estimated from an equally weighted
 * finite ensemble.
 *
 * On success, `indices`, `values`, and `standard_errors` are aligned one-to-one
 * in the order supplied by the request. Duplicate requested indices remain
 * duplicate result entries.
 */
struct SelectedEnsembleMomentEstimate
{
    // Requested raw-moment indices in the original request order.
    std::vector<qmomkit::MultiIndex> indices;
    // Estimated raw moment corresponding to each requested index.
    std::vector<double> values;
    // Estimated standard error corresponding to each requested index.
    std::vector<double> standard_errors;
    EnsembleMomentStatus status = EnsembleMomentStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == EnsembleMomentStatus::success;
    }
};

/**
 * Generalized raw moments estimated from an equally weighted finite ensemble.
 *
 * `orders` and `values` remain aligned in the input request order. The
 * exponent vectors use physical phase-space dimension order.
 */
struct GeneralizedEnsembleMomentEstimate
{
    std::vector<std::vector<double>> orders;
    std::vector<double> values;
    EnsembleMomentStatus status = EnsembleMomentStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == EnsembleMomentStatus::success;
    }
};

/**
 * Estimate one-dimensional raw moments from equally weighted samples.
 *
 * For `N` one-dimensional samples `x_p` and total number density `n0`, the
 * estimator is:
 *
 *     m_k = n0 / N * sum_p x_p^k
 *
 * Thus `m_0 == n0` for a valid non-empty ensemble. The returned standard error
 * estimates finite-sample uncertainty of the scaled sample mean for each
 * monomial. Increasing sample count usually reduces that uncertainty, while
 * high-order moments may converge slowly because they depend strongly on tail
 * samples.
 */
[[nodiscard]] OneDimensionalMomentEstimate compute_ensemble_moments(const Ensemble& ensemble, double number_density, std::size_t maximum_order);

/**
 * Estimate multidimensional raw tensor moments from equally weighted samples.
 *
 * For multi-index `alpha`, samples `xi_p`, and total number density `n0`, the
 * estimator is:
 *
 *     M_alpha = n0 / N * sum_p product_d xi[p,d]^alpha[d]
 *
 * `maximum_orders[d]` selects the largest moment order in phase-space dimension
 * `d`; tensor extents are `maximum_orders[d] + 1`. Values and standard errors
 * use QMOMKit's phase-space-order, row-major flat tensor layout. Coordinate
 * names and physical interpretations do not affect the ordering.
 */
[[nodiscard]] TensorMomentEstimate compute_ensemble_moments(const Ensemble& ensemble, double number_density,
                                                            const std::vector<std::size_t>& maximum_orders);

/**
 * Estimate only the requested multidimensional raw moments from equally
 * weighted samples.
 *
 * The total number density is applied using the same convention as the tensor
 * estimator:
 *
 *     M_alpha = n0 / N * sum_p product_d xi[p,d]^alpha[d]
 *
 * Request order and duplicate entries are preserved. The function evaluates
 * only the requested monomials and does not construct the full moment tensor.
 */
[[nodiscard]] SelectedEnsembleMomentEstimate compute_requested_ensemble_moments(const Ensemble& ensemble, double number_density,
                                                                                const qmomkit::RawMomentRequest& request);

/**
 * Estimate requested generalized raw moments from equally weighted samples.
 *
 * For each exponent vector `alpha`, the estimator is:
 *
 *     M_alpha = n0 / N * sum_p product_d xi[p,d]^alpha[d]
 *
 * Exponents must be finite, non-negative, and have the sample phase-space
 * rank. A negative coordinate with a non-integer exponent is invalid.
 */
[[nodiscard]] GeneralizedEnsembleMomentEstimate compute_generalized_ensemble_moments(const Ensemble& ensemble, double number_density,
                                                                                     const std::vector<std::vector<double>>& orders);

} // namespace qmomkit::monte_carlo
