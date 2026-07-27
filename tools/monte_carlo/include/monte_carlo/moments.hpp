#pragma once

#include "qmomkit/conditional_moments.hpp"
#include "qmomkit/moment_set.hpp"
#include "monte_carlo/ensemble.hpp"

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

} // namespace qmomkit::monte_carlo
