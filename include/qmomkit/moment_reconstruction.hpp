#pragma once

#include "qmomkit/cqmom.hpp"
#include "qmomkit/distributions.hpp"
#include "qmomkit/quadrature.hpp"

#include <cstddef>
#include <vector>

namespace qmomkit
{

/**
 * Reconstruct one-dimensional raw moments from a successful QMOM quadrature
 * rule using the convention m_k = sum_i w_i x_i^k.
 */
[[nodiscard]] std::vector<double> reconstruct_raw_moments(const Quadrature1D& quadrature, std::size_t moment_count);

/**
 * Reconstruct multidimensional raw moments from terminal CQMOM quadrature
 * nodes using phase-space-dimension moment ordering.
 */
[[nodiscard]] RawMomentTensor reconstruct_raw_moments(const CqmomQuadrature& quadrature, const std::vector<std::size_t>& extents);

/**
 * Reconstruct selected multidimensional raw moments from terminal CQMOM
 * quadrature nodes using the exact request order and phase-space indices.
 */
[[nodiscard]] SelectedRawMomentResult reconstruct_raw_moments(const CqmomQuadrature& quadrature, const RawMomentRequest& request);

/**
 * Reconstruct only the raw moments required by a fixed CQMOM moment layout
 * from terminal CQMOM quadrature nodes.
 *
 * The returned values are aligned with `layout.required_indices`, whose
 * multi-indices are ordered by phase-space dimension. Fixed rectangular and
 * adaptive ragged quadrature layouts are both reconstructed from terminal
 * nodes into the fixed required moment layout.
 */
[[nodiscard]] CqmomMomentData reconstruct_cqmom_moments(const CqmomQuadrature& quadrature, const CqmomMomentLayout& layout);

/**
 * Evaluate a one-dimensional generalized moment from an existing quadrature
 * rule using M_alpha = sum_i w_i x_i^alpha.
 *
 * The exponent may be non-integer. The exponent, weights, and abscissas must be
 * finite; invalid storage or non-finite exponents are rejected with
 * `std::invalid_argument`. Negative abscissas with non-integer exponents raise
 * `std::domain_error`. Non-finite evaluated products or sums raise
 * `std::overflow_error`.
 *
 * Generalized moments are quadrature-derived post-processing quantities. This
 * function does not create or update a generalized moment state, and the
 * returned value is not a QMOM reconstruction or transport input.
 */
[[nodiscard]] double evaluate_generalized_moment(const Quadrature1D& quadrature, double exponent);

/**
 * Evaluate a multidimensional generalized moment from terminal CQMOM
 * quadrature nodes using M_alpha = sum_j W_j product_d x_(d,j)^alpha_d.
 *
 * The exponent vector is ordered by phase-space dimension and may contain
 * non-integer values. Exponents, weights, and abscissas must be finite;
 * invalid storage or non-finite exponents are rejected with
 * `std::invalid_argument`. Negative abscissas with non-integer exponents raise
 * `std::domain_error`. Non-finite evaluated products or sums raise
 * `std::overflow_error`.
 *
 * Generalized moments are quadrature-derived post-processing quantities. This
 * function does not create or update a generalized moment state, and the
 * returned value is not a CQMOM reconstruction or transport input.
 */
[[nodiscard]] double evaluate_generalized_moment(const CqmomQuadrature& quadrature, const std::vector<double>& exponents);

} // namespace qmomkit
