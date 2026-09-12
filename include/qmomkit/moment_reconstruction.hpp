#pragma once

#include "qmomkit/cqmom.hpp"
#include "qmomkit/distributions.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/quadrature.hpp"

#include <cstddef>
#include <vector>

namespace qmomkit
{

/**
 * Advanced composable operation: reconstruct one-dimensional raw moments from
 * a successful QMOM quadrature rule using m_k = sum_i w_i x_i^k.
 *
 * An unsuccessful quadrature, a zero `moment_count`, or inconsistent/non-finite
 * quadrature storage is rejected with `std::invalid_argument`. A non-finite
 * reconstructed value raises `std::overflow_error`.
 */
[[nodiscard]] std::vector<double> reconstruct_raw_moments(const Quadrature1D& quadrature, std::size_t moment_count);

/**
 * Advanced composable operation: reconstruct multidimensional raw moments
 * from terminal CQMOM quadrature nodes using phase-space-dimension ordering.
 *
 * A rank mismatch between `extents` and the quadrature's phase-space dimension
 * count, or inconsistent/non-finite quadrature storage, is rejected with
 * `std::invalid_argument`. A non-finite reconstructed value raises
 * `std::overflow_error`.
 */
[[nodiscard]] RawMomentTensor reconstruct_raw_moments(const CqmomQuadrature& quadrature, const std::vector<std::size_t>& extents);

/**
 * Advanced composable operation: reconstruct selected multidimensional raw
 * moments from terminal CQMOM quadrature nodes using the exact request order
 * and phase-space indices.
 *
 * Inconsistent/non-finite quadrature storage, or a requested index whose rank
 * does not match the quadrature's phase-space dimension count, is rejected
 * with `std::invalid_argument`. A non-finite reconstructed value raises
 * `std::overflow_error`.
 */
[[nodiscard]] SelectedRawMomentResult reconstruct_raw_moments(const CqmomQuadrature& quadrature, const RawMomentRequest& request);

/**
 * Advanced composable operation: reconstruct only the raw moments required by
 * a fixed CQMOM moment layout from terminal CQMOM quadrature nodes.
 *
 * The returned values are aligned with `layout.required_indices`, whose
 * multi-indices are ordered by phase-space dimension. Fixed rectangular and
 * adaptive ragged quadrature layouts are both reconstructed from terminal
 * nodes into the fixed required moment layout.
 *
 * A layout hierarchy that does not match the quadrature hierarchy, an
 * inconsistent/non-finite quadrature, or a layout that is not the fixed
 * required layout for its hierarchy is rejected with `std::invalid_argument`.
 * A required index outside the hierarchy's moment extents raises
 * `std::out_of_range`. A non-finite reconstructed value raises
 * `std::overflow_error`.
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
 * returned value is not a QMOM reconstruction or phase-space evolution input.
 */
[[nodiscard]] double evaluate_generalized_moment(const Quadrature1D& quadrature, double exponent);

/** Evaluate a one-dimensional generalized moment from a public QMOM result. */
[[nodiscard]] double evaluate_generalized_moment(const QmomResult& quadrature, double exponent);

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
 * returned value is not a CQMOM reconstruction or phase-space evolution input.
 */
[[nodiscard]] double evaluate_generalized_moment(const CqmomQuadrature& quadrature, const std::vector<double>& exponents);

} // namespace qmomkit
