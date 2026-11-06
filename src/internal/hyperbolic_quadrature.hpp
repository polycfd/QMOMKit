#pragma once

#include "qmomkit/moment_set.hpp"
#include "qmomkit/quadrature.hpp"
#include "qmomkit/wheeler.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace qmomkit::internal
{

/**
 * HyQMOM reconstruction (Fox & Laurent 2022) from m_0, ..., m_{2n}, n =
 * node_count, on the real line: the n + 1 nodes of J_{n+1} when the rounding
 * test finds the moments inside the realizable set. Otherwise adaptive mode
 * returns the classical quadrature that compute_quadrature(...) returns for n
 * nodes, and fixed mode is not_realizable. rmin and eabs apply only to that
 * classical quadrature. The coordinate scale is zero, as in compute_quadrature(...).
 */
[[nodiscard]] Quadrature1D compute_hyperbolic_quadrature(const MomentSet& moments, std::size_t node_count, const WheelerOptions& options);

struct HyperbolicCharacteristicSpeeds
{
    std::vector<double> speeds;
    QuadratureStatus status = QuadratureStatus::success;
    std::string message;
};

/**
 * The characteristic speeds, in ascending order, of the moment system of free
 * transport closed by HyQMOM with node count n = node_count at the moments
 * m_0, ..., m_{2n}: those of the quadrature that compute_hyperbolic_quadrature
 * reconstructs from the same moments and options.
 * - n + 1 nodes: the 2n + 1 eigenvalues of J_n and K_{n+1}. K_{n+1} adds to
 *   J_n the closure a_n of compute_hyperbolic_quadrature on the diagonal and
 *   the square root of beta_n = (2n + 1) b_n / n beside it (Fox & Laurent
 *   2022, eqs. 4.17 and 4.18).
 * - the classical quadrature of k <= n nodes of adaptive mode: the 2k + 1
 *   speeds of the k-point system on its boundary, beta_k = 0, its k
 *   abscissas, each twice, and their arithmetic mean. For k = n, these are the
 *   boundary speeds of Fox & Laurent 2022, Remark 4.9. For k = 1, all three
 *   speeds are the abscissa.
 *
 * Where the reconstruction fails, its status and message are returned with no
 * speeds. The status is also not_realizable, with no speeds, when an
 * eigenvalue computation does not converge or a speed is too large for a
 * double.
 */
[[nodiscard]] HyperbolicCharacteristicSpeeds compute_hyperbolic_characteristic_speeds(const MomentSet& moments, std::size_t node_count,
                                                                                      const WheelerOptions& options);

} // namespace qmomkit::internal
