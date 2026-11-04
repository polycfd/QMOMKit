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
 * test finds the moments in the interior of moment space. Otherwise adaptive
 * mode returns the classical rule that compute_quadrature(...) returns for n
 * nodes, and fixed mode is not_realizable. rmin and eabs apply only to that
 * classical rule. The coordinate scale is zero, as in compute_quadrature(...).
 */
[[nodiscard]] Quadrature1D compute_hyperbolic_quadrature(const MomentSet& moments, std::size_t node_count, const WheelerOptions& options);

struct HyperbolicCharacteristicSpeeds
{
    std::vector<double> speeds;
    QuadratureStatus status = QuadratureStatus::success;
    std::string message;
};

/**
 * The 2n + 1 characteristic speeds, ascending, of the moment system of free
 * transport closed by HyQMOM at the moments m_0, ..., m_{2n} of n + 1 weighted
 * points, n + 1 = abscissas.size(): the eigenvalues of J_n and of K_{n+1}
 * (Fox & Laurent 2022, eqs. 4.17 and 4.18), from the recurrence of the points
 * and the closure of compute_hyperbolic_quadrature. The points must be at
 * least two, finite, with positive weights, as the nodes of a hyperbolic
 * state at its maximum count are. not_realizable, with no speeds, when the
 * recurrence of the points vanishes, b_k = 0, or is not finite, as for points
 * all at zero or coordinates beyond the range of a double, or when an
 * eigensolve fails.
 */
[[nodiscard]] HyperbolicCharacteristicSpeeds compute_hyperbolic_characteristic_speeds(const std::vector<double>& abscissas,
                                                                                      const std::vector<double>& weights);

} // namespace qmomkit::internal
