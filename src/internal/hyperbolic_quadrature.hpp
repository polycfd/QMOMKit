#pragma once

#include "qmomkit/moment_set.hpp"
#include "qmomkit/quadrature.hpp"
#include "qmomkit/wheeler.hpp"

#include <cstddef>

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

} // namespace qmomkit::internal
