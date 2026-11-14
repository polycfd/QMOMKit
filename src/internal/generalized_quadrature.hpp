#pragma once

#include "qmomkit/abscissa_bounds.hpp"
#include "qmomkit/moment_set.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/quadrature.hpp"
#include "qmomkit/wheeler.hpp"

#include <cstddef>

namespace qmomkit::internal
{

struct GeneralizedQuadrature
{
    Quadrature1D quadrature;
    // The number k of nodes of the classical quadrature on which the result
    // rests, 0 on failure.
    std::size_t base_node_count = 0;
};

/**
 * GQMOM reconstruction (Fox, Laurent & Passalacqua 2023) from m_0, ..., m_{2N},
 * N = node_count, with the N_G = generalized.node_count nodes of the continued
 * recurrence. Only the normal continuation reads nu. The gamma continuation,
 * which needs its bounds, goes through compute_bounded_generalized_quadrature(...)
 * instead. Adaptive mode takes the classical selection of
 * compute_quadrature(...) for N nodes, with rmin and eabs, and fixed mode the
 * classical N-node checks. The k-node quadrature found is continued to N_G nodes
 * when sigma(k, k) exceeds its rounding bound tau_k, and returned as it is
 * when sigma(k, k) is zero within tau_k. Below -tau_k, fixed mode is
 * not_realizable and adaptive mode returns the classical quadrature. A
 * non-finite sigma(k, k) or tau_k is not_realizable in fixed mode. In adaptive
 * mode, an overflowed +-inf returns the classical quadrature and a NaN is
 * not_realizable. The coordinate scale is zero, as in compute_quadrature(...).
 * With exact_base_node_count, a base node count below N is not continued, and
 * the result is not_realizable instead.
 */
[[nodiscard]] GeneralizedQuadrature compute_generalized_quadrature(const MomentSet& moments, std::size_t node_count, const WheelerOptions& options,
                                                                   const GeneralizedOptions& generalized, bool exact_base_node_count);

/**
 * compute_generalized_quadrature(...) within abscissa bounds, which give the
 * gamma continuation its lower bound. As for bounded classical QMOM, a
 * reconstruction is accepted when it succeeds with every node inside the
 * bounds. Fixed mode is not_realizable otherwise. Adaptive mode retries below
 * the base node count of a rejected reconstruction, down to one node, and
 * returns the first reconstruction accepted or the last failure. With
 * exact_base_node_count, the base node count must be N: a smaller one is not
 * continued, a rejected reconstruction is not retried, and the result is
 * not_realizable instead.
 */
[[nodiscard]] GeneralizedQuadrature compute_bounded_generalized_quadrature(const MomentSet& moments, std::size_t node_count,
                                                                           const WheelerOptions& options, const GeneralizedOptions& generalized,
                                                                           AbscissaBounds bounds, bool exact_base_node_count);

} // namespace qmomkit::internal
