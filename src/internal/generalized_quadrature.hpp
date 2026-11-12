#pragma once

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
 * N = node_count, on the real line, with the N_G = generalized.node_count nodes
 * of the continued recurrence. Adaptive mode takes the classical selection of
 * compute_quadrature(...) for N nodes, with rmin and eabs, and fixed mode the
 * classical N-node checks. The k-node quadrature found is continued to N_G nodes
 * when sigma(k, k) exceeds its rounding bound tau_k, and returned as it is
 * when sigma(k, k) is zero within tau_k. Below -tau_k, fixed mode is
 * not_realizable and adaptive mode returns the classical quadrature. A
 * non-finite sigma(k, k) or tau_k is not_realizable in fixed mode. In adaptive
 * mode, an overflowed +-inf returns the classical quadrature and a NaN is
 * not_realizable. The coordinate scale is zero, as in compute_quadrature(...).
 */
[[nodiscard]] GeneralizedQuadrature compute_generalized_quadrature(const MomentSet& moments, std::size_t node_count, const WheelerOptions& options,
                                                                   const GeneralizedOptions& generalized);

} // namespace qmomkit::internal
