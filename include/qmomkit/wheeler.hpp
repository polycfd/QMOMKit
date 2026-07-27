#pragma once

#include "qmomkit/moment_set.hpp"
#include "qmomkit/quadrature.hpp"

#include <cstddef>

namespace qmomkit
{

/**
 * Numerical controls for the one-dimensional Wheeler algorithm.
 *
 * Defaults follow the documented QMOMKit Wheeler convention.
 */
struct WheelerOptions
{
    bool adaptive = false;
    double rmin = 1.0e-10;
    double eabs = 1.0e-10;
    double cutoff = 1.0e-30;
};

/**
 * Invert a one-dimensional raw moment sequence into QMOM weights and abscissas.
 *
 * The input moments must use raw moment ordering m_k at index k and must
 * contain at least 2 * max_nodes values.
 */
[[nodiscard]] Quadrature1D compute_quadrature(const MomentSet& moments, std::size_t max_nodes, const WheelerOptions& options = {});

} // namespace qmomkit
