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
    bool adaptive = true;
    double rmin = 1.0e-10;
    double eabs = 1.0e-10;
    double cutoff = 1.0e-30;
};

/**
 * Low-level Wheeler inversion utility.
 *
 * Most applications should use `compute_qmom(...)` with `QmomCase`. Direct
 * Wheeler callers must provide raw moments with m_k at index k and at least
 * 2 * max_nodes values.
 */
[[nodiscard]] Quadrature1D compute_quadrature(const MomentSet& moments, std::size_t max_nodes, const WheelerOptions& options = {});

} // namespace qmomkit
