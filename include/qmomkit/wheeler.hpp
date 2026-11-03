#pragma once

#include "qmomkit/moment_set.hpp"
#include "qmomkit/quadrature.hpp"

#include <cstddef>
#include <string>
#include <vector>

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

/**
 * Low-level diagnostic of the Wheeler inversion: the recurrence of
 * m_0, ..., m_{2n-1}, as every QMOM reconstruction computes it. Index
 * k = 0, ..., n - 1: the recurrence coefficients a_k and b_k, the diagonal
 * sigma(k, k) of the Wheeler table and tau_k, the rounding bound of
 * sigma(k, k). Entry 0: a = m1/m0, b = 0, sigma = m0, tolerance = 0 (no test at
 * k = 0). Every value is the recurrence as computed, so any of them may be
 * non-finite, for instance for a tiny m0. At the first k with
 * sigma(k, k) <= tau_k, sigma(k, k), tau_k and b_k still hold the test; a_k,
 * which divides by sigma(k, k), and all later entries carry no information. The
 * vectors are empty unless the status is success.
 */
struct WheelerRecurrence
{
    std::vector<double> a;         // a_k
    std::vector<double> b;         // b_k
    std::vector<double> sigma;     // sigma(k, k)
    std::vector<double> tolerance; // tau_k, the rounding bound of sigma(k, k)
    QuadratureStatus status = QuadratureStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == QuadratureStatus::success;
    }
};

/**
 * The Wheeler recurrence of the first 2 * node_count moments, with the values
 * on which compute_quadrature(...) tests how many nodes the moments support. No
 * acceptance test and no cutoff on m0 are applied. The status is invalid_input
 * for a zero node_count, fewer than 2 * node_count moments or a non-finite one
 * among them, and not_realizable for m0 < 0 and m0 = 0.
 */
[[nodiscard]] WheelerRecurrence compute_wheeler_recurrence(const MomentSet& moments, std::size_t node_count);

} // namespace qmomkit
