#pragma once

#include "qmomkit/moment_set.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace qmomkit
{

enum class QmomStatus
{
    success,
    invalid_input,
    not_realizable,
    no_adaptive_quadrature
};

/**
 * Numerical controls for one-dimensional QMOM quadrature construction.
 */
struct QmomOptions
{
    bool adaptive = false;
    double rmin = 1.0e-10;
    double eabs = 1.0e-10;
    double cutoff = 1.0e-30;
};

/**
 * User-facing case for one-dimensional QMOM quadrature construction.
 *
 * `moments` stores unnormalized raw moments in ascending order:
 * `moments[k] == m_k`. A case with `node_count == n` requires at least
 * `2 * n` moments.
 */
struct QmomCase
{
    MomentSet moments;
    std::size_t node_count = 0;
    QmomOptions options;
};

struct QmomResult
{
    std::vector<double> abscissas;
    std::vector<double> weights;
    QmomStatus status = QmomStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == QmomStatus::success;
    }
};

/**
 * Construct a one-dimensional QMOM quadrature from raw moments.
 *
 * The input moments are raw moments in ascending order. Returned abscissas are
 * sorted in ascending order, and weights use the same ordering.
 */
[[nodiscard]] QmomResult compute_qmom(const MomentSet& moments, std::size_t node_count, const QmomOptions& options = {});

/**
 * Construct a one-dimensional QMOM quadrature from a case object.
 */
[[nodiscard]] QmomResult compute_qmom(const QmomCase& case_);

} // namespace qmomkit
