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
    not_realizable
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

class QmomState;

namespace internal
{
class QmomStateAccess;
}

/**
 * Accepted one-dimensional QMOM phase-space evolution state.
 *
 * The stored moments are reconstructed from the accepted quadrature at
 * construction. Users can inspect the state but cannot mutate moments and
 * quadrature independently.
 */
class QmomState
{
public:
    QmomState(const QmomState&) = default;
    QmomState& operator=(const QmomState&) = default;
    QmomState(QmomState&&) noexcept = default;
    QmomState& operator=(QmomState&&) noexcept = default;
    ~QmomState() = default;

    [[nodiscard]] const std::vector<double>& moments() const noexcept
    {
        return moments_;
    }

    [[nodiscard]] const QmomResult& quadrature() const noexcept
    {
        return quadrature_;
    }

    [[nodiscard]] std::size_t active_node_count() const noexcept
    {
        return quadrature_.weights.size();
    }

private:
    QmomState(std::vector<double> moments, QmomResult quadrature);
    friend class internal::QmomStateAccess;

    std::vector<double> moments_;
    QmomResult quadrature_;
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

    /**
     * Construct the accepted initial phase-space evolution state from the configured maximum
     * node count and canonicalize stored moments from the accepted quadrature.
     */
    [[nodiscard]] QmomState initial_state() const;
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
