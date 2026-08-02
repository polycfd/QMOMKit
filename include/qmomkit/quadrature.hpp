#pragma once

#include <string>
#include <vector>

namespace qmomkit
{

enum class QuadratureStatus
{
    success,
    invalid_input,
    not_realizable
};

/**
 * One-dimensional quadrature rule produced by Wheeler inversion.
 *
 * Abscissas are sorted in ascending order and weights use the same ordering.
 * On failure, the result contains fallback one-node values only where that
 * behavior is mathematically defined.
 */
struct Quadrature1D
{
    std::vector<double> abscissas;
    std::vector<double> weights;
    QuadratureStatus status = QuadratureStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == QuadratureStatus::success;
    }
};

} // namespace qmomkit
