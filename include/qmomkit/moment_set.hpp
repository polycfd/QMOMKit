#pragma once

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <vector>

namespace qmomkit
{

/**
 * Stores a one-dimensional raw moment sequence.
 *
 * moments[k] is the raw moment of order k, m_k = sum_i w_i x_i^k. A
 * quadrature with n nodes requires moments 0 through 2*n - 1.
 */
class MomentSet
{
public:
    MomentSet() = default;
    explicit MomentSet(std::vector<double> moments);
    MomentSet(std::initializer_list<double> moments);

    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] const std::vector<double>& values() const noexcept;

    [[nodiscard]] double operator[](std::size_t moment_order) const;

private:
    std::vector<double> moments_;
};

} // namespace qmomkit
