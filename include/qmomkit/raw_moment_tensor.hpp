#pragma once

#include <cstddef>
#include <vector>

namespace qmomkit
{

/**
 * Global raw moment tensor flattened in phase-space dimension order.
 *
 * `values[flatten_multi_index(alpha, extents)]` stores the raw mixed moment
 * with phase-space-order multi-index `alpha`.
 */
struct RawMomentTensor
{
    std::vector<std::size_t> extents;
    std::vector<double> values{};
};

} // namespace qmomkit
