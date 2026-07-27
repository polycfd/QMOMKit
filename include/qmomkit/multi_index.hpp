#pragma once

#include <cstddef>
#include <vector>

namespace qmomkit
{

using MultiIndex = std::vector<std::size_t>;

/**
 * Return the row-major product of extents.
 *
 * The last index is the fastest-varying index. All extents must be positive.
 */
[[nodiscard]] std::size_t product_extents(const std::vector<std::size_t>& extents);

/**
 * Convert a row-major multi-index to a flat index.
 *
 * The supplied multi-index and extents use the same axis ordering. The last
 * axis is fastest varying, matching C/C++ row-major storage.
 */
[[nodiscard]] std::size_t flatten_multi_index(const MultiIndex& multi_index, const std::vector<std::size_t>& extents);

/**
 * Convert a row-major flat index to a multi-index.
 */
[[nodiscard]] MultiIndex unflatten_multi_index(std::size_t flat_index, const std::vector<std::size_t>& extents);

/**
 * Enumerate all row-major multi-indices for the supplied extents.
 */
[[nodiscard]] std::vector<MultiIndex> enumerate_multi_indices(const std::vector<std::size_t>& extents);

} // namespace qmomkit
