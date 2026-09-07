#pragma once

#include <cstddef>
#include <limits>
#include <optional>

namespace qmomkit::internal
{

[[nodiscard]] inline std::optional<std::size_t> checked_size_add(const std::size_t lhs, const std::size_t rhs) noexcept
{
    if (lhs > std::numeric_limits<std::size_t>::max() - rhs)
    {
        return std::nullopt;
    }
    return lhs + rhs;
}

[[nodiscard]] inline std::optional<std::size_t> checked_size_multiply(const std::size_t lhs, const std::size_t rhs) noexcept
{
    if (rhs != 0 && lhs > std::numeric_limits<std::size_t>::max() / rhs)
    {
        return std::nullopt;
    }
    return lhs * rhs;
}

} // namespace qmomkit::internal
