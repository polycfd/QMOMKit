#pragma once

#include "checked_size.hpp"
#include "qmomkit/qmom.hpp"

#include <cstddef>
#include <optional>

namespace qmomkit::internal
{

// The moments m_0, ..., m_{M-1} a reconstruction with node count n reads:
// M = 2n for the classical family, 2n + 1 for the hyperbolic one; empty when
// M does not fit in std::size_t. 2n is even, so 2n + 1 fits once 2n does.
[[nodiscard]] inline std::optional<std::size_t> qmom_required_moment_count(const std::size_t node_count, const QmomFamily family)
{
    const std::optional<std::size_t> even_moment_count = checked_size_multiply(2, node_count);
    if (!even_moment_count.has_value())
    {
        return std::nullopt;
    }
    switch (family)
    {
    case QmomFamily::classical:
        return even_moment_count;
    case QmomFamily::hyperbolic:
        return *even_moment_count + 1;
    }
    return even_moment_count;
}

} // namespace qmomkit::internal
