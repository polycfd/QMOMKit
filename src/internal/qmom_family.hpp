#pragma once

#include "checked_size.hpp"
#include "qmomkit/qmom.hpp"

#include <cstddef>
#include <optional>

namespace qmomkit::internal
{

// The moments m_0, ..., m_{M-1} a reconstruction with node count N reads:
// M = 2N for the classical family, 2N + 1 for the hyperbolic and generalized
// ones; empty when M does not fit in std::size_t. 2N is even, so 2N + 1 fits
// once 2N does.
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
    case QmomFamily::generalized:
        return *even_moment_count + 1;
    }
    return even_moment_count;
}

// compute_qmom(...) for a growth trial of the evolution, which accepts only
// node_count as the base node count. The generalized family then neither
// continues a smaller base node count nor retries below a rejected one, and
// returns not_realizable instead. For the classical and hyperbolic families,
// it is compute_qmom(...).
[[nodiscard]] QmomResult compute_qmom_growth_trial(const MomentSet& moments, std::size_t node_count, const QmomOptions& options,
                                                   std::optional<AbscissaBounds> abscissa_bounds);

} // namespace qmomkit::internal
