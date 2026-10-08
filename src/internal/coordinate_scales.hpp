#pragma once

#include "qmomkit/cqmom.hpp"

#include <vector>

namespace qmomkit::internal
{

// The scale of each phase-space coordinate, from the marginal moments of the
// whole distribution: sqrt(m2 / m0), else |m1| / m0, else one.
[[nodiscard]] std::vector<double> cqmom_coordinate_scales(const CqmomMomentData& moments);

} // namespace qmomkit::internal
