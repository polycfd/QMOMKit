#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace qmomkit::internal
{

[[nodiscard]] inline bool all_finite(const std::vector<double>& values)
{
    return std::all_of(values.begin(), values.end(), [](const double value) { return std::isfinite(value); });
}

} // namespace qmomkit::internal
