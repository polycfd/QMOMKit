#include "internal/coordinate_scales.hpp"

#include <cmath>
#include <cstddef>
#include <limits>

namespace qmomkit::internal
{
namespace
{
[[nodiscard]] bool finite_positive(const long double value)
{
    return std::isfinite(value) && value > 0.0L;
}

[[nodiscard]] bool convert_to_double(const long double value, double& converted)
{
    if (!std::isfinite(value) || value < -static_cast<long double>(std::numeric_limits<double>::max()) ||
        value > static_cast<long double>(std::numeric_limits<double>::max()))
    {
        return false;
    }
    converted = static_cast<double>(value);
    return value == 0.0L || converted != 0.0;
}

[[nodiscard]] double coordinate_scale(const double m0, const double m1, const double m2)
{
    if (!std::isfinite(m0) || m0 <= 0.0)
    {
        return 1.0;
    }

    double converted = 1.0;
    if (std::isfinite(m2) && m2 > 0.0 &&
        convert_to_double(std::sqrt(static_cast<long double>(m2)) / std::sqrt(static_cast<long double>(m0)), converted) &&
        finite_positive(static_cast<long double>(converted)))
    {
        return converted;
    }

    if (std::isfinite(m1) && convert_to_double(std::abs(static_cast<long double>(m1)) / static_cast<long double>(m0), converted) &&
        finite_positive(static_cast<long double>(converted)))
    {
        return converted;
    }

    return 1.0;
}
} // namespace

std::vector<double> cqmom_coordinate_scales(const CqmomMomentData& moments)
{
    const std::size_t dimensions = moments.layout.hierarchy.number_of_dimensions();
    std::vector<double> scales(dimensions, 1.0);
    const MultiIndex zero(dimensions, 0);
    const double m0 = moments.contains_required_index(zero) ? moments.moment(zero) : 0.0;

    for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
    {
        MultiIndex first(dimensions, 0);
        MultiIndex second(dimensions, 0);
        first[dimension] = 1;
        second[dimension] = 2;
        const double m1 = moments.contains_required_index(first) ? moments.moment(first) : 0.0;
        const double m2 = moments.contains_required_index(second) ? moments.moment(second) : 0.0;
        scales[dimension] = coordinate_scale(m0, m1, m2);
    }
    return scales;
}

} // namespace qmomkit::internal
