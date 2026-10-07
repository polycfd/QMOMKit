#include "internal/reconstruction_scaling.hpp"

#include "internal/recursive_quadrature.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

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

[[nodiscard]] double automatic_scale(const double m0, const double m1, const double m2)
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

[[nodiscard]] bool scale_value(const double value, const double scale, const std::size_t exponent, double& scaled)
{
    long double transformed = static_cast<long double>(value);
    for (std::size_t power = 0; power < exponent; ++power)
    {
        transformed /= static_cast<long double>(scale);
        if (value != 0.0 && transformed == 0.0L)
        {
            return false;
        }
    }
    return convert_to_double(transformed, scaled);
}

[[nodiscard]] bool scale_mixed_value(const double value, const MultiIndexView index, const std::vector<double>& scales, double& scaled)
{
    if (index.size() != scales.size())
    {
        return false;
    }
    long double transformed = static_cast<long double>(value);
    for (std::size_t dimension = 0; dimension < index.size(); ++dimension)
    {
        for (std::size_t power = 0; power < index[dimension]; ++power)
        {
            transformed /= static_cast<long double>(scales[dimension]);
            if (value != 0.0 && transformed == 0.0L)
            {
                return false;
            }
        }
    }
    return convert_to_double(transformed, scaled);
}

[[nodiscard]] bool scale_bounds(AbscissaBounds& bounds, const double scale)
{
    double lower = 0.0;
    double upper = 0.0;
    if (!scale_value(bounds.lower, scale, 1, lower) || !scale_value(bounds.upper, scale, 1, upper))
    {
        return false;
    }
    bounds = AbscissaBounds{lower, upper};
    return true;
}
} // namespace

std::optional<QmomScalingContext> make_qmom_scaling_context(const std::vector<double>& moments, std::optional<AbscissaBounds> abscissa_bounds)
{
    const double scale =
        automatic_scale(moments.empty() ? 0.0 : moments[0], moments.size() > 1 ? moments[1] : 0.0, moments.size() > 2 ? moments[2] : 0.0);
    QmomScalingContext context{moments, std::move(abscissa_bounds), scale};
    for (std::size_t order = 0; order < context.moments.size(); ++order)
    {
        if (!scale_value(context.moments[order], scale, order, context.moments[order]))
        {
            return std::nullopt;
        }
    }
    if (context.abscissa_bounds.has_value() && !scale_bounds(*context.abscissa_bounds, scale))
    {
        return std::nullopt;
    }
    return context;
}

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
        scales[dimension] = automatic_scale(m0, m1, m2);
    }
    return scales;
}

std::optional<CqmomScalingContext> make_cqmom_scaling_context(const CqmomMomentData& moments, const CqmomAbscissaBounds& coordinate_support)
{
    if (moments.values.size() != moments.layout.required_indices.size())
    {
        return std::nullopt;
    }
    const std::size_t dimensions = moments.layout.hierarchy.number_of_dimensions();
    CqmomScalingContext context{moments, coordinate_support, cqmom_coordinate_scales(moments)};

    for (std::size_t entry = 0; entry < context.moments.values.size(); ++entry)
    {
        double scaled = 0.0;
        const MultiIndexView index = context.moments.layout.required_indices[entry];
        if (!scale_mixed_value(context.moments.values[entry], index, context.scales, scaled))
        {
            return std::nullopt;
        }
        context.moments.values[entry] = scaled;
    }

    if (!context.coordinate_support.empty())
    {
        if (context.coordinate_support.size() != dimensions)
        {
            return std::nullopt;
        }
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
        {
            if (context.coordinate_support[dimension].has_value() && !scale_bounds(*context.coordinate_support[dimension], context.scales[dimension]))
            {
                return std::nullopt;
            }
        }
    }
    return context;
}

bool restore_abscissas(std::vector<double>& abscissas, const double scale)
{
    for (double& abscissa : abscissas)
    {
        const long double transformed = static_cast<long double>(abscissa) * static_cast<long double>(scale);
        double physical = 0.0;
        if ((abscissa != 0.0 && scale != 0.0 && transformed == 0.0L) || !convert_to_double(transformed, physical))
        {
            return false;
        }
        abscissa = physical;
    }
    return true;
}

bool restore_cqmom_abscissas(std::vector<std::vector<double>>& abscissas, const std::vector<double>& scales)
{
    if (abscissas.size() != scales.size())
    {
        return false;
    }
    for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
    {
        if (!restore_abscissas(abscissas[dimension], scales[dimension]))
        {
            return false;
        }
    }
    return true;
}

bool restore_bounds(AbscissaBounds& bounds, const double scale)
{
    const long double lower_value = static_cast<long double>(bounds.lower) * static_cast<long double>(scale);
    const long double upper_value = static_cast<long double>(bounds.upper) * static_cast<long double>(scale);
    double lower = 0.0;
    double upper = 0.0;
    if ((bounds.lower != 0.0 && scale != 0.0 && lower_value == 0.0L) || (bounds.upper != 0.0 && scale != 0.0 && upper_value == 0.0L) ||
        !convert_to_double(lower_value, lower) || !convert_to_double(upper_value, upper))
    {
        return false;
    }
    bounds = AbscissaBounds{lower, upper};
    return true;
}

bool can_restore_cqmom_result(const RecursiveQuadratureResult& result, const std::vector<double>& scales)
{
    std::vector<std::vector<double>> abscissas = result.quadrature.abscissas();
    if (!restore_cqmom_abscissas(abscissas, scales))
    {
        return false;
    }
    for (const RecursiveQuadratureBranchDiagnostics& diagnostic : result.diagnostics)
    {
        if (diagnostic.phase_space_dimension >= scales.size())
        {
            return false;
        }
        const double scale = scales[diagnostic.phase_space_dimension];
        if (diagnostic.bounded_support_interval_present)
        {
            AbscissaBounds bounds = diagnostic.bounded_support_interval;
            if (!restore_bounds(bounds, scale))
            {
                return false;
            }
        }
        std::vector<double> nodes = diagnostic.bounded_returned_nodes;
        if (!restore_abscissas(nodes, scale))
        {
            return false;
        }
    }
    return true;
}

} // namespace qmomkit::internal
