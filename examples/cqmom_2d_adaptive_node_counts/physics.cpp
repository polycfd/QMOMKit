#include "physics.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cqmom_2d_adaptive_node_counts
{
namespace
{

constexpr double time_tolerance = 1.0e-12;
constexpr std::size_t velocity_dimension = 1;

double binomial_coefficient(const std::size_t n, const std::size_t k)
{
    if (k > n)
    {
        return 0.0;
    }

    const std::size_t reduced_k = std::min(k, n - k);
    double coefficient = 1.0;
    for (std::size_t i = 1; i <= reduced_k; ++i)
    {
        coefficient *= static_cast<double>(n + 1 - i);
        coefficient /= static_cast<double>(i);
    }
    return coefficient;
}

double raw_moment(const qmomkit::MultiIndex& order, const std::vector<double>& weights, const std::vector<double>& sizes,
                  const std::vector<double>& velocities)
{
    double moment = 0.0;
    for (std::size_t node = 0; node < weights.size(); ++node)
    {
        moment += weights[node] * std::pow(sizes[node], static_cast<double>(order[0])) * std::pow(velocities[node], static_cast<double>(order[1]));
    }
    return moment;
}

double direct_symmetric_source_rate(const qmomkit::MultiIndex& order, const VarianceInjectionParameters& parameters,
                                    const std::vector<double>& weights, const std::vector<double>& sizes, const std::vector<double>& velocities)
{
    double rate = 0.0;
    for (std::size_t node = 0; node < weights.size(); ++node)
    {
        const double u = velocities[node];
        const double velocity_order = static_cast<double>(order[1]);
        const double redistributed =
            0.5 * (std::pow(u - parameters.velocity_displacement, velocity_order) + std::pow(u + parameters.velocity_displacement, velocity_order)) -
            std::pow(u, velocity_order);
        rate += weights[node] * std::pow(sizes[node], static_cast<double>(order[0])) * redistributed;
    }
    return parameters.source_rate * rate;
}

} // namespace

VarianceInjectionParameters default_variance_injection_parameters()
{
    return {};
}

std::size_t velocity_phase_space_dimension()
{
    return velocity_dimension;
}

bool source_pulse_active(const double time, const VarianceInjectionParameters& parameters)
{
    return time + time_tolerance >= parameters.pulse_start && time < parameters.pulse_end - time_tolerance;
}

qmomkit::PhaseSpaceVelocity advection_velocity(const VarianceInjectionParameters& parameters)
{
    return qmomkit::make_phase_space_velocity(
        [parameters](double, const std::vector<double>& coordinates)
        {
            if (coordinates.size() < 2)
            {
                throw std::invalid_argument("The adaptive-node-count example expects phase-space coordinates ordered as (d,u).");
            }

            const double u = coordinates[velocity_dimension];
            return std::vector<double>{0.0, -parameters.relaxation_rate * (u - parameters.equilibrium_velocity)};
        });
}

qmomkit::CqmomMomentSource variance_injection_source(const VarianceInjectionParameters& parameters)
{
    return [parameters](const double time, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> derivative(moments.size(), 0.0);
        if (!source_pulse_active(time, parameters))
        {
            return qmomkit::CqmomMomentRate{std::move(derivative)};
        }

        for (std::size_t entry = 0; entry < moments.size(); ++entry)
        {
            const qmomkit::MultiIndex& order = moments.index(entry);
            if (order.size() != 2)
            {
                throw std::invalid_argument("The variance-injection source expects 2D raw moments ordered as (d,u).");
            }

            const std::size_t velocity_order = order[1];
            double redistribution_change = 0.0;
            for (std::size_t displacement_order = 2; displacement_order <= velocity_order; displacement_order += 2)
            {
                qmomkit::MultiIndex lower_order = order;
                lower_order[1] = velocity_order - displacement_order;
                const std::optional<double> lower_moment = moments.value(lower_order);
                if (!lower_moment.has_value())
                {
                    throw std::invalid_argument("The variance-injection source requires lower-order velocity moments from the CQMOM context.");
                }

                redistribution_change += binomial_coefficient(velocity_order, displacement_order) *
                                         std::pow(parameters.velocity_displacement, static_cast<double>(displacement_order)) * *lower_moment;
            }

            derivative[entry] = parameters.source_rate * redistribution_change;
        }

        return qmomkit::CqmomMomentRate{std::move(derivative)};
    };
}

bool validate_variance_injection_source_rates(const VarianceInjectionParameters& parameters, std::string& message)
{
    const std::vector<double> weights{0.2, 0.35, 0.45};
    const std::vector<double> sizes{1.25, 2.0, 2.75};
    const std::vector<double> velocities{-0.4, 0.15, 0.9};
    std::vector<qmomkit::MultiIndex> indices;
    indices.reserve(12);
    for (const std::size_t size_order : {std::size_t{0}, std::size_t{2}})
    {
        for (std::size_t velocity_order = 0; velocity_order <= 5; ++velocity_order)
        {
            indices.push_back(qmomkit::MultiIndex{size_order, velocity_order});
        }
    }

    std::vector<double> values;
    values.reserve(indices.size());
    for (const qmomkit::MultiIndex& index : indices)
    {
        values.push_back(raw_moment(index, weights, sizes, velocities));
    }

    const qmomkit::PhaseSpaceHierarchy hierarchy{{2, 3}};
    const qmomkit::CqmomSourceContext active_context{parameters.pulse_start, hierarchy, indices, values};
    const qmomkit::CqmomMomentRate active_rate = variance_injection_source(parameters)(parameters.pulse_start, active_context);
    if (active_rate.values.size() != indices.size())
    {
        message = "The source-rate self-check received a derivative vector with the wrong size.";
        return false;
    }

    double maximum_error = 0.0;
    for (std::size_t entry = 0; entry < indices.size(); ++entry)
    {
        const double expected = direct_symmetric_source_rate(indices[entry], parameters, weights, sizes, velocities);
        maximum_error = std::max(maximum_error, std::abs(active_rate.values[entry] - expected));
    }
    if (maximum_error > 2.0e-13)
    {
        message = "The source-rate self-check failed against direct quadrature evaluation.";
        return false;
    }

    const qmomkit::CqmomSourceContext inactive_context{parameters.pulse_start - 1.0, hierarchy, indices, values};
    const qmomkit::CqmomMomentRate inactive_rate = variance_injection_source(parameters)(parameters.pulse_start - 1.0, inactive_context);
    for (const double value : inactive_rate.values)
    {
        if (std::abs(value) > std::numeric_limits<double>::epsilon())
        {
            message = "The source-rate self-check found a nonzero rate outside the pulse.";
            return false;
        }
    }

    std::ostringstream stream;
    stream << "max error: " << std::scientific << maximum_error;
    message = stream.str();
    return true;
}

} // namespace cqmom_2d_adaptive_node_counts
