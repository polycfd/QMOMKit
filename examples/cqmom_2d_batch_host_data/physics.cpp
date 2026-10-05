#include "physics.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cqmom_2d_batch_host_data
{

qmomkit::BatchPhaseSpaceVelocity advection_velocity(const LocalHostData& host_data)
{
    // Capturing the reference parameter by reference refers to the host's data,
    // not to the parameter, so the capture stays valid after this function
    // returns. The same holds for the moment source below.
    return [&host_data](const std::size_t location, const double time, const std::vector<double>& coordinates, std::vector<double>& velocities,
                        const std::size_t dimensions)
    {
        if (dimensions != 2)
        {
            throw std::invalid_argument("The 2D CQMOM advection velocity expects coordinates ordered as (d, u).");
        }

        // The carrier flow depends on the location and the time only, so it is
        // reconstructed once for all nodes of this location.
        const double carrier_velocity = carrier_velocity_at(host_data, location, time);
        for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
        {
            const std::size_t first = node * dimensions;
            const double d = coordinates[first];
            const double u = coordinates[first + 1];
            velocities[first] = 1.0 - d;
            velocities[first + 1] = carrier_velocity - u;
        }
    };
}

qmomkit::BatchCqmomMomentSource binary_fragmentation_source(const LocalHostData& host_data)
{
    return [&host_data](const std::size_t location, double, const qmomkit::CqmomSourceContext& moments)
    {
        const double fragmentation_rate = host_data.fragmentation_rate[location];
        std::vector<double> source(moments.size(), 0.0);
        for (std::size_t index = 0; index < moments.size(); ++index)
        {
            const qmomkit::MultiIndexView order = moments.index(index);
            if (order.size() != 2)
            {
                throw std::invalid_argument("The 2D CQMOM fragmentation source received a moment index with the wrong rank.");
            }

            const double diameter_order = static_cast<double>(order[0]);
            const double daughter_factor = std::pow(2.0, 1.0 - diameter_order);
            source[index] = fragmentation_rate * (daughter_factor - 1.0) * moments.value(index);
        }
        return qmomkit::CqmomMomentRate{std::move(source)};
    };
}

double exact_number_density(const double fragmentation_rate, const double initial_time, const double initial_number_density, const double time)
{
    return initial_number_density * std::exp(fragmentation_rate * (time - initial_time));
}

double exact_mean_velocity(const double position, const double initial_time, const double initial_mean_velocity, const double time)
{
    // With u_f = A sin(phi - omega t), phi = 2 pi x and omega = 2 pi c, the
    // periodic part of the solution is
    //     m_p(t) = A [sin(phi - omega t) + omega cos(phi - omega t)] / (1 + omega^2),
    // and the start-up transient decays as exp(-(t - t0)).
    const double phase = 2.0 * pi * position;
    const double omega = 2.0 * pi * carrier_wave_speed;
    const auto periodic_part = [phase, omega](const double t)
    { return carrier_amplitude * (std::sin(phase - omega * t) + omega * std::cos(phase - omega * t)) / (1.0 + omega * omega); };
    return periodic_part(time) + (initial_mean_velocity - periodic_part(initial_time)) * std::exp(-(time - initial_time));
}

} // namespace cqmom_2d_batch_host_data
