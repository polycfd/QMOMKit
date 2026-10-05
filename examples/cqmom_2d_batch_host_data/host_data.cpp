#include "host_data.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

namespace cqmom_2d_batch_host_data
{
namespace
{
// The fragmentation rate grows linearly across the locations.
constexpr double minimum_fragmentation_rate = 0.25;
constexpr double maximum_fragmentation_rate = 0.75;

[[nodiscard]] double carrier_flow_velocity(const double position, const double time)
{
    return carrier_amplitude * std::sin(2.0 * pi * (position - carrier_wave_speed * time));
}
} // namespace

LocalHostData make_local_host_data(const std::size_t location_count, const double initial_time)
{
    LocalHostData host_data;
    host_data.position.resize(location_count);
    host_data.fragmentation_rate.resize(location_count);
    host_data.carrier_velocity_current.resize(location_count);
    host_data.carrier_velocity_previous.resize(location_count);

    for (std::size_t location = 0; location < location_count; ++location)
    {
        const double fraction = static_cast<double>(location) / static_cast<double>(location_count);
        host_data.position[location] = fraction;
        host_data.fragmentation_rate[location] = minimum_fragmentation_rate + fraction * (maximum_fragmentation_rate - minimum_fragmentation_rate);
        host_data.carrier_velocity_current[location] = carrier_flow_velocity(host_data.position[location], initial_time);
    }
    host_data.carrier_velocity_previous = host_data.carrier_velocity_current;
    host_data.sample_time = initial_time;
    host_data.sample_interval = 0.0;
    return host_data;
}

void publish_carrier_flow_level(LocalHostData& host_data, const double time)
{
    host_data.carrier_velocity_previous = host_data.carrier_velocity_current;
    for (std::size_t location = 0; location < host_data.position.size(); ++location)
    {
        host_data.carrier_velocity_current[location] = carrier_flow_velocity(host_data.position[location], time);
    }
    host_data.sample_interval = time - host_data.sample_time;
    host_data.sample_time = time;
}

double carrier_velocity_at(const LocalHostData& host_data, const std::size_t location, const double time)
{
    const double current = host_data.carrier_velocity_current[location];
    if (!(host_data.sample_interval > 0.0))
    {
        return current;
    }
    const double s = (time - host_data.sample_time) / host_data.sample_interval;
    return current + s * (current - host_data.carrier_velocity_previous[location]);
}

} // namespace cqmom_2d_batch_host_data
