#pragma once

#include <cstddef>
#include <vector>

namespace cqmom_2d_batch_host_data
{

// The host's carrier flow is a traveling wave, u_f(x, t) = A sin(2 pi (x - c t)),
// so neighboring locations see the same flow at different phases.
inline constexpr double pi = 3.14159265358979323846;
inline constexpr double carrier_amplitude = 1.0;
inline constexpr double carrier_wave_speed = 0.5;

// The data a host makes available to its QMOMKit callbacks for the locations
// it owns, for example the mesh nodes held by one processor. QMOMKit itself
// never sees this struct; only the callbacks read it.
//
// Every array is indexed by the local location index i. The host stores the
// CqmomState of location i at position i of the state vector it passes to the
// batch calls, and QMOMKit passes that same position to the callbacks as
// batch_index. A callback therefore reads entry batch_index of these arrays
// directly, and one index runs through the host data and the QMOMKit states.
struct LocalHostData
{
    // Fixed per-location data. The position is used only by the host's own
    // flow solver; the fragmentation rate is read by the moment source.
    std::vector<double> position;
    std::vector<double> fragmentation_rate;

    // Time-varying carrier-flow velocity at every location: the level the host
    // computed last and the level before it. The velocity callback
    // reconstructs the flow linearly in time from these two levels at the
    // stage time QMOMKit passes it, because the host cannot evaluate its field
    // at QMOMKit's internal stage times.
    std::vector<double> carrier_velocity_current;
    std::vector<double> carrier_velocity_previous;
    double sample_time = 0.0;
    // Time between the two stored levels. Zero means no previous level exists
    // yet, which is the case on the first step only; the reconstruction then
    // holds the current level constant for that step.
    double sample_interval = 0.0;
};

// Builds the host data for location_count locations spread over [0, 1) and
// publishes the carrier-flow level at initial_time.
[[nodiscard]] LocalHostData make_local_host_data(std::size_t location_count, double initial_time);

// Stands in for the host's flow solver: computes the carrier flow at every
// location at the given time and publishes it as the newest level, keeping the
// level it replaces.
void publish_carrier_flow_level(LocalHostData& host_data, double time);

// Carrier-flow velocity at one location, reconstructed linearly in time from
// the two stored levels.
[[nodiscard]] double carrier_velocity_at(const LocalHostData& host_data, std::size_t location, double time);

} // namespace cqmom_2d_batch_host_data
