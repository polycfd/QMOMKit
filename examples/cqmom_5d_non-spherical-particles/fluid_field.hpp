#pragma once

#include <cmath>

namespace cqmom_5d_non_spherical_particles
{

struct BackgroundFluidVelocity
{
    double u = 0.0;
    double v = 0.0;
    double w = 0.0;
};

// What the driver owns and the particle callback reads.
//
// A coupled host cannot evaluate its fields at a library's internal stage times,
// and it does not know them at t + dt either: all it has is the level it just
// computed and the level it computed before. Keeping both is enough. The
// callback reconstructs the field linearly in time from them, which is
// second-order accurate, so the coupling stops being the limiting error of the
// second-order Strang-split step. Handing over a single frozen value instead
// would make the coupling first-order.
struct CarrierFlowSamples
{
    BackgroundFluidVelocity current;  // host field at sample_time
    BackgroundFluidVelocity previous; // host field one outer step earlier
    double sample_time = 0.0;
    // Time between the two stored levels. Zero means no previous level exists
    // yet, which is the case on the first step only; the reconstruction then
    // holds current constant, costing first-order accuracy for that one step.
    double sample_interval = 0.0;
};

// Reconstructs the carrier flow at the stage time QMOMKit passes the callback by
// extending the straight line through the two stored host levels.
[[nodiscard]] inline BackgroundFluidVelocity carrier_flow_at(const CarrierFlowSamples& samples, const double time)
{
    if (!(samples.sample_interval > 0.0))
    {
        return samples.current;
    }

    const double s = (time - samples.sample_time) / samples.sample_interval;

    return {samples.current.u + s * (samples.current.u - samples.previous.u), samples.current.v + s * (samples.current.v - samples.previous.v),
            samples.current.w + s * (samples.current.w - samples.previous.w)};
}

// Stands in for the host's flow solver. Only the driver calls this, once per
// outer step, exactly as a host would publish one new field level per step.
[[nodiscard]] inline BackgroundFluidVelocity carrier_flow_velocity(const double time)
{
    constexpr double amplitude_x = 2.0;
    constexpr double amplitude_y = 2.0;
    constexpr double amplitude_z = 0.0;
    constexpr double angular_frequency = 1.0;

    return {amplitude_x * std::cos(angular_frequency * time), amplitude_y * std::sin(angular_frequency * time),
            amplitude_z * std::cos(angular_frequency * time)};
}

} // namespace cqmom_5d_non_spherical_particles
