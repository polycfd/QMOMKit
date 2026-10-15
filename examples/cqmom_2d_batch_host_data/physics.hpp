#pragma once

#include "host_data.hpp"
#include "qmomkit/phase_space_evolution.hpp"

namespace cqmom_2d_batch_host_data
{

// Both callbacks capture the host data by reference and read it on every
// invocation, so the host publishes new levels between steps without
// rebuilding the callbacks. The host data must outlive them, so binding a
// temporary is rejected at compile time.

// Phase-space velocity for states ordered as (V, u). The volume relaxes toward
// 1 and the particle velocity toward the local carrier-flow velocity:
//
//     dV/dt = 1 - V,    du/dt = u_f(x_i, t) - u.
[[nodiscard]] qmomkit::BatchPhaseSpaceVelocity advection_velocity(const LocalHostData& host_data);
qmomkit::BatchPhaseSpaceVelocity advection_velocity(LocalHostData&& host_data) = delete;

// Binary fragmentation V -> V/2 at the local fragmentation rate lambda_i,
// leaving u unchanged: dM_jk/dt = lambda_i (2^(1-j) - 1) M_jk.
[[nodiscard]] qmomkit::BatchCqmomMomentSource binary_fragmentation_source(const LocalHostData& host_data);
qmomkit::BatchCqmomMomentSource binary_fragmentation_source(LocalHostData&& host_data) = delete;

// Exact solutions of this model at one location, which the console summary
// prints next to the computed values.
//
// Advection leaves the number density unchanged and fragmentation multiplies it
// by exp(lambda_i (t - t0)).
[[nodiscard]] double exact_number_density(double fragmentation_rate, double initial_time, double initial_number_density, double time);

// Fragmentation leaves u unchanged, so the mean particle velocity
// m = M[0,1] / M[0,0] obeys dm/dt = u_f(x_i, t) - m. For the host's
// traveling-wave carrier flow this has a closed-form solution.
[[nodiscard]] double exact_mean_velocity(double position, double initial_time, double initial_mean_velocity, double time);

} // namespace cqmom_2d_batch_host_data
