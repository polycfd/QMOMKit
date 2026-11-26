#pragma once

#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

// Particles in a flow, in SI units. The phase-space coordinates are, in this
// order, the particle velocity u in m/s and the particle diameter d in m.
namespace cqmom_2d_hyperbolic_drag_dispersion
{

// Density of the particles in kg/m^3.
constexpr double particle_density = 1000.0;
// Dynamic viscosity of the fluid in Pa s.
constexpr double fluid_viscosity = 1.8e-5;
// Velocity of the flow in m/s.
constexpr double flow_velocity = 0.5;
// Standard deviation in m/s of the particle velocity about the flow velocity
// once the turbulent dispersion has balanced the drag.
constexpr double turbulent_velocity_scale = 0.1;

/**
 * The Stokes response time in s of a particle of diameter `diameter` in m,
 * rho_p d^2 / (18 mu_f).
 */
[[nodiscard]] double response_time(double diameter);

/**
 * The velocity diffusion coefficient in m^2/s^3 of a particle of diameter
 * `diameter` in m, sigma_t^2 / tau_p(d), with which the velocity of every size
 * relaxes to the variance sigma_t^2 about the flow velocity.
 */
[[nodiscard]] double dispersion_coefficient(double diameter);

/**
 * Stokes drag: du/dt = (U - u) / tau_p(d), dd/dt = 0.
 */
[[nodiscard]] qmomkit::PhaseSpaceVelocity drag_velocity();

/**
 * Turbulent dispersion as velocity diffusion: dM[k,l]/dt =
 * k (k - 1) times the integral of D(d) u^(k-2) d^l n, evaluated over the
 * nodes of the source-stage quadrature, since D is not a polynomial in d.
 */
[[nodiscard]] qmomkit::CqmomMomentSource dispersion_source();

} // namespace cqmom_2d_hyperbolic_drag_dispersion
