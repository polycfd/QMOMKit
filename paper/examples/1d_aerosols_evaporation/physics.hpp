#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/node_advection.hpp"
#include "qmomkit/quadrature.hpp"

namespace paper_1d_aerosols_evaporation
{

[[nodiscard]] qmomkit::PhaseSpaceVelocity evaporation_velocity();

[[nodiscard]] bool is_within_modeled_support(double diameter);

void enforce_residual_diameter(qmomkit::Quadrature1D& quadrature);

void enforce_residual_diameter(qmomkit::monte_carlo::Ensemble& ensemble);

} // namespace paper_1d_aerosols_evaporation
