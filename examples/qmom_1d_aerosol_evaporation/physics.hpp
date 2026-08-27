#pragma once

#include "qmomkit/node_advection.hpp"
#include "qmomkit/quadrature.hpp"

namespace qmom_1d_aerosol_evaporation
{

[[nodiscard]] qmomkit::PhaseSpaceVelocity evaporation_velocity();

void enforce_residual_diameter(qmomkit::Quadrature1D& quadrature);

} // namespace qmom_1d_aerosol_evaporation
