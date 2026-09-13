#include "physics.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace qmom_1d_aerosol_evaporation
{
namespace
{

constexpr double residual_diameter = 0.8;         // micrometres
constexpr double evaporation_coefficient = 500.0; // micrometres squared per second

} // namespace

qmomkit::PhaseSpaceVelocity evaporation_velocity()
{
    // This example deliberately keeps the per-node convenience form to show the
    // alternative; every other example writes a `PhaseSpaceVelocity` directly.
    // `make_phase_space_velocity` calls the per-node callback once per node and
    // allocates one returned vector per call, so it costs more per
    // right-hand-side evaluation than a directly written packed callback. Prefer
    // the packed form when advancing many nodes; reach for this one when an
    // existing per-node model is not worth rewriting.
    return qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>& coordinates)
        {
            if (coordinates.size() != 1)
            {
                throw std::invalid_argument("Aerosol evaporation requires one diameter coordinate.");
            }

            // The phase-space coordinate is droplet diameter d in micrometres. The
            // d-squared law, d(d^2)/dt = -evaporation_coefficient, gives dd/dt =
            // -evaporation_coefficient/(2d), so evaporation moves nodes to smaller diameters.
            const double diameter = coordinates[0];
            if (diameter <= residual_diameter)
            {
                return std::vector<double>{0.0};
            }

            return std::vector<double>{-evaporation_coefficient / (2.0 * diameter)};
        });
}

void enforce_residual_diameter(qmomkit::Quadrature1D& quadrature)
{
    // Retain the 0.8-micrometre residual diameter after the native node-advection step.
    for (double& diameter : quadrature.abscissas)
    {
        diameter = std::max(diameter, residual_diameter);
    }
}

} // namespace qmom_1d_aerosol_evaporation
