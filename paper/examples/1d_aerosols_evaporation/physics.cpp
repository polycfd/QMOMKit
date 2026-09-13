#include "physics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace paper_1d_aerosols_evaporation
{
namespace
{

constexpr double residual_diameter = 0.8;
constexpr double maximum_initial_diameter = 1000.0;
constexpr double evaporation_coefficient = 500.0;

} // namespace

qmomkit::PhaseSpaceVelocity evaporation_velocity()
{
    return [](double, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
    {
        if (dimensions != 1)
        {
            throw std::invalid_argument("Aerosol evaporation requires one diameter coordinate.");
        }

        for (std::size_t node = 0; node < coordinates.size(); ++node)
        {
            const double diameter = coordinates[node];
            velocities[node] = diameter <= residual_diameter ? 0.0 : -evaporation_coefficient / (2.0 * diameter);
        }
    };
}

bool is_within_modeled_support(const double diameter)
{
    return std::isfinite(diameter) && diameter >= residual_diameter && diameter <= maximum_initial_diameter;
}

void enforce_residual_diameter(qmomkit::Quadrature1D& quadrature)
{
    for (double& diameter : quadrature.abscissas)
    {
        diameter = std::max(diameter, residual_diameter);
    }
}

void enforce_residual_diameter(qmomkit::monte_carlo::Ensemble& ensemble)
{
    for (qmomkit::monte_carlo::Sample& sample : ensemble)
    {
        if (sample.size() != 1)
        {
            throw std::invalid_argument("Aerosol evaporation requires one diameter coordinate per sample.");
        }
        sample[0] = std::max(sample[0], residual_diameter);
    }
}

} // namespace paper_1d_aerosols_evaporation
