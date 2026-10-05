#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/moment_source.hpp"

#include <random>

namespace cqmom_5d_non_spherical_particles
{

constexpr unsigned long long breakage_seed = 67890;

struct BinaryFragmentationParameters
{
    double fragmentation_rate_constant = 1.0e8;
    double fragmentation_power = 2.0;
    double minimum_fragmentation_diameter = 1.0e-6;
    double shape_relaxation_factor = 0.5;
    double recoil_velocity_standard_deviation = 0.5;
};

[[nodiscard]] BinaryFragmentationParameters binary_fragmentation_parameters();

// Builds the CQMOM moment source for binary fragmentation from the example
// parameters. Terminal nodes are supplied by the source context.
[[nodiscard]] qmomkit::CqmomMomentSource binary_breakage_source(const BinaryFragmentationParameters& parameters);

// Applies stochastic binary-breakage events to the independent particle
// reference over one source interval.
void apply_binary_breakage(qmomkit::monte_carlo::Ensemble& ensemble, double step_size, std::mt19937_64& generator,
                           const BinaryFragmentationParameters& parameters);

} // namespace cqmom_5d_non_spherical_particles
