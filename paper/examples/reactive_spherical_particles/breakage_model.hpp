#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/moment_source.hpp"

#include <random>

namespace cqmom_reactive_spherical_particles
{

constexpr unsigned long long breakage_seed = 67890;

struct BinaryFragmentationParameters
{
    double fragmentation_rate_constant = 1.0e8;
    double fragmentation_power = 2.0;
    double minimum_fragmentation_diameter = 1.0e-6;
    double recoil_velocity_standard_deviation = 0.5;
};

[[nodiscard]] BinaryFragmentationParameters binary_fragmentation_parameters();

void validate_binary_fragmentation_parameters(const BinaryFragmentationParameters& parameters);

// Builds the CQMOM moment source for binary fragmentation, capturing the
// example parameters and CQMOM numerical options by value.
[[nodiscard]] qmomkit::CqmomMomentSource binary_breakage_source(const qmomkit::CqmomCase& cqmom_case,
                                                                const BinaryFragmentationParameters& parameters);

// Applies one stochastic moment-source interval in place; the generator advances
// by reference so the example seed controls reproducibility.
void apply_binary_breakage(qmomkit::monte_carlo::Ensemble& ensemble, double step_size, std::mt19937_64& generator,
                           const BinaryFragmentationParameters& parameters);

} // namespace cqmom_reactive_spherical_particles
