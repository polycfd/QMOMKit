#pragma once

#include "monte_carlo/ensemble.hpp"

#include <cstddef>
#include <random>
#include <vector>

namespace qmom_1d_breakage_advection_mc_local_rates
{

constexpr unsigned long long breakage_seed = 67890;
constexpr unsigned long long nucleation_seed = 24680;

struct Particle
{
    std::size_t id = 0;
    double diameter = 0.0;
};

[[nodiscard]] std::vector<Particle> make_particles(const qmomkit::monte_carlo::Ensemble& ensemble);

[[nodiscard]] qmomkit::monte_carlo::Ensemble to_ensemble(const std::vector<Particle>& particles);

void update_particle_diameters(std::vector<Particle>& particles, const qmomkit::monte_carlo::Ensemble& ensemble);

[[nodiscard]] std::size_t apply_binary_breakage(std::vector<Particle>& particles, double step_size, std::mt19937_64& generator,
                                                std::size_t& next_particle_id);

[[nodiscard]] double expected_nucleation_particle_count(double step_size, double particle_weight);

[[nodiscard]] std::size_t apply_nucleation(std::vector<Particle>& particles, double step_size, double particle_weight, std::mt19937_64& generator,
                                           std::size_t& next_particle_id);

[[nodiscard]] std::vector<double> compute_particle_moments(const std::vector<Particle>& particles, double particle_weight, std::size_t maximum_order);

} // namespace qmom_1d_breakage_advection_mc_local_rates
