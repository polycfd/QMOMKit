#pragma once

#include "monte_carlo/birth.hpp"
#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

namespace paper_3d_droplets_condensation_nucleation
{

struct NucleationParameters
{
    double radius = 50.0e-6;
    double vertical_velocity = 0.0;
    double temperature = 298.15;
    double start_time = 0.20;
    double end_time = 0.40;
    double relative_population_increase = 1.0;
};

inline constexpr NucleationParameters nucleation_parameters{};

struct NucleationModel
{
    double flux = 0.0;
};

[[nodiscard]] qmomkit::PhaseSpaceVelocity droplet_phase_space_velocity();

[[nodiscard]] NucleationModel make_nucleation_model(double initial_number_density);

[[nodiscard]] double limit_step_to_nucleation_interval(double time, double maximum_step);

[[nodiscard]] bool nucleation_is_active(double time) noexcept;

[[nodiscard]] qmomkit::CqmomMomentBoundaryContribution nucleation_boundary_contribution(const NucleationModel& model);

[[nodiscard]] qmomkit::monte_carlo::PopulationBirthResult inject_nucleated_droplets(const NucleationModel& model,
                                                                                    const qmomkit::monte_carlo::Ensemble& ensemble,
                                                                                    double number_density, double time_step,
                                                                                    unsigned long long random_seed);

} // namespace paper_3d_droplets_condensation_nucleation
