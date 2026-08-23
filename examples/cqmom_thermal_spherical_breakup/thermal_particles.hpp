#pragma once

#include "monte_carlo/ensemble.hpp"
#include "monte_carlo/sampled_initialization.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <random>
#include <string>
#include <vector>

namespace cqmom_thermal_spherical_breakup
{

constexpr std::size_t initialization_sample_count = 16000;
constexpr unsigned long long initialization_seed = 2026082311ULL;
constexpr unsigned long long breakup_seed = 2026082312ULL;
constexpr double initial_time = 0.0;
constexpr double final_time = 0.02;
constexpr double time_step = 0.001;

struct Diagnostics
{
    double number_density = 0.0;
    double volume = 0.0;
    double volume_temperature = 0.0;
    double momentum_u = 0.0;
    double momentum_v = 0.0;
    double momentum_w = 0.0;
    double volume_u_square = 0.0;
};

[[nodiscard]] std::optional<std::string> read_configuration_text(const std::filesystem::path& path);
[[nodiscard]] std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
sample_initial_population(const qmomkit::config::InitialConditionSpecification& initial_condition, const qmomkit::PhaseSpaceHierarchy& hierarchy);
[[nodiscard]] qmomkit::PhaseSpaceVelocity particle_velocity();
[[nodiscard]] qmomkit::CqmomMomentSource breakup_source(std::size_t& generalized_request_count);
void apply_breakup(qmomkit::monte_carlo::Ensemble& ensemble, double interval, std::mt19937_64& generator);
[[nodiscard]] qmomkit::RawMomentTensor ensemble_moments(const qmomkit::monte_carlo::Ensemble& ensemble, double particle_weight,
                                                        const std::vector<std::size_t>& extents);
[[nodiscard]] double maximum_tensor_difference(const qmomkit::RawMomentTensor& left, const qmomkit::RawMomentTensor& right);
[[nodiscard]] qmomkit::RawMomentRequest diagnostic_request();
[[nodiscard]] Diagnostics diagnostics(const qmomkit::SelectedRawMomentResult& moments);
[[nodiscard]] Diagnostics diagnostics(const qmomkit::RawMomentTensor& moments);
void print_diagnostics(const char* label, const Diagnostics& values);
void ensure_output_directory(const std::filesystem::path& path);

} // namespace cqmom_thermal_spherical_breakup
