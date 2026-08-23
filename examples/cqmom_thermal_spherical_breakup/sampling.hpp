#pragma once

#include "monte_carlo/sampled_initialization.hpp"
#include "qmomkit/configuration.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>

namespace cqmom_thermal_spherical_breakup
{

constexpr std::size_t initialization_sample_count = 10000;
constexpr unsigned long long initialization_seed = 123;

[[nodiscard]] std::optional<std::string> read_configuration_text(const std::filesystem::path& path);
[[nodiscard]] std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
sample_initial_population(const qmomkit::config::InitialConditionSpecification& initial_condition, const qmomkit::PhaseSpaceHierarchy& hierarchy);

} // namespace cqmom_thermal_spherical_breakup
