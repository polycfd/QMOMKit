#pragma once

#include "monte_carlo/sampled_initialization.hpp"
#include "qmomkit/configuration.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>

namespace cqmom_spherical_breakup_recoil
{

constexpr std::size_t initialization_sample_count = 16000;
constexpr unsigned long long initialization_seed = 2026082301ULL;

[[nodiscard]] std::optional<std::string> read_configuration_text(const std::filesystem::path& path);
[[nodiscard]] std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
sample_initial_population(const qmomkit::config::InitialConditionSpecification& initial_condition, const qmomkit::PhaseSpaceHierarchy& hierarchy);

} // namespace cqmom_spherical_breakup_recoil
