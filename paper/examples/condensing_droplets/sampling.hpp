#pragma once

#include "monte_carlo/sampled_initialization.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/phase_space_hierarchy.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>

namespace cqmom_condensing_droplets
{

constexpr std::size_t initialization_sample_count = 10000;
constexpr unsigned long long initialization_seed = 20260199;

[[nodiscard]] std::optional<std::string> read_configuration_text(const std::filesystem::path& path);

[[nodiscard]] std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
sample_initial_distribution(const qmomkit::config::InitialConditionSpecification& initial_condition, const qmomkit::PhaseSpaceHierarchy& hierarchy);

} // namespace cqmom_condensing_droplets
