#pragma once

#include "monte_carlo/sampled_initialization.hpp"
#include "qmomkit/configuration.hpp"

#include <filesystem>
#include <optional>
#include <string>

namespace qmom_1d_breakage_advection_with_mc
{

constexpr std::size_t initialization_sample_count = 10000;
constexpr unsigned long long initialization_seed = 12345;

using SamplingSetup = qmomkit::monte_carlo::OneDimensionalSampledInitialization;

[[nodiscard]] std::optional<std::string> read_configuration_text(const std::filesystem::path& path);

[[nodiscard]] std::optional<SamplingSetup> sample_initial_distribution(const qmomkit::config::InitialConditionSpecification& initial_condition);

} // namespace qmom_1d_breakage_advection_with_mc
