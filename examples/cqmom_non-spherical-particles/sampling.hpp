#pragma once

#include "monte_carlo/sampled_initialization.hpp"
#include "qmomkit/configuration.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>

namespace cqmom_non_spherical_particles
{

constexpr std::size_t initialization_sample_count = 10000;
constexpr unsigned long long initialization_seed = 12345;

[[nodiscard]] std::optional<std::string> read_configuration_text(const std::filesystem::path& path);

// Samples the independent-product marginals with the fixed example seed and
// returns dense raw moments over maximum_orders from the same finite ensemble.
[[nodiscard]] std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
sample_initial_distribution(const qmomkit::config::InitialConditionSpecification& initial_condition);

} // namespace cqmom_non_spherical_particles
