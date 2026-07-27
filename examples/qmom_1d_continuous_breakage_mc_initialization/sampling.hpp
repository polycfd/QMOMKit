#pragma once

#include "qmomkit/configuration.hpp"
#include "monte_carlo/moments.hpp"
#include "monte_carlo/sampling.hpp"

#include <filesystem>
#include <optional>
#include <string>

namespace qmom_1d_continuous_breakage_mc_initialization
{

constexpr std::size_t initialization_sample_count = 10000;
constexpr unsigned long long initialization_seed = 12345;

struct SamplingSetup
{
    qmomkit::monte_carlo::SamplingResult sampling;
    std::size_t maximum_order = 0;
};

[[nodiscard]] std::optional<std::string> read_configuration_text(const std::filesystem::path& path);

[[nodiscard]] std::optional<SamplingSetup> sample_initial_distribution(const qmomkit::config::InitialConditionSpecification& initial_condition);

[[nodiscard]] qmomkit::monte_carlo::OneDimensionalMomentEstimate compute_sampled_moments(const SamplingSetup& sampled);

} // namespace qmom_1d_continuous_breakage_mc_initialization
