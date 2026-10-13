#pragma once

#include "monte_carlo/sampled_initialization.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>

namespace qmom_1d_breakage_advection_with_mc
{

constexpr std::size_t initialization_sample_count = 10000;
constexpr unsigned long long initialization_seed = 12345;

using SamplingSetup = qmomkit::monte_carlo::OneDimensionalSampledInitialization;

/**
 * Samples the initial distribution configured in `configuration_file` and
 * estimates the moments of the sample. Both executables start from this sample,
 * so their runs begin from the same finite population. A failure is reported on
 * standard error and returns no value.
 */
[[nodiscard]] std::optional<SamplingSetup> sample_initial_distribution(const std::filesystem::path& configuration_file);

} // namespace qmom_1d_breakage_advection_with_mc
