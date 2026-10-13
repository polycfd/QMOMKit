#pragma once

#include "monte_carlo/sampled_initialization.hpp"
#include "qmomkit/phase_space_hierarchy.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>

namespace cqmom_2d_bubbles_cavitation
{

constexpr std::size_t initialization_sample_count = 10000;
constexpr unsigned long long initialization_seed = 12345;

/**
 * Samples the initial distribution configured in `configuration_file` and
 * estimates the moment tensor the hierarchy needs from the sample. Both
 * executables start from this sample, so their runs begin from the same finite
 * population. A failure is reported on standard error and returns no value.
 */
[[nodiscard]] std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
sample_initial_distribution(const std::filesystem::path& configuration_file, const qmomkit::PhaseSpaceHierarchy& hierarchy);

} // namespace cqmom_2d_bubbles_cavitation
