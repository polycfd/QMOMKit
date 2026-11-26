#pragma once

#include "monte_carlo/sampled_initialization.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>

namespace cqmom_2d_hyperbolic_drag_dispersion
{

constexpr std::size_t initialization_sample_count = 10000;
constexpr unsigned long long initialization_seed = 12345;

/**
 * Samples the independent-product initial distribution configured in
 * `configuration_file` and estimates its moment tensor from the sample, up to
 * the orders that `initial.maximum_orders` gives. All three executables start
 * from this sample, so their runs begin from the same finite population. A
 * failure is reported on standard error and returns no value.
 */
[[nodiscard]] std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
sample_initial_distribution(const std::filesystem::path& configuration_file);

} // namespace cqmom_2d_hyperbolic_drag_dispersion
