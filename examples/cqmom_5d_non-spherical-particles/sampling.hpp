#pragma once

#include "monte_carlo/sampled_initialization.hpp"
#include "qmomkit/distributions.hpp"
#include "qmomkit/phase_space_hierarchy.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <vector>

namespace cqmom_5d_non_spherical_particles
{

constexpr std::size_t initialization_sample_count = 10000;
constexpr unsigned long long initialization_seed = 12345;

// One draw of the shared initial ensemble, with the tensor extents the
// Monte Carlo moment history is reported over.
struct SampledInitialEnsemble
{
    qmomkit::monte_carlo::SamplingResult sampling;
    std::vector<std::size_t> moment_extents;
};

// Samples the initial distribution configured in `configuration_file` with the
// fixed example seed and resolves the tensor extents, without estimating any
// moments: the Monte Carlo driver computes them from the evolving ensemble at
// every output time anyway. A failure is reported on standard error and returns
// no value.
[[nodiscard]] std::optional<SampledInitialEnsemble> sample_initial_ensemble(const std::filesystem::path& configuration_file,
                                                                            const qmomkit::PhaseSpaceHierarchy& hierarchy);

// Samples the same distribution, with the same seed and sample count, but
// estimates only the requested raw moments. The CQMOM run asks for the moments
// its hierarchy requires and never evaluates the rest of the full tensor.
// Sharing the seed and sample count keeps both runs on one realization, so the
// CQMOM and Monte Carlo results remain comparable.
[[nodiscard]] std::optional<qmomkit::monte_carlo::SelectedSampledInitialization>
sample_requested_initial_moments(const std::filesystem::path& configuration_file, const qmomkit::RawMomentRequest& request);

} // namespace cqmom_5d_non_spherical_particles
