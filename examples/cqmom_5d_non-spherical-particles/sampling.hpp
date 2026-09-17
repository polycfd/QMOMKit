#pragma once

#include "monte_carlo/sampled_initialization.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/phase_space_hierarchy.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace cqmom_5d_non_spherical_particles
{

constexpr std::size_t initialization_sample_count = 10000;
constexpr unsigned long long initialization_seed = 12345;

[[nodiscard]] std::optional<std::string> read_configuration_text(const std::filesystem::path& path);

// One draw of the shared initial ensemble, with the dense moment extents the
// Monte Carlo moment history is reported over.
struct SampledInitialEnsemble
{
    qmomkit::monte_carlo::SamplingResult sampling;
    std::vector<std::size_t> moment_extents;
};

// Samples the independent-product marginals with the fixed example seed and
// resolves the dense moment extents, without estimating any moments: the Monte
// Carlo driver recomputes dense moments from the evolving ensemble at every
// output time anyway, so estimating them once more here would be discarded work.
[[nodiscard]] std::optional<SampledInitialEnsemble> sample_initial_ensemble(const qmomkit::config::InitialConditionSpecification& initial_condition,
                                                                            const qmomkit::PhaseSpaceHierarchy& hierarchy);

// Samples the same marginals, with the same fixed seed and sample count, but
// estimates only the requested raw moments. The CQMOM run asks for the moments
// its hierarchy requires and never evaluates the rest of the dense rectangle.
// Sharing the seed and sample count keeps both runs on one realization, so the
// CQMOM and Monte Carlo results remain comparable.
[[nodiscard]] std::optional<qmomkit::monte_carlo::SelectedSampledInitialization>
sample_requested_initial_moments(const qmomkit::config::InitialConditionSpecification& initial_condition, const qmomkit::RawMomentRequest& request);

} // namespace cqmom_5d_non_spherical_particles
