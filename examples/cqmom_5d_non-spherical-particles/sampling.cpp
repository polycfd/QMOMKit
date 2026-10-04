#include "sampling.hpp"

#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <variant>
#include <vector>

namespace cqmom_5d_non_spherical_particles
{

std::optional<std::string> read_configuration_text(const std::filesystem::path& path)
{
    std::ifstream input(path);
    if (!input)
    {
        std::cerr << "Could not open configuration file: " << path << '\n';
        return std::nullopt;
    }

    std::ostringstream contents;
    contents << input.rdbuf();
    return contents.str();
}

namespace
{

// Both entry points below sample the same thing and differ only in which
// moments they estimate from the draw, so the configuration check and the
// sampling options live here once.
[[nodiscard]] const qmomkit::config::IndependentProductInitialization*
sampleable_specification(const qmomkit::config::InitialConditionSpecification& initial_condition)
{
    const auto* const specification = std::get_if<qmomkit::config::IndependentProductInitialization>(&initial_condition);
    if (specification == nullptr)
    {
        std::cerr << "Monte Carlo initialization requires a sampleable independent-product CQMOM configuration.\n";
    }
    return specification;
}

[[nodiscard]] qmomkit::monte_carlo::SamplingOptions example_sampling_options()
{
    return {initialization_sample_count, initialization_seed};
}

} // namespace

std::optional<SampledInitialEnsemble> sample_initial_ensemble(const qmomkit::config::InitialConditionSpecification& initial_condition,
                                                              const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    const qmomkit::config::IndependentProductInitialization* const specification = sampleable_specification(initial_condition);
    if (specification == nullptr)
    {
        return std::nullopt;
    }

    // Configured maximum orders are inclusive, so the reported extents are one
    // larger; without them the hierarchy's own extents are already the 2*N span.
    std::vector<std::size_t> moment_extents;
    if (specification->maximum_orders.has_value())
    {
        moment_extents.reserve(specification->maximum_orders->size());
        for (const std::size_t maximum_order : *specification->maximum_orders)
        {
            moment_extents.push_back(maximum_order + 1);
        }
    }
    else
    {
        moment_extents = hierarchy.moment_extents();
    }

    return SampledInitialEnsemble{qmomkit::monte_carlo::sample_independent_product(specification->distribution, example_sampling_options()),
                                  std::move(moment_extents)};
}

std::optional<qmomkit::monte_carlo::SelectedSampledInitialization>
sample_requested_initial_moments(const qmomkit::config::InitialConditionSpecification& initial_condition, const qmomkit::RawMomentRequest& request)
{
    const qmomkit::config::IndependentProductInitialization* const specification = sampleable_specification(initial_condition);
    if (specification == nullptr)
    {
        return std::nullopt;
    }

    // maximum_orders only sizes a full tensor; the requested indices come from the
    // CQMOM hierarchy instead, so it has no effect on this path.
    return qmomkit::monte_carlo::sample_initialization(specification->distribution, request, example_sampling_options());
}

} // namespace cqmom_5d_non_spherical_particles
