#include "sampling.hpp"

#include "qmomkit/configuration.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>
#include <variant>

namespace cqmom_5d_non_spherical_particles
{
namespace
{

// Both entry points below sample the distribution of the [initial] section and
// differ only in which moments they estimate from the draw.
[[nodiscard]] std::optional<qmomkit::config::IndependentProductInitialization>
configured_distribution(const std::filesystem::path& configuration_file)
{
    std::ifstream input(configuration_file);
    if (!input)
    {
        std::cerr << "Could not open configuration file: " << configuration_file << '\n';
        return std::nullopt;
    }
    std::ostringstream text;
    text << input.rdbuf();

    auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(text.str(), configuration_file.string());
    if (!initial_condition)
    {
        std::cerr << initial_condition.error() << '\n';
        return std::nullopt;
    }
    auto* const specification = std::get_if<qmomkit::config::IndependentProductInitialization>(&initial_condition.value());
    if (specification == nullptr)
    {
        std::cerr << "Monte Carlo initialization requires a sampleable independent-product CQMOM configuration.\n";
        return std::nullopt;
    }
    return std::move(*specification);
}

} // namespace

std::optional<SampledInitialEnsemble> sample_initial_ensemble(const std::filesystem::path& configuration_file,
                                                              const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    const std::optional<qmomkit::config::IndependentProductInitialization> specification = configured_distribution(configuration_file);
    if (!specification)
    {
        return std::nullopt;
    }

    // Configured maximum orders are inclusive, so the reported extents are one
    // larger; without them the hierarchy's own extents are already the 2*N span.
    std::vector<std::size_t> moment_extents;
    if (specification->maximum_orders.has_value())
    {
        for (const std::size_t maximum_order : *specification->maximum_orders)
        {
            moment_extents.push_back(maximum_order + 1);
        }
    }
    else
    {
        moment_extents = hierarchy.moment_extents();
    }

    SampledInitialEnsemble sampled{
        qmomkit::monte_carlo::sample_independent_product(specification->distribution, {initialization_sample_count, initialization_seed}),
        std::move(moment_extents)};
    if (!sampled.sampling.ok())
    {
        std::cerr << sampled.sampling.message << '\n';
        return std::nullopt;
    }
    return sampled;
}

std::optional<qmomkit::monte_carlo::SelectedSampledInitialization> sample_requested_initial_moments(const std::filesystem::path& configuration_file,
                                                                                                    const qmomkit::RawMomentRequest& request)
{
    const std::optional<qmomkit::config::IndependentProductInitialization> specification = configured_distribution(configuration_file);
    if (!specification)
    {
        return std::nullopt;
    }

    // maximum_orders only sizes a full tensor; the requested indices come from the
    // CQMOM hierarchy instead, so it has no effect on this path.
    qmomkit::monte_carlo::SelectedSampledInitialization sampled =
        qmomkit::monte_carlo::sample_initialization(specification->distribution, request, {initialization_sample_count, initialization_seed});
    if (!sampled.ok())
    {
        std::cerr << (sampled.sampling.ok() ? sampled.moments.message : sampled.sampling.message) << '\n';
        return std::nullopt;
    }
    return sampled;
}

} // namespace cqmom_5d_non_spherical_particles
