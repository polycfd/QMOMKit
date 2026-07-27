#include "sampling.hpp"

#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <type_traits>
#include <variant>

namespace qmom_1d_continuous_breakage_mc_initialization
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

std::optional<SamplingSetup> sample_initial_distribution(const qmomkit::config::InitialConditionSpecification& initial_condition)
{
    const qmomkit::monte_carlo::SamplingOptions sampling_options{initialization_sample_count, initialization_seed};

    return std::visit(
        [&](const auto& specification) -> std::optional<SamplingSetup>
        {
            using Specification = std::decay_t<decltype(specification)>;
            if constexpr (std::is_same_v<Specification, qmomkit::config::DistributionInitialization>)
            {
                return SamplingSetup{
                    qmomkit::monte_carlo::sample_distribution(specification.distribution, specification.number_density, sampling_options),
                    specification.maximum_order};
            }
            else if constexpr (std::is_same_v<Specification, qmomkit::config::MixtureInitialization>)
            {
                return SamplingSetup{qmomkit::monte_carlo::sample_mixture(specification.mixture, sampling_options), specification.maximum_order};
            }
            else
            {
                std::cerr << "Monte Carlo initialization requires a sampleable distribution or mixture configuration.\n";
                return std::nullopt;
            }
        },
        initial_condition);
}

qmomkit::monte_carlo::OneDimensionalMomentEstimate compute_sampled_moments(const SamplingSetup& sampled)
{
    return qmomkit::monte_carlo::compute_ensemble_moments(sampled.sampling.ensemble, sampled.sampling.number_density, sampled.maximum_order);
}

} // namespace qmom_1d_continuous_breakage_mc_initialization
