#include "sampling.hpp"

#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <type_traits>
#include <variant>

namespace qmom_1d_breakage_advection_mc_initialization
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
                return qmomkit::monte_carlo::sample_initialization(specification.distribution, specification.number_density,
                                                                   specification.maximum_order, sampling_options);
            }
            else if constexpr (std::is_same_v<Specification, qmomkit::config::MixtureInitialization>)
            {
                return qmomkit::monte_carlo::sample_initialization(specification.mixture, specification.maximum_order, sampling_options);
            }
            else
            {
                std::cerr << "Monte Carlo initialization requires a sampleable distribution or mixture configuration.\n";
                return std::nullopt;
            }
        },
        initial_condition);
}

} // namespace qmom_1d_breakage_advection_mc_initialization
