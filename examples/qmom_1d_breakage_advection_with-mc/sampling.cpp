#include "sampling.hpp"

#include "qmomkit/configuration.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

namespace qmom_1d_breakage_advection_with_mc
{

std::optional<SamplingSetup> sample_initial_distribution(const std::filesystem::path& configuration_file)
{
    std::ifstream input(configuration_file);
    if (!input)
    {
        std::cerr << "Could not open configuration file: " << configuration_file << '\n';
        return std::nullopt;
    }
    std::ostringstream text;
    text << input.rdbuf();

    // The [initial] section describes the distribution analytically; the Monte
    // Carlo companion samples it instead of computing its moments exactly.
    auto initial_condition = qmomkit::config::parse_qmom_initial_condition(text.str(), configuration_file.string());
    if (!initial_condition)
    {
        std::cerr << initial_condition.error() << '\n';
        return std::nullopt;
    }

    const qmomkit::monte_carlo::SamplingOptions sampling_options{initialization_sample_count, initialization_seed};
    std::optional<SamplingSetup> sampled = std::visit(
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
        std::move(initial_condition).value());
    if (sampled && !sampled->ok())
    {
        std::cerr << (sampled->sampling.ok() ? sampled->moments.message : sampled->sampling.message) << '\n';
        return std::nullopt;
    }
    return sampled;
}

} // namespace qmom_1d_breakage_advection_with_mc
