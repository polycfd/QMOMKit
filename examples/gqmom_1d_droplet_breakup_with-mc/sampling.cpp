#include "sampling.hpp"

#include "qmomkit/configuration.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <variant>

namespace gqmom_1d_droplet_breakup_with_mc
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

    auto initial_condition = qmomkit::config::parse_qmom_initial_condition(text.str(), configuration_file.string());
    if (!initial_condition)
    {
        std::cerr << initial_condition.error() << '\n';
        return std::nullopt;
    }

    // The configuration defines one distribution, which the Monte Carlo
    // companion samples instead of computing its moments exactly.
    const auto& distribution = std::get<qmomkit::config::DistributionInitialization>(initial_condition.value());
    const qmomkit::monte_carlo::SamplingOptions sampling_options{initialization_sample_count, initialization_seed};
    SamplingSetup sampled = qmomkit::monte_carlo::sample_initialization(distribution.distribution, distribution.number_density,
                                                                        distribution.maximum_order, sampling_options);
    if (!sampled.ok())
    {
        std::cerr << (sampled.sampling.ok() ? sampled.moments.message : sampled.sampling.message) << '\n';
        return std::nullopt;
    }
    return sampled;
}

} // namespace gqmom_1d_droplet_breakup_with_mc
