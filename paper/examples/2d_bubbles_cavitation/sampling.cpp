#include "sampling.hpp"

#include "qmomkit/configuration.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace paper_2d_bubbles_cavitation
{

std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sample_initial_distribution(const std::filesystem::path& configuration_file,
                                                                                             const qmomkit::PhaseSpaceHierarchy& hierarchy)
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
    auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(text.str(), configuration_file.string());
    if (!initial_condition)
    {
        std::cerr << initial_condition.error() << '\n';
        return std::nullopt;
    }

    const qmomkit::monte_carlo::SamplingOptions sampling_options{initialization_sample_count, initialization_seed};
    std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled = std::visit(
        [&](const auto& specification) -> std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
        {
            using Specification = std::decay_t<decltype(specification)>;
            if constexpr (std::is_same_v<Specification, qmomkit::config::IndependentProductInitialization>)
            {
                // Without configured maximum orders, sample the orders the hierarchy needs.
                std::vector<std::size_t> maximum_orders;
                if (specification.maximum_orders.has_value())
                {
                    maximum_orders = *specification.maximum_orders;
                }
                else
                {
                    for (const std::size_t extent : hierarchy.moment_extents())
                    {
                        maximum_orders.push_back(extent - 1);
                    }
                }
                return qmomkit::monte_carlo::sample_initialization(specification.distribution, maximum_orders, sampling_options);
            }
            else
            {
                std::cerr << "Monte Carlo initialization requires a sampleable independent-product CQMOM configuration.\n";
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

} // namespace paper_2d_bubbles_cavitation
