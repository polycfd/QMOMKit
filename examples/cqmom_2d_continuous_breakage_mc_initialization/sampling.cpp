#include "sampling.hpp"

#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <type_traits>
#include <variant>

namespace cqmom_2d_continuous_breakage
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

std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
sample_initial_distribution(const qmomkit::config::InitialConditionSpecification& initial_condition)
{
    const qmomkit::monte_carlo::SamplingOptions sampling_options{initialization_sample_count, initialization_seed};

    return std::visit(
        [&](const auto& specification) -> std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
        {
            using Specification = std::decay_t<decltype(specification)>;
            if constexpr (std::is_same_v<Specification, qmomkit::config::IndependentProductInitialization>)
            {
                if (!specification.maximum_orders.has_value())
                {
                    std::cerr << "Monte Carlo initialization requires initial.maximum_orders for dense tensor moment sampling.\n";
                    return std::nullopt;
                }
                return qmomkit::monte_carlo::sample_initialization(specification.distribution, *specification.maximum_orders, sampling_options);
            }
            else
            {
                std::cerr << "Monte Carlo initialization requires a sampleable independent-product CQMOM configuration.\n";
                return std::nullopt;
            }
        },
        initial_condition);
}

} // namespace cqmom_2d_continuous_breakage
