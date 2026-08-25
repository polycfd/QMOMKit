#include "sampling.hpp"

#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <type_traits>
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

std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
sample_initial_distribution(const qmomkit::config::InitialConditionSpecification& initial_condition, const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    const qmomkit::monte_carlo::SamplingOptions sampling_options{initialization_sample_count, initialization_seed};

    return std::visit(
        [&](const auto& specification) -> std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
        {
            using Specification = std::decay_t<decltype(specification)>;
            if constexpr (std::is_same_v<Specification, qmomkit::config::IndependentProductInitialization>)
            {
                std::vector<std::size_t> maximum_orders;
                if (specification.maximum_orders.has_value())
                {
                    maximum_orders = *specification.maximum_orders;
                }
                else
                {
                    const std::vector<std::size_t> extents = hierarchy.moment_extents();
                    maximum_orders.reserve(extents.size());
                    for (const std::size_t extent : extents)
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
        initial_condition);
}

} // namespace cqmom_5d_non_spherical_particles
