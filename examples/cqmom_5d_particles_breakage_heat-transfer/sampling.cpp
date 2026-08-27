#include "sampling.hpp"

#include "utilities.hpp"

#include <fstream>
#include <iostream>
#include <type_traits>
#include <variant>

namespace cqmom_5d_particles_breakage_heat_transfer
{

std::optional<std::string> read_configuration_text(const std::filesystem::path& path)
{
    std::ifstream input(path);
    if (!input)
    {
        std::cerr << "Could not open configuration file: " << path << '\n';
        return std::nullopt;
    }
    return std::string(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
}

std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
sample_initial_population(const qmomkit::config::InitialConditionSpecification& initial_condition, const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    return std::visit(
        [&hierarchy](const auto& specification) -> std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
        {
            using Specification = std::decay_t<decltype(specification)>;
            if constexpr (std::is_same_v<Specification, qmomkit::config::IndependentProductInitialization>)
            {
                const std::vector<std::size_t> orders = specification.maximum_orders.value_or(maximum_orders(hierarchy.moment_extents()));
                return qmomkit::monte_carlo::sample_initialization(specification.distribution, orders,
                                                                   {initialization_sample_count, initialization_seed});
            }
            std::cerr << "This example requires a sampleable independent-product initial condition.\n";
            return std::nullopt;
        },
        initial_condition);
}

} // namespace cqmom_5d_particles_breakage_heat_transfer
