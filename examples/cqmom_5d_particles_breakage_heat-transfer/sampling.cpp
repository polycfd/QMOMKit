#include "sampling.hpp"

#include "qmomkit/configuration.hpp"
#include "utilities.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace cqmom_5d_particles_breakage_heat_transfer
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

    std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled = std::visit(
        [&hierarchy](const auto& specification) -> std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
        {
            using Specification = std::decay_t<decltype(specification)>;
            if constexpr (std::is_same_v<Specification, qmomkit::config::IndependentProductInitialization>)
            {
                // Without configured maximum orders, sample the orders the hierarchy needs.
                const std::vector<std::size_t> orders = specification.maximum_orders.value_or(maximum_orders(hierarchy.moment_extents()));
                return qmomkit::monte_carlo::sample_initialization(specification.distribution, orders,
                                                                   {initialization_sample_count, initialization_seed});
            }
            std::cerr << "Monte Carlo initialization requires a sampleable independent-product CQMOM configuration.\n";
            return std::nullopt;
        },
        std::move(initial_condition).value());
    if (sampled && !sampled->ok())
    {
        std::cerr << (sampled->sampling.ok() ? sampled->moments.message : sampled->sampling.message) << '\n';
        return std::nullopt;
    }
    return sampled;
}

} // namespace cqmom_5d_particles_breakage_heat_transfer
