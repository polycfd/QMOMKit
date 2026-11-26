#include "sampling.hpp"

#include "qmomkit/configuration.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <variant>

namespace cqmom_2d_hyperbolic_drag_dispersion
{

std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sample_initial_distribution(const std::filesystem::path& configuration_file)
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
    const auto* specification = std::get_if<qmomkit::config::IndependentProductInitialization>(&initial_condition.value());
    if (specification == nullptr || !specification->maximum_orders.has_value())
    {
        std::cerr << "Monte Carlo initialization requires an independent-product [initial] section with initial.maximum_orders.\n";
        return std::nullopt;
    }

    const qmomkit::monte_carlo::SamplingOptions sampling_options{initialization_sample_count, initialization_seed};
    qmomkit::monte_carlo::TensorSampledInitialization sampled =
        qmomkit::monte_carlo::sample_initialization(specification->distribution, *specification->maximum_orders, sampling_options);
    if (!sampled.ok())
    {
        std::cerr << (sampled.sampling.ok() ? sampled.moments.message : sampled.sampling.message) << '\n';
        return std::nullopt;
    }
    return sampled;
}

} // namespace cqmom_2d_hyperbolic_drag_dispersion
