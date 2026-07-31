#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "monte_carlo/sampled_initialization.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/multi_index.hpp"
#include "qmomkit/node_advection.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

namespace
{

constexpr std::size_t sample_count = 4096;
constexpr unsigned long long random_seed = 73001;
constexpr double initial_time = 0.0;
constexpr double final_time = 0.1;
constexpr double output_time_step = 0.05;
constexpr double density_tolerance = 1.0e-12;

enum class Formulation
{
    qmom,
    cqmom
};

struct DriverPaths
{
    std::filesystem::path qmomkit_moments;
    std::filesystem::path monte_carlo_moments;
    std::filesystem::path monte_carlo_statistics;
};

class TextTableWriter
{
public:
    TextTableWriter(const std::filesystem::path& path, const std::vector<std::string>& labels) : output_(path)
    {
        if (!output_)
        {
            throw std::runtime_error("Could not open output file: " + path.string());
        }
        output_ << std::setprecision(17);
        output_ << "# time";
        for (const std::string& label : labels)
        {
            output_ << ' ' << label;
        }
        output_ << '\n';
    }

    void write_row(const double time, const std::vector<double>& values)
    {
        output_ << time;
        for (const double value : values)
        {
            output_ << ' ' << value;
        }
        output_ << '\n';
    }

private:
    std::ofstream output_;
};

class MonteCarloStatisticsWriter
{
public:
    MonteCarloStatisticsWriter(const std::filesystem::path& path, const std::vector<std::string>& moment_labels) : output_(path)
    {
        if (!output_)
        {
            throw std::runtime_error("Could not open output file: " + path.string());
        }
        output_ << std::setprecision(17);
        output_ << "# time sample_count seed number_density";
        for (const std::string& label : moment_labels)
        {
            output_ << " SE_" << label;
        }
        output_ << '\n';
    }

    void write_row(const double time, const std::size_t samples, const unsigned long long seed, const double number_density,
                   const std::vector<double>& standard_errors)
    {
        output_ << time << ' ' << samples << ' ' << seed << ' ' << number_density;
        for (const double standard_error : standard_errors)
        {
            output_ << ' ' << standard_error;
        }
        output_ << '\n';
    }

private:
    std::ofstream output_;
};

[[nodiscard]] std::optional<std::string> read_file(const std::filesystem::path& path)
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

[[nodiscard]] std::string trim(const std::string& input)
{
    const std::size_t first = input.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
    {
        return {};
    }
    const std::size_t last = input.find_last_not_of(" \t\r\n");
    return input.substr(first, last - first + 1);
}

[[nodiscard]] bool parse_positive_size(const std::string& text)
{
    if (text.empty() || text.front() == '-')
    {
        return false;
    }
    try
    {
        std::size_t consumed = 0;
        const unsigned long long value = std::stoull(text, &consumed);
        return consumed == text.size() && value > 0;
    }
    catch (const std::exception&)
    {
        return false;
    }
}

[[nodiscard]] std::optional<Formulation> infer_formulation_from_phase_space(const std::string& text)
{
    std::istringstream input(text);
    std::string section;
    std::string line;
    while (std::getline(input, line))
    {
        const std::size_t comment = line.find_first_of("#;");
        if (comment != std::string::npos)
        {
            line.erase(comment);
        }
        line = trim(line);
        if (line.empty())
        {
            continue;
        }
        if (line.front() == '[' && line.back() == ']')
        {
            section = trim(line.substr(1, line.size() - 2));
            continue;
        }
        const std::size_t equals = line.find('=');
        if (equals == std::string::npos)
        {
            continue;
        }
        const std::string key = trim(line.substr(0, equals));
        const std::string value = trim(line.substr(equals + 1));
        if (section == "phase_space" && key == "node_counts")
        {
            std::istringstream values(value);
            std::string token;
            std::size_t count = 0;
            while (std::getline(values, token, ','))
            {
                if (!parse_positive_size(trim(token)))
                {
                    return std::nullopt;
                }
                ++count;
            }
            if (count == 1 && value.back() != ',')
            {
                return Formulation::qmom;
            }
            if (count > 1 && value.back() != ',')
            {
                return Formulation::cqmom;
            }
            return std::nullopt;
        }
    }
    return std::nullopt;
}

void print_configuration_errors(const std::string& errors)
{
    std::cerr << errors << '\n';
}

[[nodiscard]] qmomkit::AdaptiveRk23Options rk23_options()
{
    qmomkit::AdaptiveRk23Options options;
    options.relative_tolerance = 1.0e-8;
    options.absolute_tolerance = 1.0e-11;
    options.max_step = 0.02;
    return options;
}

[[nodiscard]] qmomkit::PhaseSpaceVelocity validation_velocity()
{
    return [](const double, const std::vector<double>& coordinates)
    {
        std::vector<double> velocity(coordinates.size(), 0.0);
        for (std::size_t dimension = 0; dimension < coordinates.size(); ++dimension)
        {
            velocity[dimension] = 0.04 * static_cast<double>(dimension + 1) + 0.02 * coordinates[dimension];
        }
        return velocity;
    };
}

[[nodiscard]] std::vector<std::string> one_dimensional_labels(const std::size_t moment_count)
{
    std::vector<std::string> labels;
    labels.reserve(moment_count);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        labels.push_back("M" + std::to_string(order));
    }
    return labels;
}

[[nodiscard]] std::string tensor_label(const qmomkit::MultiIndex& index)
{
    std::string label = "M";
    const bool needs_separator = std::any_of(index.begin(), index.end(), [](const std::size_t component) { return component >= 10; });
    for (std::size_t dimension = 0; dimension < index.size(); ++dimension)
    {
        if (needs_separator && dimension > 0)
        {
            label += '_';
        }
        label += std::to_string(index[dimension]);
    }
    return label;
}

[[nodiscard]] std::vector<std::string> tensor_labels(const std::vector<std::size_t>& extents)
{
    const std::vector<qmomkit::MultiIndex> indices = qmomkit::enumerate_multi_indices(extents);
    std::vector<std::string> labels;
    labels.reserve(indices.size());
    for (const qmomkit::MultiIndex& index : indices)
    {
        labels.push_back(tensor_label(index));
    }
    return labels;
}

[[nodiscard]] DriverPaths output_paths(const std::filesystem::path& output_directory, const Formulation formulation)
{
    DriverPaths paths;
    paths.monte_carlo_moments = output_directory / "sample_based_mc_moments.txt";
    paths.monte_carlo_statistics = output_directory / "sample_based_mc_statistics.txt";
    if (formulation == Formulation::qmom)
    {
        paths.qmomkit_moments = output_directory / "sample_based_qmom_moments.txt";
    }
    else
    {
        paths.qmomkit_moments = output_directory / "sample_based_cqmom_moments.txt";
    }
    return paths;
}

[[nodiscard]] bool close_to_number_density(const double value, const double number_density)
{
    return std::abs(value - number_density) <= density_tolerance * std::max(1.0, std::abs(number_density));
}

bool check_m0(const double value, const double number_density, const std::string& label)
{
    if (close_to_number_density(value, number_density))
    {
        return true;
    }
    std::cerr << label << " M0 mismatch: expected " << number_density << ", got " << value << '\n';
    return false;
}

[[nodiscard]] qmomkit::monte_carlo::SamplingOptions sampling_options()
{
    return qmomkit::monte_carlo::SamplingOptions{sample_count, random_seed};
}

std::optional<qmomkit::monte_carlo::OneDimensionalSampledInitialization>
sample_qmom_initial_condition(const qmomkit::config::InitialConditionSpecification& initial_condition)
{
    return std::visit(
        [](const auto& specification) -> std::optional<qmomkit::monte_carlo::OneDimensionalSampledInitialization>
        {
            using Specification = std::decay_t<decltype(specification)>;
            if constexpr (std::is_same_v<Specification, qmomkit::config::DistributionInitialization>)
            {
                return qmomkit::monte_carlo::sample_initialization(specification.distribution, specification.number_density,
                                                                   specification.maximum_order, sampling_options());
            }
            else if constexpr (std::is_same_v<Specification, qmomkit::config::MixtureInitialization>)
            {
                return qmomkit::monte_carlo::sample_initialization(specification.mixture, specification.maximum_order, sampling_options());
            }
            else
            {
                std::cerr << "Sample-based QMOM validation requires a distribution or mixture initial condition.\n";
                std::cerr << "Explicit moments do not define a unique sampleable ensemble.\n";
                return std::nullopt;
            }
        },
        initial_condition);
}

std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
sample_cqmom_initial_condition(const qmomkit::config::InitialConditionSpecification& initial_condition, const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
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
                return qmomkit::monte_carlo::sample_initialization(specification.distribution, maximum_orders, sampling_options());
            }
            else
            {
                std::cerr << "Sample-based CQMOM validation requires an independent-product initial condition.\n";
                std::cerr << "Explicit tensors do not define a unique sampleable ensemble.\n";
                return std::nullopt;
            }
        },
        initial_condition);
}

int run_qmom_case(const std::filesystem::path& configuration_path, const std::string& configuration_text,
                  const std::filesystem::path& output_directory)
{
    const auto configured = qmomkit::config::load_qmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured.error());
        return 1;
    }
    const auto initial_condition = qmomkit::config::parse_qmom_initial_condition(configuration_text, configuration_path.string());
    if (!initial_condition.ok())
    {
        print_configuration_errors(initial_condition.error());
        return 1;
    }

    const std::optional<qmomkit::monte_carlo::OneDimensionalSampledInitialization> sampled = sample_qmom_initial_condition(initial_condition.value());
    if (!sampled.has_value())
    {
        return 1;
    }
    if (!sampled->ok())
    {
        std::cerr << sampled->sampling.message << sampled->moments.message << '\n';
        return 1;
    }

    // The deterministic path starts from moments estimated from the same
    // finite realization used by the Monte Carlo path.
    std::vector<double> qmomkit_moments = sampled->moments.moments.values();
    qmomkit::QmomResult qmom = qmomkit::compute_qmom(sampled->moments.moments, configured.value().node_count, configured.value().options);
    if (!qmom.ok())
    {
        std::cerr << "Initial sample-based QMOM inversion failed: " << qmom.message << '\n';
        return 1;
    }

    const DriverPaths paths = output_paths(output_directory, Formulation::qmom);
    const std::vector<std::string> labels = one_dimensional_labels(qmomkit_moments.size());
    TextTableWriter qmom_output(paths.qmomkit_moments, labels);
    TextTableWriter mc_output(paths.monte_carlo_moments, labels);
    MonteCarloStatisticsWriter stats_output(paths.monte_carlo_statistics, labels);

    const double number_density = sampled->sampling.number_density;
    if (!check_m0(qmomkit_moments.front(), number_density, "QMOM") || !check_m0(sampled->moments.moments[0], number_density, "Monte Carlo"))
    {
        return 1;
    }

    qmom_output.write_row(initial_time, qmomkit_moments);
    mc_output.write_row(initial_time, sampled->moments.moments.values());
    stats_output.write_row(initial_time, sampled->sampling.ensemble.size(), sampled->sampling.seed, number_density, sampled->moments.standard_errors);

    const qmomkit::PhaseSpaceVelocity velocity = validation_velocity();
    qmomkit::NodeAdvectionOptions node_options;
    node_options.rk23_options = rk23_options();
    qmomkit::monte_carlo::Ensemble ensemble = sampled->sampling.ensemble;
    double time = initial_time;

    while (time < final_time)
    {
        const double step = (time + output_time_step > final_time) ? (final_time - time) : output_time_step;
        const qmomkit::QmomNodeAdvectionResult advanced_qmom = qmomkit::advance_nodes(qmom, time, step, velocity, node_options);
        if (!advanced_qmom.ok())
        {
            std::cerr << "QMOM node advection failed: " << advanced_qmom.message << '\n';
            return 1;
        }
        qmom = qmomkit::QmomResult{advanced_qmom.quadrature.abscissas, advanced_qmom.quadrature.weights, qmomkit::QmomStatus::success, ""};
        // Compare raw moments after each continuous node/sample update; node
        // coordinates themselves are diagnostics, not the validation quantity.
        qmomkit_moments = qmomkit::reconstruct_raw_moments(advanced_qmom.quadrature, qmomkit_moments.size());

        const qmomkit::monte_carlo::EnsembleEvolutionResult advanced_mc =
            qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + step, node_options.rk23_options);
        if (!advanced_mc.ok())
        {
            std::cerr << "Monte Carlo ensemble evolution failed: " << advanced_mc.message << '\n';
            return 1;
        }
        ensemble = advanced_mc.ensemble;
        const qmomkit::monte_carlo::OneDimensionalMomentEstimate mc_moments =
            qmomkit::monte_carlo::compute_ensemble_moments(ensemble, number_density, qmomkit_moments.size() - 1);
        if (!mc_moments.ok())
        {
            std::cerr << "Monte Carlo moment extraction failed: " << mc_moments.message << '\n';
            return 1;
        }

        time += step;
        if (!check_m0(qmomkit_moments.front(), number_density, "QMOM") || !check_m0(mc_moments.moments[0], number_density, "Monte Carlo"))
        {
            return 1;
        }
        qmom_output.write_row(time, qmomkit_moments);
        mc_output.write_row(time, mc_moments.moments.values());
        stats_output.write_row(time, ensemble.size(), sampled->sampling.seed, number_density, mc_moments.standard_errors);
    }

    std::cerr << "Wrote sample-based QMOM history to " << paths.qmomkit_moments << '\n';
    std::cerr << "Wrote Monte Carlo history to " << paths.monte_carlo_moments << '\n';
    std::cerr << "Wrote Monte Carlo statistics to " << paths.monte_carlo_statistics << '\n';
    return 0;
}

int run_cqmom_case(const std::filesystem::path& configuration_path, const std::string& configuration_text,
                   const std::filesystem::path& output_directory)
{
    const auto configured = qmomkit::config::load_cqmom_dense(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured.error());
        return 1;
    }
    const auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(configuration_text, configuration_path.string());
    if (!initial_condition.ok())
    {
        print_configuration_errors(initial_condition.error());
        return 1;
    }

    const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
        sample_cqmom_initial_condition(initial_condition.value(), configured.value().hierarchy);
    if (!sampled.has_value())
    {
        return 1;
    }
    if (!sampled->ok())
    {
        std::cerr << sampled->sampling.message << sampled->moments.message << '\n';
        return 1;
    }

    // The deterministic path starts from the tensor moments estimated from the
    // same finite realization used by the Monte Carlo path.
    qmomkit::RawMomentTensor qmomkit_moments = sampled->moments.moments;
    qmomkit::CqmomDenseCase dense_case = configured.value();
    dense_case.moments = qmomkit_moments;
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(dense_case);
    if (!initial.ok())
    {
        std::cerr << "Initial sample-based CQMOM inversion failed: " << initial.message << '\n';
        return 1;
    }
    qmomkit::CqmomQuadrature cqmom = initial.quadrature;

    const DriverPaths paths = output_paths(output_directory, Formulation::cqmom);
    const std::vector<std::string> labels = tensor_labels(qmomkit_moments.extents);
    TextTableWriter cqmom_output(paths.qmomkit_moments, labels);
    TextTableWriter mc_output(paths.monte_carlo_moments, labels);
    MonteCarloStatisticsWriter stats_output(paths.monte_carlo_statistics, labels);

    const double number_density = sampled->sampling.number_density;
    if (!check_m0(qmomkit_moments.values.front(), number_density, "CQMOM") ||
        !check_m0(sampled->moments.moments.values.front(), number_density, "Monte Carlo"))
    {
        return 1;
    }

    cqmom_output.write_row(initial_time, qmomkit_moments.values);
    mc_output.write_row(initial_time, sampled->moments.moments.values);
    stats_output.write_row(initial_time, sampled->sampling.ensemble.size(), sampled->sampling.seed, number_density, sampled->moments.standard_errors);

    const qmomkit::PhaseSpaceVelocity velocity = validation_velocity();
    qmomkit::NodeAdvectionOptions node_options;
    node_options.rk23_options = rk23_options();
    qmomkit::monte_carlo::Ensemble ensemble = sampled->sampling.ensemble;
    double time = initial_time;

    while (time < final_time)
    {
        const double step = (time + output_time_step > final_time) ? (final_time - time) : output_time_step;
        const qmomkit::CqmomNodeAdvectionResult advanced_cqmom = qmomkit::advance_nodes(cqmom, time, step, velocity, node_options);
        if (!advanced_cqmom.ok())
        {
            std::cerr << "CQMOM node advection failed: " << advanced_cqmom.message << '\n';
            return 1;
        }
        cqmom = advanced_cqmom.quadrature;
        // Compare reconstructed tensor moments after each continuous update;
        // CQMOM branch structure is diagnostic rather than a comparison target.
        qmomkit_moments = qmomkit::reconstruct_raw_moments(cqmom, qmomkit_moments.extents);

        const qmomkit::monte_carlo::EnsembleEvolutionResult advanced_mc =
            qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + step, node_options.rk23_options);
        if (!advanced_mc.ok())
        {
            std::cerr << "Monte Carlo ensemble evolution failed: " << advanced_mc.message << '\n';
            return 1;
        }
        ensemble = advanced_mc.ensemble;
        std::vector<std::size_t> maximum_orders(qmomkit_moments.extents.size(), 0);
        for (std::size_t dimension = 0; dimension < qmomkit_moments.extents.size(); ++dimension)
        {
            maximum_orders[dimension] = qmomkit_moments.extents[dimension] - 1;
        }
        const qmomkit::monte_carlo::TensorMomentEstimate mc_moments =
            qmomkit::monte_carlo::compute_ensemble_moments(ensemble, number_density, maximum_orders);
        if (!mc_moments.ok())
        {
            std::cerr << "Monte Carlo tensor moment extraction failed: " << mc_moments.message << '\n';
            return 1;
        }

        time += step;
        if (!check_m0(qmomkit_moments.values.front(), number_density, "CQMOM") ||
            !check_m0(mc_moments.moments.values.front(), number_density, "Monte Carlo"))
        {
            return 1;
        }
        cqmom_output.write_row(time, qmomkit_moments.values);
        mc_output.write_row(time, mc_moments.moments.values);
        stats_output.write_row(time, ensemble.size(), sampled->sampling.seed, number_density, mc_moments.standard_errors);
    }

    std::cerr << "Wrote sample-based CQMOM history to " << paths.qmomkit_moments << '\n';
    std::cerr << "Wrote Monte Carlo history to " << paths.monte_carlo_moments << '\n';
    std::cerr << "Wrote Monte Carlo statistics to " << paths.monte_carlo_statistics << '\n';
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    const std::filesystem::path configuration_path = (argc > 1) ? argv[1] : "configuration.qmom";
    const std::filesystem::path output_directory = (argc > 2) ? argv[2] : ".";
    std::filesystem::create_directories(output_directory);

    const std::optional<std::string> configuration_text = read_file(configuration_path);
    if (!configuration_text.has_value())
    {
        return 1;
    }
    const std::optional<Formulation> formulation = infer_formulation_from_phase_space(*configuration_text);
    if (!formulation.has_value())
    {
        std::cerr << "Configuration phase_space.node_counts must define one positive dimension for QMOM or multiple positive dimensions for CQMOM.\n";
        return 1;
    }

    try
    {
        if (*formulation == Formulation::qmom)
        {
            return run_qmom_case(configuration_path, *configuration_text, output_directory);
        }
        return run_cqmom_case(configuration_path, *configuration_text, output_directory);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
