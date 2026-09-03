#include "qmomkit/configuration.hpp"
#include "qmomkit/moment_reconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace
{

constexpr double tolerance = 2.0e-10;

struct Terminal
{
    std::vector<double> coordinates;
    double weight;
    qmomkit::ConditioningPath path;
};

struct Errors
{
    double node = 0.0;
    double weight = 0.0;
    double moment = 0.0;
};

std::string path_key(const qmomkit::ConditioningPath& path)
{
    std::string key;
    for (std::size_t level = 0; level < path.size(); ++level)
    {
        if (level != 0)
        {
            key += ',';
        }
        key += std::to_string(path.node_index_at_level(level));
    }
    return key;
}

double relative_error(const double actual, const double expected)
{
    return std::abs(actual - expected) / std::max(1.0, std::abs(expected));
}

std::vector<Terminal> physical_reference(const std::vector<double>& units)
{
    // This five-terminal rule is deliberately ragged in conditioning order
    // (xi2, xi0, xi1): root xi2 has two nodes; its xi0 children have one and
    // two nodes; the final xi1 children have two, one, and two nodes. Raw
    // moments are prescribed independently in the .qmom fixtures as
    // M_alpha = sum W xi0^a0 xi1^a1 xi2^a2.  The rescaled fixture uses
    // xi' = (1e6 xi0, 1e-6 xi1, 1e3 xi2), so its mixed moments obey
    // M'_alpha = 1e(6 a0 - 6 a1 + 3 a2) M_alpha irrespective of recursive
    // conditioning order.  Comparing recovered coordinates after division by
    // these physical unit factors makes escaped normalized coordinates obvious.
    const std::vector<Terminal> canonical{{{0.5, -2.0, -1.0}, 0.10, qmomkit::ConditioningPath({0, 0, 0})},
                                          {{0.5, 1.0, -1.0}, 0.30, qmomkit::ConditioningPath({0, 0, 1})},
                                          {{-1.0, 0.25, 2.0}, 0.18, qmomkit::ConditioningPath({1, 0, 0})},
                                          {{1.5, -0.5, 2.0}, 0.168, qmomkit::ConditioningPath({1, 1, 0})},
                                          {{1.5, 2.0, 2.0}, 0.252, qmomkit::ConditioningPath({1, 1, 1})}};
    std::vector<Terminal> scaled = canonical;
    for (Terminal& terminal : scaled)
    {
        for (std::size_t dimension = 0; dimension < units.size(); ++dimension)
        {
            terminal.coordinates[dimension] *= units[dimension];
        }
    }
    return scaled;
}

std::vector<Terminal> recovered_terminals(const qmomkit::CqmomQuadrature& quadrature)
{
    std::vector<Terminal> terminals;
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        std::vector<double> coordinates(quadrature.number_of_dimensions(), 0.0);
        for (std::size_t dimension = 0; dimension < coordinates.size(); ++dimension)
        {
            coordinates[dimension] = quadrature.abscissas()[dimension][node];
        }
        terminals.push_back({std::move(coordinates), quadrature.weights()[node], quadrature.terminal_paths()[node]});
    }
    std::sort(terminals.begin(), terminals.end(), [](const Terminal& lhs, const Terminal& rhs) { return lhs.path < rhs.path; });
    return terminals;
}

void check_result(const std::string& route, const qmomkit::CqmomResult& result, const qmomkit::RawMomentTensor& configured,
                  const std::vector<double>& units, Errors& errors, bool& passed)
{
    if (!result.ok())
    {
        std::cerr << route << " reconstruction failed: " << result.message << '\n';
        passed = false;
        return;
    }
    if (result.quadrature.has_rectangular_storage() || !result.quadrature.adaptive_branch_counts_enabled())
    {
        std::cerr << route << " did not return ragged adaptive CQMOM storage.\n";
        passed = false;
    }

    const std::vector<Terminal> expected = physical_reference(units);
    const std::vector<Terminal> recovered = recovered_terminals(result.quadrature);
    if (recovered.size() != expected.size())
    {
        std::cerr << route << " terminal count differs from the prescribed ragged rule.\n";
        passed = false;
        return;
    }
    for (std::size_t terminal = 0; terminal < expected.size(); ++terminal)
    {
        if (recovered[terminal].path != expected[terminal].path)
        {
            std::cerr << route << " terminal path differs from the prescribed ragged rule.\n";
            passed = false;
        }
        errors.weight = std::max(errors.weight, std::abs(recovered[terminal].weight - expected[terminal].weight));
        for (std::size_t dimension = 0; dimension < units.size(); ++dimension)
        {
            errors.node = std::max(errors.node, relative_error(recovered[terminal].coordinates[dimension] / units[dimension],
                                                               expected[terminal].coordinates[dimension] / units[dimension]));
        }
    }

    const std::map<std::string, std::size_t> expected_counts{{"0:", 2}, {"1:0", 1}, {"1:1", 2}, {"2:0,0", 2}, {"2:1,0", 1}, {"2:1,1", 2}};
    std::map<std::string, std::size_t> actual_counts;
    for (const qmomkit::CqmomDiagnostics& diagnostic : result.diagnostics)
    {
        actual_counts.emplace(std::to_string(diagnostic.conditioning_level) + ":" + path_key(diagnostic.parent_path), diagnostic.actual_node_count);
        if (diagnostic.phase_space_dimension != result.quadrature.hierarchy().conditioning_order()[diagnostic.conditioning_level] ||
            diagnostic.requested_node_count != 2 || !diagnostic.adaptive_branch_counts_enabled)
        {
            std::cerr << route << " diagnostics do not retain the configured recursive topology.\n";
            passed = false;
        }
    }
    if (actual_counts != expected_counts)
    {
        std::cerr << route << " active branch counts differ from the prescribed ragged topology.\n";
        passed = false;
    }

    const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(result.quadrature, configured.extents);
    for (std::size_t index = 0; index < configured.values.size(); ++index)
    {
        errors.moment = std::max(errors.moment, relative_error(reconstructed.values[index], configured.values[index]));
    }
}

} // namespace

int main(const int argc, char** argv)
{
    if (argc != 3)
    {
        std::cerr << "Usage: automatic_scaling <configuration.qmom> <reference|rescaled>\n";
        return 2;
    }

    const std::string representation = argv[2];
    const std::vector<double> units = representation == "reference"  ? std::vector<double>{1.0, 1.0, 1.0}
                                      : representation == "rescaled" ? std::vector<double>{1.0e6, 1.0e-6, 1.0e3}
                                                                     : std::vector<double>{};
    if (units.empty())
    {
        std::cerr << "Unknown unit representation '" << representation << "'.\n";
        return 2;
    }

    const std::filesystem::path configuration_path = argv[1];
    const auto dense_configured = qmomkit::config::load_cqmom_dense(configuration_path);
    const auto lean_configured = qmomkit::config::load_cqmom_lean(configuration_path);
    if (!dense_configured.ok() || !lean_configured.ok())
    {
        std::cerr << (!dense_configured.ok() ? dense_configured.error() : lean_configured.error()) << '\n';
        return 1;
    }

    const qmomkit::CqmomDenseCase& dense_case = dense_configured.value();
    const qmomkit::CqmomLeanCase& lean_case = lean_configured.value();
    bool passed = dense_case.hierarchy.conditioning_order() == std::vector<std::size_t>({2, 0, 1}) && dense_case.adaptive_branch_counts &&
                  dense_case.qmom_options.adaptive && dense_case.qmom_options.automatic_scaling && lean_case.qmom_options.automatic_scaling;
    if (!passed)
    {
        std::cerr << "The configured automatic-scaling, adaptive, or non-identity conditioning policy was not propagated.\n";
    }

    // The old layout test compared two freshly-created layouts without a CQMOM
    // workflow.  Here both configured routes actually reconstruct the scaled
    // physical rule, while their case-level layouts are also compared to an
    // otherwise identical explicit automatic_scaling=false case.
    const qmomkit::CqmomMomentLayout required_layout = qmomkit::create_cqmom_moment_layout(dense_case.hierarchy);
    qmomkit::CqmomLeanCase unscaled_case = lean_case;
    unscaled_case.qmom_options.automatic_scaling = false;
    if (lean_case.moments.layout.required_indices != required_layout.required_indices ||
        unscaled_case.moments.layout.required_indices != lean_case.moments.layout.required_indices)
    {
        std::cerr << "automatic scaling changed the required CQMOM moment layout.\n";
        passed = false;
    }

    Errors errors;
    check_result("dense", qmomkit::compute_cqmom(dense_case), dense_case.moments, units, errors, passed);
    check_result("lean", qmomkit::compute_cqmom(lean_case), dense_case.moments, units, errors, passed);
    passed = passed && errors.node <= tolerance && errors.weight <= tolerance && errors.moment <= tolerance;

    std::cout << std::scientific << std::setprecision(4) << "representation=" << representation << " max_normalized_node_error=" << errors.node
              << " max_weight_error=" << errors.weight << " max_relative_moment_error=" << errors.moment << " tolerance=" << tolerance
              << " verdict=" << (passed ? "PASS" : "FAIL") << '\n';
    return passed ? 0 : 1;
}
