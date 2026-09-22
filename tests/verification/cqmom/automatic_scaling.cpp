// CQMOM with automatic scaling enabled, on a problem whose coordinates span
// very different magnitudes.
//
// What is verified
//   That switching automatic scaling on does not change the answer: the
//   quadrature recovered from a rescaled problem, once mapped back, must match
//   the one recovered from the canonical problem, and both must reproduce the
//   moments they were given.
//
// How the reference is constructed
//   Chosen, not computed. A canonical rule of weighted points is written down
//   and a second copy is produced by rescaling each coordinate by a fixed
//   factor. The raw moments of both are evaluated from their own point sets, so
//   the two cases have independent references that are related by a known
//   transformation.
//
// Error metric and pass criteria
//   Maximum absolute deviation of positions, weights and moments from the
//   prescribed rules, at 2.0e-10, with 1.0e-10 where the comparison is
//   sensitive to the scaling step. These sit above the 1.0e-12 used by the
//   unscaled studies because scaling multiplies and divides the data by the
//   chosen factors, and each pass costs a little accuracy.
//
// Reproduce
//   ctest -R verification_cqmom_.*scaling
//
// Build modes
//   Both. Failures are reported through explicit comparisons and an exit code,
//   not assert().

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
constexpr double feature_sensitive_tolerance = 1.0e-10;

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
        terminals.push_back({std::move(coordinates), quadrature.weights()[node], quadrature.terminal_paths()[node].to_conditioning_path()});
    }
    std::sort(terminals.begin(), terminals.end(), [](const Terminal& lhs, const Terminal& rhs) { return lhs.path < rhs.path; });
    return terminals;
}

qmomkit::RawMomentTensor moments_from_terminals(const std::vector<Terminal>& terminals, const std::vector<std::size_t>& extents)
{
    qmomkit::RawMomentTensor moments{extents, {}};
    for (std::size_t a0 = 0; a0 < extents[0]; ++a0)
    {
        for (std::size_t a1 = 0; a1 < extents[1]; ++a1)
        {
            for (std::size_t a2 = 0; a2 < extents[2]; ++a2)
            {
                double value = 0.0;
                for (const Terminal& terminal : terminals)
                {
                    value += terminal.weight * std::pow(terminal.coordinates[0], static_cast<int>(a0)) *
                             std::pow(terminal.coordinates[1], static_cast<int>(a1)) * std::pow(terminal.coordinates[2], static_cast<int>(a2));
                }
                moments.values.push_back(value);
            }
        }
    }
    return moments;
}

void check_result(const std::string& route, const qmomkit::CqmomResult& result, const qmomkit::RawMomentTensor& configured,
                  const std::vector<double>& units, Errors& errors, bool& passed, const bool require_diagnostic_counts = true)
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
    if (require_diagnostic_counts && actual_counts != expected_counts)
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

int run_feature_sensitive_case()
{
    // The same prescribed ragged terminal rule as the covariance interaction
    // case is represented in deliberately disparate physical units. The
    // factors create mixed raw moments from approximately 1e-24 to 1e36 while
    // the normalized terminal rule and its moments remain O(1). The neighboring
    // probe at (1e6, 1e-6, 1e3) reconstructs without scaling; this stronger,
    // but still finite and non-singular, representation establishes that the
    // automatic-scaling path changes the outcome.
    const std::vector<double> units{1.0e8, 1.0e-8, 1.0e4};
    const std::vector<Terminal> expected = physical_reference(units);
    const qmomkit::RawMomentTensor moments = moments_from_terminals(expected, {4, 4, 4});
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {2, 0, 1});

    qmomkit::QmomOptions scaled_options;
    scaled_options.adaptive = true;
    scaled_options.automatic_scaling = true;
    const qmomkit::CqmomDenseCase scaled_case{hierarchy, moments, {}, scaled_options, true};

    bool passed = true;
    Errors scaled_errors;
    const qmomkit::CqmomResult scaled = qmomkit::compute_cqmom(scaled_case);
    check_result("feature-sensitive scaled", scaled, moments, units, scaled_errors, passed, false);
    if (scaled_errors.node > feature_sensitive_tolerance || scaled_errors.weight > feature_sensitive_tolerance ||
        scaled_errors.moment > feature_sensitive_tolerance)
    {
        std::cerr << "automatic scaling did not recover the independently prescribed physical rule.\n";
        passed = false;
    }

    qmomkit::CqmomDenseCase unscaled_case = scaled_case;
    unscaled_case.qmom_options.automatic_scaling = false;
    const qmomkit::CqmomResult unscaled = qmomkit::compute_cqmom(unscaled_case);
    const std::size_t unscaled_terminal_count = unscaled.ok() ? unscaled.quadrature.weights().size() : 0;
    if (!unscaled.ok() || unscaled_terminal_count >= expected.size())
    {
        std::cerr << "The deliberately difficult unscaled representation did not lose the prescribed five-terminal topology.\n";
        passed = false;
    }

    double unscaled_moment_error = 0.0;
    if (unscaled.ok())
    {
        const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(unscaled.quadrature, moments.extents);
        for (std::size_t index = 0; index < moments.values.size(); ++index)
        {
            unscaled_moment_error = std::max(unscaled_moment_error, relative_error(reconstructed.values[index], moments.values[index]));
        }
    }

    // Silent on success. These numbers exist to diagnose a failure, and CTest
    // surfaces stderr through --output-on-failure.
    if (!passed)
    {
        std::cerr << std::scientific << std::setprecision(4)
                  << "feature_sensitive units=(1e8,1e-8,1e4) scaled_terminals=" << (scaled.ok() ? scaled.quadrature.weights().size() : 0)
                  << " scaled_node_error=" << scaled_errors.node << " scaled_weight_error=" << scaled_errors.weight
                  << " scaled_moment_error=" << scaled_errors.moment << " scaled_tolerance=" << feature_sensitive_tolerance
                  << " unscaled_terminals=" << unscaled_terminal_count << " unscaled_moment_error=" << unscaled_moment_error << '\n';
    }
    return passed ? 0 : 1;
}

} // namespace

int run_configured_case(const std::filesystem::path& configuration_path, const std::string& representation)
{
    const std::vector<double> units = representation == "reference"  ? std::vector<double>{1.0, 1.0, 1.0}
                                      : representation == "rescaled" ? std::vector<double>{1.0e6, 1.0e-6, 1.0e3}
                                                                     : std::vector<double>{};
    if (units.empty())
    {
        std::cerr << "Unknown unit representation '" << representation << "'.\n";
        return 2;
    }

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

    if (!passed)
    {
        std::cerr << std::scientific << std::setprecision(4) << "representation=" << representation << " max_normalized_node_error=" << errors.node
                  << " max_weight_error=" << errors.weight << " max_relative_moment_error=" << errors.moment << " tolerance=" << tolerance << '\n';
    }
    return passed ? 0 : 1;
}

int main(const int argc, char** argv)
{
    if (argc == 2 && std::string(argv[1]) == "feature-sensitive")
    {
        return run_feature_sensitive_case();
    }
    if (argc == 2 && std::string(argv[1]) == "unit-invariance")
    {
        const std::filesystem::path case_directory = std::filesystem::path(__FILE__).parent_path();
        const int reference = run_configured_case(case_directory / "automatic_scaling_reference_3d.qmom", "reference");
        const int rescaled = run_configured_case(case_directory / "automatic_scaling_rescaled_3d.qmom", "rescaled");
        return reference == 0 && rescaled == 0 ? 0 : 1;
    }
    if (argc == 3)
    {
        return run_configured_case(argv[1], argv[2]);
    }
    std::cerr << "Usage: automatic_scaling <configuration.qmom> <reference|rescaled> | unit-invariance | feature-sensitive\n";
    return 2;
}
