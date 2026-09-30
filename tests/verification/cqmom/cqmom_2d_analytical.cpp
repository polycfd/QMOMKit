// Two-dimensional CQMOM reconstruction against discrete measures whose exact
// quadrature is known in closed form.
//
// What is verified
//   That compute_cqmom recovers the parent nodes and weights on the first
//   conditioning axis, the conditional nodes and weights on the second, the
//   terminal paths and global weights, and the adaptive branch node counts,
//   and that reconstructing raw moments from the recovered quadrature returns
//   the moments the case started from.
//
// How the reference is constructed
//   The reference is chosen, not computed. Each case below writes down the
//   answer directly: a handful of points in the (x, y) plane, each carrying a
//   weight, which is exactly what a CQMOM quadrature is. The parent list gives
//   the x positions and what share of the population sits at each; for every
//   parent, a second list gives the y positions reached from it and how that
//   parent's share splits between them. The moments in the matching .qmom file
//   are then the raw moments of that set of weighted points, worked out by
//   hand,
//
//       M[i,j] = sum_p w_p * x_p^i * sum_c v_(p,c) * y_(p,c)^j,
//
//   so the exact CQMOM answer is the set of points we started from. Nothing
//   in the reference is read back from QMOMKit, which is what makes this a
//   verification rather than a change detector. A two-node axis is exactly
//   determined by the four moments the hierarchy requires, so there is only
//   one quadrature that fits, and the comparison is against that one answer.
//
// Assumptions
//   Both axes use two requested nodes and extents {4, 4}; conditioning order
//   is {0, 1}, so x is the parent axis. The adaptive cases expect Wheeler to
//   return fewer nodes on a branch that reaches fewer y positions than the two
//   requested, which is why their branch counts differ from {2, 2}.
//
// Error metric
//   Maximum absolute error over abscissas, maximum absolute error over
//   weights, and maximum moment error relative to max(|expected|, 1). The
//   absolute moment error is reported for diagnosis but deliberately not
//   gated: moment magnitudes span 1 to ~5.6e3 across these cases, so an
//   absolute bound would be a different criterion per case. deterministic_affine
//   shows this directly — its absolute moment error of ~6e-12 exceeds the
//   tolerance while its relative error is ~1e-15.
//
// Tolerance
//   1.0e-12 on abscissas, weights and relative moment error. The decimal
//   moment literals are themselves only accurate to ~2.6e-16 relative, and the
//   observed errors are 3.1e-15 to 6.8e-15 on abscissas, 2.0e-15 to 2.9e-15 on
//   weights, and 3.8e-16 to 1.1e-15 on relative moments. The tolerance
//   therefore sits roughly two to three orders of magnitude above the
//   measured margin: loose enough not to pin a specific rounding path, tight
//   enough that a genuine inversion error cannot hide beneath it.
//
// Reproduce
//   ctest -R "verification_cqmom_(independent_product|deterministic_affine|branch_dependent)"
//   or run the executable directly with a .qmom path and the matching case
//   identifier.
//
// Build modes
//   Both. This study uses explicit comparisons and an exit code rather than
//   assert(), so unlike the unit tier it is not weakened by NDEBUG.

#include "qmomkit/configuration.hpp"
#include "qmomkit/moment_reconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace
{

// Gates abscissas, weights and the unit-floored relative moment error. The
// absolute moment error is reported but not gated; see the header.
constexpr double tolerance = 1.0e-12;

struct ErrorMetrics
{
    double node = 0.0;
    double weight = 0.0;
    double moment_absolute = 0.0;
    double moment_relative = 0.0;
};

struct ChildRule
{
    std::vector<double> y_nodes;
    std::vector<double> local_weights;
};

struct ExpectedCase
{
    std::string name;
    bool adaptive = false;
    std::vector<double> parent_nodes;
    std::vector<double> parent_weights;
    std::vector<ChildRule> children;
    std::vector<qmomkit::ConditioningPath> terminal_paths;
};

struct TerminalNode
{
    qmomkit::ConditioningPath path;
    double x = 0.0;
    double y = 0.0;
    double global_weight = 0.0;
};

struct VerificationResult
{
    bool passed = true;
    std::string message;
    ErrorMetrics errors;
    std::vector<std::size_t> recovered_branch_counts;
};

void append_failure(VerificationResult& result, const std::string& message)
{
    result.passed = false;
    if (!result.message.empty())
    {
        result.message += '\n';
    }
    result.message += message;
}

void update_error(double& maximum, const double actual, const double expected)
{
    maximum = std::max(maximum, std::abs(actual - expected));
}

void update_relative_error(double& maximum, const double actual, const double expected)
{
    const double absolute = std::abs(actual - expected);
    maximum = std::max(maximum, absolute / std::max(std::abs(expected), 1.0));
}

std::vector<std::pair<double, double>> sort_rule(std::vector<std::pair<double, double>> rule)
{
    std::sort(rule.begin(), rule.end(), [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });
    return rule;
}

// y is spread the same way whatever x is, so both parents lead to the same
// pair of y positions. Non-adaptive, so Wheeler must return the full {2, 2}
// branch counts and rectangular storage.
//
//   0.4*d(x-1)[0.25*d(y+1) + 0.75*d(y-2)] + 0.6*d(x-3)[0.25*d(y+1) + 0.75*d(y-2)]
ExpectedCase independent_product()
{
    return ExpectedCase{
        "independent_product",
        false,
        {1.0, 3.0},
        {0.4, 0.6},
        {{{-1.0, 2.0}, {0.25, 0.75}}, {{-1.0, 2.0}, {0.25, 0.75}}},
        {qmomkit::ConditioningPath({0, 0}), qmomkit::ConditioningPath({0, 1}), qmomkit::ConditioningPath({1, 0}), qmomkit::ConditioningPath({1, 1})}};
}

// Each parent leads to exactly one y position, lying on the straight line
// y = 2x + 1: parent x=1 gives y=3, parent x=3 gives y=7. Adaptive, so each
// branch must come back with one node instead of the two requested, giving
// branch counts {1, 1} and ragged storage. A branch that failed to reduce
// would have to invent a second y position that is not in the data.
//
//   0.4*d(x-1)d(y-3) + 0.6*d(x-3)d(y-7)
ExpectedCase deterministic_affine()
{
    return ExpectedCase{"deterministic_affine",
                        true,
                        {1.0, 3.0},
                        {0.4, 0.6},
                        {{{3.0}, {1.0}}, {{7.0}, {1.0}}},
                        {qmomkit::ConditioningPath({0, 0}), qmomkit::ConditioningPath({1, 0})}};
}

// The two parents lead to different numbers of y positions - one from x=1,
// two from x=3 - which is the shape rectangular storage cannot hold. Adaptive,
// so the expected branch counts are {1, 2}: one branch reduced, the other not,
// within a single reconstruction.
//
//   0.35*d(x-1)d(y+0.5) + 0.65*d(x-3)[0.3*d(y-1) + 0.7*d(y-4)]
ExpectedCase branch_dependent()
{
    return ExpectedCase{"branch_dependent",
                        true,
                        {1.0, 3.0},
                        {0.35, 0.65},
                        {{{-0.5}, {1.0}}, {{1.0, 4.0}, {0.3, 0.7}}},
                        {qmomkit::ConditioningPath({0, 0}), qmomkit::ConditioningPath({1, 0}), qmomkit::ConditioningPath({1, 1})}};
}

const ExpectedCase* find_case(const std::string& name)
{
    static const std::vector<ExpectedCase> cases{independent_product(), deterministic_affine(), branch_dependent()};
    for (const ExpectedCase& test_case : cases)
    {
        if (test_case.name == name)
        {
            return &test_case;
        }
    }
    return nullptr;
}

std::vector<TerminalNode> terminal_nodes(const qmomkit::CqmomQuadrature& quadrature)
{
    std::vector<TerminalNode> nodes;
    if (quadrature.has_rectangular_storage())
    {
        for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(quadrature.hierarchy()))
        {
            nodes.push_back({qmomkit::ConditioningPath(node_index), quadrature.abscissa(0, node_index), quadrature.abscissa(1, node_index),
                             quadrature.weight(node_index)});
        }
        return nodes;
    }

    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        nodes.push_back({quadrature.terminal_paths()[node].to_conditioning_path(), quadrature.abscissas()[0][node], quadrature.abscissas()[1][node],
                         quadrature.weights()[node]});
    }
    return nodes;
}

std::map<std::size_t, std::vector<TerminalNode>> group_by_parent(const std::vector<TerminalNode>& nodes)
{
    std::map<std::size_t, std::vector<TerminalNode>> grouped;
    for (const TerminalNode& node : nodes)
    {
        grouped[node.path.node_index_at_level(0)].push_back(node);
    }
    return grouped;
}

void verify_configuration_fields(const ExpectedCase& expected, const qmomkit::CqmomLeanCase& configured_case, VerificationResult& result)
{
    if (configured_case.moments.layout.hierarchy.node_counts() != std::vector<std::size_t>{2, 2})
    {
        append_failure(result, "Configured node_counts are not {2, 2}.");
    }
    if (configured_case.moments.layout.hierarchy.conditioning_order() != std::vector<std::size_t>{0, 1})
    {
        append_failure(result, "Configured conditioning_order is not {0, 1}.");
    }
    if (configured_case.moments.layout.hierarchy.moment_extents() != std::vector<std::size_t>{4, 4})
    {
        append_failure(result, "Configured moment extents are not {4, 4}.");
    }
    if (configured_case.adaptive_branch_counts != expected.adaptive)
    {
        append_failure(result, "Configured adaptive_branch_counts does not match the analytical case.");
    }
    if (configured_case.qmom_options.adaptive != expected.adaptive)
    {
        append_failure(result, "Configured qmom.adaptive does not match the analytical case.");
    }
}

void verify_diagnostics(const ExpectedCase& expected, const qmomkit::CqmomResult& cqmom, VerificationResult& result)
{
    const std::map<std::size_t, std::size_t> expected_child_counts = [&expected]
    {
        std::map<std::size_t, std::size_t> counts;
        for (std::size_t parent = 0; parent < expected.children.size(); ++parent)
        {
            counts[parent] = expected.children[parent].y_nodes.size();
        }
        return counts;
    }();

    bool saw_root = false;
    std::map<std::size_t, bool> saw_child_branch;
    for (const qmomkit::CqmomDiagnostics& diagnostics : cqmom.diagnostics)
    {
        if (diagnostics.qmom_status != qmomkit::QmomStatus::success ||
            diagnostics.conditional_moment_status != qmomkit::ConditionalMomentStatus::success)
        {
            append_failure(result, "CQMOM diagnostics contain a non-success status.");
        }
        if (diagnostics.requested_node_count != 2)
        {
            append_failure(result, "CQMOM diagnostics did not report requested_node_count = 2.");
        }
        if (diagnostics.adaptive_branch_counts_enabled != expected.adaptive)
        {
            append_failure(result, "CQMOM diagnostics adaptive flag does not match the analytical case.");
        }
        if (diagnostics.conditioning_level == 0 && diagnostics.parent_path.empty())
        {
            saw_root = true;
            if (diagnostics.actual_node_count != expected.parent_nodes.size())
            {
                append_failure(result, "Parent diagnostics did not recover two x nodes.");
            }
        }
        if (diagnostics.conditioning_level == 1 && diagnostics.parent_path.size() == 1)
        {
            const std::size_t parent = diagnostics.parent_path.node_index_at_level(0);
            const auto expected_count = expected_child_counts.find(parent);
            if (expected_count == expected_child_counts.end())
            {
                append_failure(result, "Diagnostics reported an unexpected parent branch.");
            }
            else
            {
                saw_child_branch[parent] = true;
                if (diagnostics.actual_node_count != expected_count->second)
                {
                    append_failure(result, "Diagnostics actual_node_count does not match expected child count.");
                }
                if (diagnostics.adaptive_node_count_reduced != (expected.adaptive && expected_count->second < diagnostics.requested_node_count))
                {
                    append_failure(result, "Diagnostics adaptive reduction flag is inconsistent.");
                }
            }
        }
    }
    if (!saw_root)
    {
        append_failure(result, "Did not find root parent diagnostics.");
    }
    for (const auto& expected_count : expected_child_counts)
    {
        if (!saw_child_branch[expected_count.first])
        {
            append_failure(result, "Did not find diagnostics for every parent branch.");
        }
    }
}

void verify_parent_and_children(const ExpectedCase& expected, const std::vector<TerminalNode>& nodes, VerificationResult& result)
{
    const std::map<std::size_t, std::vector<TerminalNode>> grouped = group_by_parent(nodes);
    result.recovered_branch_counts.clear();
    for (std::size_t parent = 0; parent < expected.parent_nodes.size(); ++parent)
    {
        const auto branch = grouped.find(parent);
        if (branch == grouped.end())
        {
            append_failure(result, "Missing expected parent branch.");
            continue;
        }

        const std::vector<TerminalNode>& branch_nodes = branch->second;
        result.recovered_branch_counts.push_back(branch_nodes.size());
        double parent_weight = 0.0;
        for (const TerminalNode& node : branch_nodes)
        {
            parent_weight += node.global_weight;
            update_error(result.errors.node, node.x, expected.parent_nodes[parent]);
        }
        update_error(result.errors.weight, parent_weight, expected.parent_weights[parent]);

        std::vector<std::pair<double, double>> actual_children;
        for (const TerminalNode& node : branch_nodes)
        {
            actual_children.emplace_back(node.y, node.global_weight / expected.parent_weights[parent]);
        }
        actual_children = sort_rule(std::move(actual_children));

        std::vector<std::pair<double, double>> expected_children;
        for (std::size_t child = 0; child < expected.children[parent].y_nodes.size(); ++child)
        {
            expected_children.emplace_back(expected.children[parent].y_nodes[child], expected.children[parent].local_weights[child]);
        }
        expected_children = sort_rule(std::move(expected_children));

        if (actual_children.size() != expected_children.size())
        {
            append_failure(result, "Recovered child count does not match expected conditional quadrature.");
            continue;
        }
        for (std::size_t child = 0; child < expected_children.size(); ++child)
        {
            update_error(result.errors.node, actual_children[child].first, expected_children[child].first);
            update_error(result.errors.weight, actual_children[child].second, expected_children[child].second);
        }
    }
}

void verify_terminal_nodes(const ExpectedCase& expected, std::vector<TerminalNode> nodes, VerificationResult& result)
{
    std::sort(nodes.begin(), nodes.end(), [](const TerminalNode& lhs, const TerminalNode& rhs) { return lhs.path < rhs.path; });
    if (nodes.size() != expected.terminal_paths.size())
    {
        append_failure(result, "Recovered terminal node count does not match expected terminal count.");
        return;
    }

    double total_weight = 0.0;
    for (std::size_t terminal = 0; terminal < nodes.size(); ++terminal)
    {
        const qmomkit::ConditioningPath& path = expected.terminal_paths[terminal];
        if (nodes[terminal].path != path)
        {
            append_failure(result, "Recovered terminal path does not match expected path.");
            continue;
        }
        const std::size_t parent = path.node_index_at_level(0);
        const std::size_t child = path.node_index_at_level(1);
        const double expected_weight = expected.parent_weights[parent] * expected.children[parent].local_weights[child];
        update_error(result.errors.weight, nodes[terminal].global_weight, expected_weight);
        update_error(result.errors.node, nodes[terminal].x, expected.parent_nodes[parent]);
        update_error(result.errors.node, nodes[terminal].y, expected.children[parent].y_nodes[child]);
        total_weight += nodes[terminal].global_weight;
        if (!(nodes[terminal].global_weight > 0.0))
        {
            append_failure(result, "Recovered non-positive terminal weight.");
        }
    }
    update_error(result.errors.weight, total_weight, 1.0);
}

void verify_reconstructed_moments(const qmomkit::CqmomQuadrature& quadrature, const qmomkit::CqmomMomentData& configured, VerificationResult& result)
{
    const qmomkit::CqmomMomentData reconstructed = qmomkit::reconstruct_cqmom_moments(quadrature, configured.layout);
    if (reconstructed.values.size() != configured.values.size())
    {
        append_failure(result, "Reconstructed CQMOM moment set has wrong size.");
        return;
    }
    for (std::size_t index = 0; index < configured.values.size(); ++index)
    {
        update_error(result.errors.moment_absolute, reconstructed.values[index], configured.values[index]);
        update_relative_error(result.errors.moment_relative, reconstructed.values[index], configured.values[index]);
    }
}

VerificationResult verify_case(const ExpectedCase& expected, const qmomkit::CqmomLeanCase& configured_case)
{
    VerificationResult result;
    verify_configuration_fields(expected, configured_case, result);

    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(configured_case);
    if (!cqmom.ok())
    {
        append_failure(result, "compute_cqmom failed: " + cqmom.message);
        return result;
    }
    if (cqmom.quadrature.has_rectangular_storage() == expected.adaptive)
    {
        append_failure(result, "Quadrature rectangular/adaptive storage mode does not match the analytical case.");
    }
    if (cqmom.quadrature.adaptive_branch_counts_enabled() != expected.adaptive)
    {
        append_failure(result, "Quadrature adaptive flag does not match the analytical case.");
    }

    const std::vector<TerminalNode> nodes = terminal_nodes(cqmom.quadrature);
    verify_diagnostics(expected, cqmom, result);
    verify_parent_and_children(expected, nodes, result);
    verify_terminal_nodes(expected, nodes, result);
    verify_reconstructed_moments(cqmom.quadrature, configured_case.moments, result);
    return result;
}

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::CqmomLeanCase>& configured)
{
    std::cerr << configured.error() << '\n';
}

} // namespace

int main(const int argc, char** argv)
{
    if (argc != 3)
    {
        std::cerr << "Usage: cqmom_2d_analytical_verification <configuration.qmom> <case_identifier>\n";
        return 2;
    }

    const std::filesystem::path configuration_path = argv[1];
    const std::string case_identifier = argv[2];
    const ExpectedCase* const expected = find_case(case_identifier);
    if (expected == nullptr)
    {
        std::cerr << "Unknown analytical CQMOM case identifier '" << case_identifier << "'.\n";
        return 2;
    }

    const auto configured = qmomkit::config::load_cqmom_lean(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        return 1;
    }

    const VerificationResult result = verify_case(*expected, configured.value());
    const bool passed =
        result.passed && result.errors.node <= tolerance && result.errors.weight <= tolerance && result.errors.moment_relative <= tolerance;
    // Silent on success. These numbers exist to diagnose a failure, and CTest
    // surfaces stderr through --output-on-failure.
    if (!(passed))
    {
        std::cerr << "QMOMKit two-dimensional CQMOM analytical verification\n";
        std::cerr << "Configuration: " << configuration_path << '\n';
        std::cerr << "Tolerance: " << std::scientific << tolerance << " (absolute, with unit-floored relative moment error)\n\n";
        std::cerr << std::left << std::setw(28) << "case" << std::right << std::setw(14) << "node abs" << std::setw(14) << "weight abs"
                  << std::setw(14) << "moment abs" << std::setw(14) << "moment rel" << std::setw(16) << "branch counts" << std::setw(9) << "status"
                  << '\n';
        std::cerr << std::string(109, '-') << '\n';
        std::cerr << std::left << std::setw(28) << expected->name << std::right << std::scientific << std::setprecision(4) << std::setw(14)
                  << result.errors.node << std::setw(14) << result.errors.weight << std::setw(14) << result.errors.moment_absolute << std::setw(14)
                  << result.errors.moment_relative;
        std::cerr << std::setw(8) << "[";
        for (std::size_t index = 0; index < result.recovered_branch_counts.size(); ++index)
        {
            if (index != 0)
            {
                std::cerr << ',';
            }
            std::cerr << result.recovered_branch_counts[index];
        }
        std::cerr << "]" << std::setw(9) << (passed ? "PASS" : "FAIL") << '\n';

        if (!result.passed)
        {
            std::cerr << result.message << '\n';
        }
    }

    return passed ? 0 : 1;
}
