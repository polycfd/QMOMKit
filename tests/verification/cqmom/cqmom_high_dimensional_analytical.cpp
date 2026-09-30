// CQMOM reconstruction in three to six dimensions, against sets of weighted
// points chosen in advance so the right answer is known exactly.
//
// What is verified
//   That compute_cqmom recovers every terminal point, its local and global
//   weight, and the branch node counts, in the original phase-space dimension
//   order rather than the conditioning order; and that raw moments rebuilt
//   from the recovered quadrature return the moments the case started from.
//
// How the reference is constructed
//   Chosen, not computed. Two families:
//
//   independent_product_Nd - each coordinate gets its own pair of positions
//     and weights, independent of the others, so the full set of points is
//     every combination of them and the weight of a point is the product of
//     the per-coordinate weights.
//
//   affine_chain_5d and mixed_affine_independent_4d - the first coordinate
//     gets two positions; each later coordinate is fixed by a straight-line
//     relation to the previous one (x1 = 2*x0 + 1, and so on). Each branch
//     therefore reaches exactly one position on the next coordinate, and the
//     adaptive path must return one node there instead of two.
//
//   The moments in the matching .qmom file are the raw moments of those
//   weighted points, so nothing in the reference comes from QMOMKit. Checked:
//   all 1024 moments of affine_chain_5d and all 64 of independent_product_3d
//   agree with the hand-evaluated values to 3.9e-16 relative or better.
//
// Error metric
//   Maximum absolute error over coordinates, over local weights, and over
//   global weights, plus maximum moment error relative to max(|expected|, 1).
//   The absolute moment error is printed but not gated, because moment
//   magnitudes vary by orders of magnitude across six dimensions.
//
// Tolerance
//   1.0e-11 on coordinates, local and global weights, and relative moment
//   error. This is set by the six-dimensional case, not chosen loosely: the
//   error grows steeply with dimension, from 9.4e-16 on coordinates at three
//   dimensions to 5.1e-13 at six, because each extra conditioning level feeds
//   one inversion into the next. The tolerance leaves roughly a factor of
//   twenty over the worst case; the three- and four-dimensional cases sit
//   three to four orders below it.
//
// Reproduce
//   ctest -R "verification_cqmom_(affine_chain_5d|independent_product_[3-6]d|mixed_affine_independent_4d)"
//
// Build modes
//   Both. Explicit comparisons and an exit code, not assert(), so NDEBUG does
//   not weaken this study.

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
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace
{

constexpr double tolerance = 1.0e-11;

struct ErrorMetrics
{
    double node = 0.0;
    double local_weight = 0.0;
    double global_weight = 0.0;
    double moment_absolute = 0.0;
    double moment_relative = 0.0;
};

struct TerminalPoint
{
    std::vector<double> coordinates;
    double global_weight = 0.0;
    qmomkit::ConditioningPath path;
    std::vector<double> local_weights;
};

struct AffineRelation
{
    std::size_t target_dimension = 0;
    std::vector<std::pair<std::size_t, double>> terms;
    double intercept = 0.0;
};

struct ExpectedCase
{
    std::string name;
    std::vector<std::size_t> node_counts;
    std::vector<std::size_t> conditioning_order;
    std::vector<std::size_t> extents;
    bool adaptive_wheeler = false;
    bool adaptive_branch_counts = false;
    std::vector<TerminalPoint> terminals;
    std::vector<AffineRelation> affine_relations;
};

struct RecoveredTerminal
{
    std::vector<double> coordinates;
    double global_weight = 0.0;
    qmomkit::ConditioningPath path;
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
    maximum = std::max(maximum, std::abs(actual - expected) / std::max(1.0, std::abs(expected)));
}

double integer_power(const double base, const std::size_t exponent)
{
    double value = 1.0;
    for (std::size_t order = 0; order < exponent; ++order)
    {
        value *= base;
    }
    return value;
}

std::string path_key(const qmomkit::ConditioningPath& path)
{
    std::ostringstream output;
    for (std::size_t index = 0; index < path.size(); ++index)
    {
        if (index != 0)
        {
            output << ',';
        }
        output << path.node_index_at_level(index);
    }
    return output.str();
}

bool lexicographic_coordinates_less(const std::vector<double>& lhs, const std::vector<double>& rhs)
{
    for (std::size_t dimension = 0; dimension < lhs.size(); ++dimension)
    {
        if (lhs[dimension] < rhs[dimension])
        {
            return true;
        }
        if (rhs[dimension] < lhs[dimension])
        {
            return false;
        }
    }
    return false;
}

qmomkit::ConditioningPath rectangular_path_from_node_index(const qmomkit::MultiIndex& node_index, const std::vector<std::size_t>& conditioning_order)
{
    std::vector<std::size_t> path(conditioning_order.size(), 0);
    for (std::size_t level = 0; level < conditioning_order.size(); ++level)
    {
        path[level] = node_index[conditioning_order[level]];
    }
    return qmomkit::ConditioningPath(path);
}

std::vector<RecoveredTerminal> recovered_terminals(const qmomkit::CqmomQuadrature& quadrature)
{
    std::vector<RecoveredTerminal> terminals;
    if (quadrature.has_rectangular_storage())
    {
        for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(quadrature.hierarchy()))
        {
            std::vector<double> coordinates(quadrature.number_of_dimensions(), 0.0);
            for (std::size_t dimension = 0; dimension < coordinates.size(); ++dimension)
            {
                coordinates[dimension] = quadrature.abscissa(dimension, node_index);
            }
            terminals.push_back({std::move(coordinates), quadrature.weight(node_index),
                                 rectangular_path_from_node_index(node_index, quadrature.hierarchy().conditioning_order())});
        }
        return terminals;
    }

    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        std::vector<double> coordinates(quadrature.number_of_dimensions(), 0.0);
        for (std::size_t dimension = 0; dimension < coordinates.size(); ++dimension)
        {
            coordinates[dimension] = quadrature.abscissas()[dimension][node];
        }
        terminals.push_back({std::move(coordinates), quadrature.weights()[node], quadrature.terminal_paths()[node].to_conditioning_path()});
    }
    return terminals;
}

std::vector<TerminalPoint> sort_expected_by_path(std::vector<TerminalPoint> terminals)
{
    std::sort(terminals.begin(), terminals.end(), [](const TerminalPoint& lhs, const TerminalPoint& rhs) { return lhs.path < rhs.path; });
    return terminals;
}

std::vector<RecoveredTerminal> sort_recovered_by_path(std::vector<RecoveredTerminal> terminals)
{
    std::sort(terminals.begin(), terminals.end(), [](const RecoveredTerminal& lhs, const RecoveredTerminal& rhs) { return lhs.path < rhs.path; });
    return terminals;
}

std::vector<std::pair<double, double>> sort_rule(std::vector<std::pair<double, double>> rule)
{
    std::sort(rule.begin(), rule.end(), [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });
    return rule;
}

std::vector<TerminalPoint> sort_expected_by_coordinates(std::vector<TerminalPoint> terminals)
{
    std::sort(terminals.begin(), terminals.end(),
              [](const TerminalPoint& lhs, const TerminalPoint& rhs) { return lexicographic_coordinates_less(lhs.coordinates, rhs.coordinates); });
    return terminals;
}

std::vector<RecoveredTerminal> sort_recovered_by_coordinates(std::vector<RecoveredTerminal> terminals)
{
    std::sort(terminals.begin(), terminals.end(), [](const RecoveredTerminal& lhs, const RecoveredTerminal& rhs)
              { return lexicographic_coordinates_less(lhs.coordinates, rhs.coordinates); });
    return terminals;
}

std::vector<std::size_t> prefix_vector(const qmomkit::ConditioningPath& path, const std::size_t size)
{
    std::vector<std::size_t> prefix(size, 0);
    for (std::size_t index = 0; index < size; ++index)
    {
        prefix[index] = path.node_index_at_level(index);
    }
    return prefix;
}

qmomkit::ConditioningPath prefix_path(const qmomkit::ConditioningPath& path, const std::size_t size)
{
    return qmomkit::ConditioningPath(prefix_vector(path, size));
}

std::vector<TerminalPoint> affine_chain_terminals()
{
    // Manufactured adaptive case: every later coordinate is an affine function
    // of the previous one, so each conditional branch after the root has only
    // one mathematically distinct child. Adaptive branch counts should collapse
    // those deterministic branches while preserving terminal weights and raw
    // moments in the original phase-space dimension order.
    const std::vector<double> x0_nodes{-0.5, 1.5};
    const std::vector<double> x0_weights{0.4, 0.6};
    std::vector<TerminalPoint> terminals;
    for (std::size_t node = 0; node < x0_nodes.size(); ++node)
    {
        std::vector<double> x(5, 0.0);
        x[0] = x0_nodes[node];
        x[1] = 2.0 * x[0] + 1.0;
        x[2] = -0.5 * x[1] + 0.75;
        x[3] = 1.25 * x[2] - 0.2;
        x[4] = -0.8 * x[3] + 0.5;
        terminals.push_back({std::move(x), x0_weights[node], qmomkit::ConditioningPath({node, 0, 0, 0, 0}), {x0_weights[node], 1.0, 1.0, 1.0, 1.0}});
    }
    return terminals;
}

struct Rule1D
{
    std::vector<double> nodes;
    std::vector<double> weights;
};

std::vector<Rule1D> independent_rules(const std::size_t dimensions)
{
    const std::vector<Rule1D> rules{{{-1.0, 1.5}, {0.35, 0.65}}, {{0.2, 1.2}, {0.6, 0.4}},    {{-0.5, 0.75}, {0.45, 0.55}},
                                    {{2.0, 3.5}, {0.7, 0.3}},    {{-2.0, 0.5}, {0.25, 0.75}}, {{1.0, 2.25}, {0.8, 0.2}}};
    return {rules.begin(), rules.begin() + static_cast<std::ptrdiff_t>(dimensions)};
}

bool is_independent_product_case(const ExpectedCase& expected)
{
    return expected.name.rfind("independent_product_", 0) == 0;
}

std::vector<TerminalPoint> independent_product_terminals(const std::size_t dimensions, const std::vector<std::size_t>& conditioning_order)
{
    // Independent-product cases use tensor products of simple 1D two-node
    // rules. They are deliberately conditioned in nontrivial dimension orders
    // so the verification separates phase-space ordering from recursive CQMOM
    // conditioning order.
    const std::vector<Rule1D> rules = independent_rules(dimensions);
    std::vector<TerminalPoint> terminals;
    const std::size_t count = std::size_t{1} << dimensions;
    for (std::size_t flat = 0; flat < count; ++flat)
    {
        std::vector<std::size_t> node_index(dimensions, 0);
        std::vector<double> coordinates(dimensions, 0.0);
        double weight = 1.0;
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
        {
            node_index[dimension] = (flat >> dimension) & std::size_t{1};
            coordinates[dimension] = rules[dimension].nodes[node_index[dimension]];
            weight *= rules[dimension].weights[node_index[dimension]];
        }

        std::vector<std::size_t> path(dimensions, 0);
        std::vector<double> local_weights(dimensions, 0.0);
        for (std::size_t level = 0; level < dimensions; ++level)
        {
            const std::size_t dimension = conditioning_order[level];
            path[level] = node_index[dimension];
            local_weights[level] = rules[dimension].weights[node_index[dimension]];
        }
        terminals.push_back({std::move(coordinates), weight, qmomkit::ConditioningPath(path), std::move(local_weights)});
    }
    return terminals;
}

std::vector<TerminalPoint> mixed_terminals()
{
    // Mixed manufactured case: dimensions 0 and 2 are independent two-node
    // sources, while dimensions 1 and 3 are deterministic affine descendants.
    // This exercises adaptive branch reduction in only part of a 4D hierarchy.
    const Rule1D xi0{{-0.5, 1.0}, {0.45, 0.55}};
    const Rule1D xi2{{-1.0, 2.0}, {0.3, 0.7}};
    std::vector<TerminalPoint> terminals;
    for (std::size_t i = 0; i < xi0.nodes.size(); ++i)
    {
        for (std::size_t k = 0; k < xi2.nodes.size(); ++k)
        {
            std::vector<double> x(4, 0.0);
            x[0] = xi0.nodes[i];
            x[1] = 2.0 * x[0] + 1.0;
            x[2] = xi2.nodes[k];
            x[3] = -x[1] + 0.5 * x[2];
            terminals.push_back(
                {std::move(x), xi0.weights[i] * xi2.weights[k], qmomkit::ConditioningPath({i, k, 0, 0}), {xi0.weights[i], xi2.weights[k], 1.0, 1.0}});
        }
    }
    return terminals;
}

ExpectedCase affine_chain_5d()
{
    return {"affine_chain_5d",
            {2, 2, 2, 2, 2},
            {0, 1, 2, 3, 4},
            {4, 4, 4, 4, 4},
            true,
            true,
            affine_chain_terminals(),
            {{1, {{0, 2.0}}, 1.0}, {2, {{1, -0.5}}, 0.75}, {3, {{2, 1.25}}, -0.2}, {4, {{3, -0.8}}, 0.5}}};
}

ExpectedCase independent_product_case(const std::size_t dimensions)
{
    const std::map<std::size_t, std::vector<std::size_t>> orders{{3, {1, 0, 2}}, {4, {2, 0, 3, 1}}, {5, {4, 1, 3, 0, 2}}, {6, {5, 2, 0, 4, 1, 3}}};
    const std::vector<std::size_t> order = orders.at(dimensions);
    return {"independent_product_" + std::to_string(dimensions) + "d",
            std::vector<std::size_t>(dimensions, 2),
            order,
            std::vector<std::size_t>(dimensions, 4),
            false,
            false,
            independent_product_terminals(dimensions, order),
            {}};
}

ExpectedCase mixed_affine_independent_4d()
{
    return {"mixed_affine_independent_4d",
            {2, 2, 2, 2},
            {0, 2, 1, 3},
            {4, 4, 4, 4},
            true,
            true,
            mixed_terminals(),
            {{1, {{0, 2.0}}, 1.0}, {3, {{1, -1.0}, {2, 0.5}}, 0.0}}};
}

const ExpectedCase* find_case(const std::string& name)
{
    static const std::vector<ExpectedCase> cases{affine_chain_5d(),           independent_product_case(3), independent_product_case(4),
                                                 independent_product_case(5), independent_product_case(6), mixed_affine_independent_4d()};
    for (const ExpectedCase& test_case : cases)
    {
        if (test_case.name == name)
        {
            return &test_case;
        }
    }
    return nullptr;
}

qmomkit::RawMomentTensor expected_raw_moments(const ExpectedCase& expected)
{
    // Build the dense reference tensor directly from the manufactured terminal
    // quadrature. This keeps the expected raw moments independent of the CQMOM
    // reconstruction routine under test.
    qmomkit::RawMomentTensor moments{expected.extents, std::vector<double>(qmomkit::product_extents(expected.extents), 0.0)};
    for (std::size_t flat = 0; flat < moments.values.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, expected.extents);
        for (const TerminalPoint& terminal : expected.terminals)
        {
            double contribution = terminal.global_weight;
            for (std::size_t dimension = 0; dimension < expected.extents.size(); ++dimension)
            {
                contribution *= integer_power(terminal.coordinates[dimension], index[dimension]);
            }
            moments.values[flat] += contribution;
        }
    }
    return moments;
}

void verify_configuration(const ExpectedCase& expected, const qmomkit::CqmomLeanCase& configured_case, VerificationResult& result)
{
    // The .qmom file is part of the reference definition. Before running the
    // algorithm, verify that it encodes the intended hierarchy, adaptive policy,
    // and analytical moment tensor.
    if (configured_case.moments.layout.hierarchy.node_counts() != expected.node_counts)
    {
        append_failure(result, "Configured node_counts do not match expected hierarchy.");
    }
    if (configured_case.moments.layout.hierarchy.conditioning_order() != expected.conditioning_order)
    {
        append_failure(result, "Configured conditioning_order does not match expected case.");
    }
    if (configured_case.moments.layout.hierarchy.moment_extents() != expected.extents)
    {
        append_failure(result, "Configured extents do not match expected case.");
    }
    if (configured_case.qmom_options.adaptive != expected.adaptive_wheeler)
    {
        append_failure(result, "Configured qmom.adaptive does not match expected case.");
    }
    if (configured_case.adaptive_branch_counts != expected.adaptive_branch_counts)
    {
        append_failure(result, "Configured adaptive_branch_counts does not match expected case.");
    }

    const qmomkit::CqmomMomentData analytical = qmomkit::extract_cqmom_moments(expected_raw_moments(expected), configured_case.moments.layout);
    if (analytical.values.size() != configured_case.moments.values.size())
    {
        append_failure(result, "Configured moment count does not match the analytical reference.");
        return;
    }
    for (std::size_t index = 0; index < analytical.values.size(); ++index)
    {
        update_error(result.errors.moment_absolute, configured_case.moments.values[index], analytical.values[index]);
        update_relative_error(result.errors.moment_relative, configured_case.moments.values[index], analytical.values[index]);
    }
}

std::map<std::string, std::size_t> expected_branch_counts(const ExpectedCase& expected)
{
    std::map<std::string, std::size_t> counts;
    for (std::size_t level = 0; level < expected.conditioning_order.size(); ++level)
    {
        std::map<std::string, std::vector<std::size_t>> children;
        for (const TerminalPoint& terminal : expected.terminals)
        {
            const qmomkit::ConditioningPath parent = prefix_path(terminal.path, level);
            const std::size_t child = terminal.path.node_index_at_level(level);
            std::vector<std::size_t>& values = children[path_key(parent)];
            if (std::find(values.begin(), values.end(), child) == values.end())
            {
                values.push_back(child);
            }
        }
        for (const auto& item : children)
        {
            counts[std::to_string(level) + ":" + item.first] = item.second.size();
        }
    }
    return counts;
}

void verify_diagnostics(const ExpectedCase& expected, const qmomkit::CqmomResult& cqmom, VerificationResult& result)
{
    // Diagnostics are checked as executable documentation of the recursive
    // branch structure: each public branch record should follow the configured
    // conditioning order and, when adaptive branch counts are enabled, report
    // the mathematically reduced child count.
    std::map<std::string, std::size_t> expected_counts = expected_branch_counts(expected);
    result.recovered_branch_counts.clear();
    for (const qmomkit::CqmomDiagnostics& diagnostics : cqmom.diagnostics)
    {
        if (diagnostics.qmom_status != qmomkit::QmomStatus::success ||
            diagnostics.conditional_moment_status != qmomkit::ConditionalMomentStatus::success)
        {
            append_failure(result, "CQMOM diagnostics contain a non-success status.");
        }
        if (diagnostics.phase_space_dimension != expected.conditioning_order[diagnostics.conditioning_level])
        {
            append_failure(result, "Diagnostic phase-space dimension does not follow conditioning_order.");
        }
        if (diagnostics.requested_node_count != expected.node_counts[diagnostics.phase_space_dimension])
        {
            append_failure(result, "Diagnostic requested_node_count does not match configured hierarchy.");
        }
        if (diagnostics.adaptive_branch_counts_enabled != expected.adaptive_branch_counts)
        {
            append_failure(result, "Diagnostic adaptive_branch_counts flag does not match expected case.");
        }
        const std::string key = std::to_string(diagnostics.conditioning_level) + ":" + path_key(diagnostics.parent_path);
        const auto found = expected_counts.find(key);
        if (found == expected_counts.end())
        {
            append_failure(result, "Diagnostics reported an unexpected branch path.");
            continue;
        }
        result.recovered_branch_counts.push_back(diagnostics.actual_node_count);
        if (diagnostics.actual_node_count != found->second)
        {
            append_failure(result, "Diagnostic actual_node_count does not match expected branch count.");
        }
        const bool expected_reduced = expected.adaptive_branch_counts && found->second < diagnostics.requested_node_count;
        if (diagnostics.adaptive_node_count_reduced != expected_reduced)
        {
            append_failure(result, "Diagnostic adaptive reduction flag is inconsistent.");
        }
        expected_counts.erase(found);
    }
    if (!expected_counts.empty())
    {
        append_failure(result, "Diagnostics did not cover every expected branch.");
    }
    std::sort(result.recovered_branch_counts.begin(), result.recovered_branch_counts.end());
}

void verify_terminals(const ExpectedCase& expected, std::vector<RecoveredTerminal> recovered, VerificationResult& result)
{
    // Terminal coordinates and global weights are the primary CQMOM result.
    // Adaptive cases compare by conditioning path because reduced branches have
    // meaningful public paths; fixed rectangular cases match by coordinates so
    // equivalent tensor-product ordering differences do not hide correctness.
    std::vector<TerminalPoint> analytical = expected.adaptive_branch_counts ? sort_expected_by_path(expected.terminals) : expected.terminals;
    recovered = expected.adaptive_branch_counts ? sort_recovered_by_path(std::move(recovered)) : std::move(recovered);
    if (recovered.size() != analytical.size())
    {
        append_failure(result, "Recovered terminal count does not match expected terminal count.");
        return;
    }

    double recovered_total_weight = 0.0;
    double expected_total_weight = 0.0;
    std::vector<bool> matched(recovered.size(), false);
    for (std::size_t terminal = 0; terminal < analytical.size(); ++terminal)
    {
        std::size_t recovered_index = terminal;
        if (!expected.adaptive_branch_counts)
        {
            double best_error = std::numeric_limits<double>::infinity();
            for (std::size_t candidate = 0; candidate < recovered.size(); ++candidate)
            {
                if (matched[candidate])
                {
                    continue;
                }
                double candidate_error = 0.0;
                for (std::size_t dimension = 0; dimension < analytical[terminal].coordinates.size(); ++dimension)
                {
                    candidate_error = std::max(candidate_error,
                                               std::abs(recovered[candidate].coordinates[dimension] - analytical[terminal].coordinates[dimension]));
                }
                if (candidate_error < best_error)
                {
                    best_error = candidate_error;
                    recovered_index = candidate;
                }
            }
            matched[recovered_index] = true;
        }

        if (expected.adaptive_branch_counts && recovered[recovered_index].path != analytical[terminal].path)
        {
            append_failure(result, "Recovered terminal path does not match analytical path.");
        }
        if (recovered[recovered_index].coordinates.size() != analytical[terminal].coordinates.size())
        {
            append_failure(result, "Recovered coordinate dimension does not match analytical dimension.");
            continue;
        }
        for (std::size_t dimension = 0; dimension < analytical[terminal].coordinates.size(); ++dimension)
        {
            update_error(result.errors.node, recovered[recovered_index].coordinates[dimension], analytical[terminal].coordinates[dimension]);
        }
        update_error(result.errors.global_weight, recovered[recovered_index].global_weight, analytical[terminal].global_weight);
        if (!(recovered[recovered_index].global_weight > 0.0))
        {
            append_failure(result, "Recovered non-positive terminal weight.");
        }
        for (const AffineRelation& relation : expected.affine_relations)
        {
            double expected_value = relation.intercept;
            for (const auto& term : relation.terms)
            {
                expected_value += term.second * recovered[recovered_index].coordinates[term.first];
            }
            update_error(result.errors.node, recovered[recovered_index].coordinates[relation.target_dimension], expected_value);
        }
        recovered_total_weight += recovered[recovered_index].global_weight;
        expected_total_weight += analytical[terminal].global_weight;
    }
    update_error(result.errors.global_weight, recovered_total_weight, expected_total_weight);
}

std::map<std::string, double> branch_weights_from_terminals(const std::vector<RecoveredTerminal>& recovered)
{
    std::map<std::string, double> branch_weight;
    branch_weight[""] = 1.0;
    for (const RecoveredTerminal& terminal : recovered)
    {
        for (std::size_t prefix_size = 1; prefix_size <= terminal.path.size(); ++prefix_size)
        {
            branch_weight[path_key(prefix_path(terminal.path, prefix_size))] += terminal.global_weight;
        }
    }
    return branch_weight;
}

void verify_local_weights(const ExpectedCase& expected, const std::vector<RecoveredTerminal>& recovered, VerificationResult& result)
{
    // Global terminal weights alone can mask an incorrect conditional split.
    // Reconstructing local branch weights from terminal masses verifies the
    // CQMOM conditional probabilities exposed by the recursive hierarchy.
    const std::map<std::string, double> branch_weight = branch_weights_from_terminals(recovered);

    std::map<std::string, double> expected_local_weight;
    for (const TerminalPoint& terminal : expected.terminals)
    {
        for (std::size_t level = 0; level < expected.conditioning_order.size(); ++level)
        {
            expected_local_weight.emplace(path_key(prefix_path(terminal.path, level + 1)), terminal.local_weights[level]);
        }
    }

    for (std::size_t level = 0; level < expected.conditioning_order.size(); ++level)
    {
        std::vector<std::string> checked_children;
        for (const RecoveredTerminal& terminal : recovered)
        {
            const std::string parent_key = path_key(prefix_path(terminal.path, level));
            const std::string child_key = path_key(prefix_path(terminal.path, level + 1));
            if (std::find(checked_children.begin(), checked_children.end(), child_key) != checked_children.end())
            {
                continue;
            }
            checked_children.push_back(child_key);
            const double parent_weight = branch_weight.at(parent_key);
            if (parent_weight == 0.0)
            {
                append_failure(result, "Recovered parent branch has zero accumulated weight.");
                continue;
            }
            const auto expected_local = expected_local_weight.find(child_key);
            if (expected_local == expected_local_weight.end())
            {
                append_failure(result, "Recovered unexpected child branch while checking local weights.");
                continue;
            }
            update_error(result.errors.local_weight, branch_weight.at(child_key) / parent_weight, expected_local->second);
        }
    }
}

void verify_independent_conditional_rules(const ExpectedCase& expected, const std::vector<RecoveredTerminal>& recovered, VerificationResult& result)
{
    // For independent-product references, every conditional branch should see
    // the same 1D marginal rule for the conditioned dimension. This protects
    // branch-local Vandermonde behavior separately from terminal reconstruction.
    if (!is_independent_product_case(expected))
    {
        return;
    }

    const std::vector<Rule1D> rules = independent_rules(expected.conditioning_order.size());
    const std::map<std::string, double> branch_weight = branch_weights_from_terminals(recovered);

    for (std::size_t level = 0; level < expected.conditioning_order.size(); ++level)
    {
        const std::size_t dimension = expected.conditioning_order[level];
        std::map<std::string, std::vector<std::pair<double, double>>> branch_rules;
        for (const RecoveredTerminal& terminal : recovered)
        {
            const std::string parent_key = path_key(prefix_path(terminal.path, level));
            const std::string child_key = path_key(prefix_path(terminal.path, level + 1));
            std::vector<std::pair<double, double>>& rule = branch_rules[parent_key];
            const bool already_recorded = std::any_of(rule.begin(), rule.end(), [&terminal, dimension](const auto& item)
                                                      { return std::abs(item.first - terminal.coordinates[dimension]) <= tolerance; });
            if (already_recorded)
            {
                continue;
            }
            const double parent_weight = branch_weight.at(parent_key);
            if (parent_weight == 0.0)
            {
                append_failure(result, "Independent-product branch has zero parent weight.");
                continue;
            }
            rule.emplace_back(terminal.coordinates[dimension], branch_weight.at(child_key) / parent_weight);
        }

        std::vector<std::pair<double, double>> expected_rule;
        for (std::size_t node = 0; node < rules[dimension].nodes.size(); ++node)
        {
            expected_rule.emplace_back(rules[dimension].nodes[node], rules[dimension].weights[node]);
        }
        expected_rule = sort_rule(std::move(expected_rule));

        for (auto& item : branch_rules)
        {
            item.second = sort_rule(std::move(item.second));
            if (item.second.size() != expected_rule.size())
            {
                append_failure(result, "Independent-product conditional rule has the wrong node count.");
                continue;
            }
            for (std::size_t node = 0; node < expected_rule.size(); ++node)
            {
                update_error(result.errors.node, item.second[node].first, expected_rule[node].first);
                update_error(result.errors.local_weight, item.second[node].second, expected_rule[node].second);
            }
        }
    }
}

void verify_reconstructed_moments(const qmomkit::CqmomQuadrature& quadrature, const qmomkit::RawMomentTensor& configured, VerificationResult& result)
{
    // The final quadrature must reproduce the full configured raw moment tensor
    // across all high-dimensional multi-indices, not only the minimal recursive
    // subset used internally by CQMOM.
    const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(quadrature, configured.extents);
    if (reconstructed.values.size() != configured.values.size())
    {
        append_failure(result, "Reconstructed moment tensor has wrong size.");
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
    // Verify one high-dimensional CQMOM workflow from configured dense moments
    // through public compute_cqmom output. The reference is the manufactured
    // terminal quadrature plus moments assembled above, so failures point to
    // hierarchy, branch-count, terminal, or reconstruction behavior.
    VerificationResult result;
    verify_configuration(expected, configured_case, result);

    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(configured_case);
    if (!cqmom.ok())
    {
        append_failure(result, "compute_cqmom failed: " + cqmom.message);
        return result;
    }
    if (cqmom.quadrature.has_rectangular_storage() == expected.adaptive_branch_counts)
    {
        append_failure(result, "Quadrature storage mode does not match fixed/adaptive expectation.");
    }
    if (cqmom.quadrature.adaptive_branch_counts_enabled() != expected.adaptive_branch_counts)
    {
        append_failure(result, "Quadrature adaptive_branch_counts flag does not match expected case.");
    }

    const std::vector<RecoveredTerminal> terminals = recovered_terminals(cqmom.quadrature);
    verify_diagnostics(expected, cqmom, result);
    verify_terminals(expected, terminals, result);
    verify_local_weights(expected, terminals, result);
    verify_independent_conditional_rules(expected, terminals, result);
    verify_reconstructed_moments(cqmom.quadrature, expected_raw_moments(expected), result);
    return result;
}

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::CqmomLeanCase>& configured)
{
    std::cerr << configured.error() << '\n';
}

// Called only from the failure report below, so it writes to stderr.
void print_counts(const std::vector<std::size_t>& counts)
{
    std::cerr << '[';
    for (std::size_t index = 0; index < counts.size(); ++index)
    {
        if (index != 0)
        {
            std::cerr << ',';
        }
        std::cerr << counts[index];
    }
    std::cerr << ']';
}

} // namespace

int main(const int argc, char** argv)
{
    // Command-line driver shared by 3D-6D independent-product cases and
    // adaptive affine cases. Each case identifier selects both the mathematical
    // reference and the expected diagnostic/adaptive behavior.
    if (argc != 3)
    {
        std::cerr << "Usage: cqmom_high_dimensional_analytical <configuration.qmom> <case_identifier>\n";
        return 2;
    }

    const std::filesystem::path configuration_path = argv[1];
    const std::string case_identifier = argv[2];
    const ExpectedCase* const expected = find_case(case_identifier);
    if (expected == nullptr)
    {
        std::cerr << "Unknown high-dimensional CQMOM case identifier '" << case_identifier << "'.\n";
        return 2;
    }

    const auto configured = qmomkit::config::load_cqmom_lean(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        return 1;
    }

    const VerificationResult result = verify_case(*expected, configured.value());
    const bool passed = result.passed && result.errors.node <= tolerance && result.errors.local_weight <= tolerance &&
                        result.errors.global_weight <= tolerance && result.errors.moment_relative <= tolerance;

    // Silent on success. These numbers exist to diagnose a failure, and CTest
    // surfaces stderr through --output-on-failure.
    if (!(passed))
    {
        std::cerr << "QMOMKit high-dimensional CQMOM analytical verification\n";
        std::cerr << "Configuration: " << configuration_path << '\n';
        std::cerr << "Tolerance: " << std::scientific << tolerance << " (absolute nodes/weights, unit-floored relative moments)\n\n";
        std::cerr << std::left << std::setw(34) << "case" << std::right << std::setw(14) << "node abs" << std::setw(14) << "local w abs"
                  << std::setw(14) << "global w abs" << std::setw(14) << "moment abs" << std::setw(14) << "moment rel" << std::setw(18)
                  << "branch counts" << std::setw(9) << "status" << '\n';
        std::cerr << std::string(131, '-') << '\n';
        std::cerr << std::left << std::setw(34) << expected->name << std::right << std::scientific << std::setprecision(4) << std::setw(14)
                  << result.errors.node << std::setw(14) << result.errors.local_weight << std::setw(14) << result.errors.global_weight
                  << std::setw(14) << result.errors.moment_absolute << std::setw(14) << result.errors.moment_relative << std::setw(8);
        print_counts(result.recovered_branch_counts);
        std::cerr << std::setw(9) << (passed ? "PASS" : "FAIL") << '\n';

        if (!result.passed)
        {
            std::cerr << result.message << '\n';
        }
    }

    return passed ? 0 : 1;
}
