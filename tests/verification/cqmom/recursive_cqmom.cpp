#include "../verification_helpers.hpp"
#include "internal/recursive_quadrature.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace
{

struct ValidationMetrics
{
    double max_abs_moment_error = 0.0;
    double max_rel_moment_error = 0.0;
    double max_weight_error = 0.0;
    double max_abscissa_error = 0.0;
};

std::vector<double> raw_moments_from_quadrature(const qmomkit::internal::RecursiveQuadrature& quadrature, const std::vector<std::size_t>& extents)
{
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    const qmomkit::PhaseSpaceHierarchy& hierarchy = quadrature.hierarchy();
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        double value = 0.0;
        for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(hierarchy))
        {
            double contribution = quadrature.weight(node_index);
            for (std::size_t dimension = 0; dimension < hierarchy.number_of_dimensions(); ++dimension)
            {
                contribution *= std::pow(quadrature.abscissa(dimension, node_index), static_cast<int>(moment_index[dimension]));
            }
            value += contribution;
        }
        moments[flat_moment] = value;
    }
    return moments;
}

std::vector<double> raw_moments_from_terminal_quadrature(const std::vector<std::vector<double>>& abscissas, const std::vector<double>& weights,
                                                         const std::vector<std::size_t>& extents)
{
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        double value = 0.0;
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            double contribution = weights[node];
            for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
            {
                contribution *= std::pow(abscissas[dimension][node], static_cast<int>(moment_index[dimension]));
            }
            value += contribution;
        }
        moments[flat_moment] = value;
    }
    return moments;
}

ValidationMetrics validate_case(const std::string& name, const qmomkit::internal::RecursiveQuadrature& expected,
                                const std::vector<std::size_t>& extents, const double moment_tolerance, const double node_tolerance)
{
    const std::vector<double> raw_moments = raw_moments_from_quadrature(expected, extents);
    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature(
        {expected.hierarchy(), {extents, raw_moments}, {}, {}, false, qmomkit::CqmomDiagnosticLevel::basic});
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(!result.diagnostics.empty());
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostic : result.diagnostics)
    {
        QMOMKIT_VERIFY_CHECK(diagnostic.conditioning_level < expected.number_of_dimensions());
        QMOMKIT_VERIFY_CHECK(diagnostic.phase_space_dimension == expected.hierarchy().conditioning_order()[diagnostic.conditioning_level]);
        QMOMKIT_VERIFY_CHECK(diagnostic.wheeler_status == qmomkit::QuadratureStatus::success);
        QMOMKIT_VERIFY_CHECK(diagnostic.conditional_moment_status == qmomkit::ConditionalMomentStatus::success);
    }

    const std::vector<double> reconstructed = raw_moments_from_quadrature(result.quadrature, extents);
    QMOMKIT_VERIFY_CHECK(reconstructed.size() == raw_moments.size());

    ValidationMetrics metrics;
    for (std::size_t i = 0; i < raw_moments.size(); ++i)
    {
        const double abs_error = std::abs(reconstructed[i] - raw_moments[i]);
        const double scale = std::max(1.0, std::abs(raw_moments[i]));
        metrics.max_abs_moment_error = std::max(metrics.max_abs_moment_error, abs_error);
        metrics.max_rel_moment_error = std::max(metrics.max_rel_moment_error, abs_error / scale);
    }

    for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(expected.hierarchy()))
    {
        metrics.max_weight_error = std::max(metrics.max_weight_error, std::abs(result.quadrature.weight(node_index) - expected.weight(node_index)));
        for (std::size_t dimension = 0; dimension < expected.number_of_dimensions(); ++dimension)
        {
            metrics.max_abscissa_error = std::max(
                metrics.max_abscissa_error, std::abs(result.quadrature.abscissa(dimension, node_index) - expected.abscissa(dimension, node_index)));
        }
    }

    std::cout << std::setprecision(17) << name << " max_abs_moment_error=" << metrics.max_abs_moment_error
              << " max_rel_moment_error=" << metrics.max_rel_moment_error << " max_weight_error=" << metrics.max_weight_error
              << " max_abscissa_error=" << metrics.max_abscissa_error << '\n';

    QMOMKIT_VERIFY_CHECK(metrics.max_abs_moment_error <= moment_tolerance);
    QMOMKIT_VERIFY_CHECK(metrics.max_rel_moment_error <= moment_tolerance);
    QMOMKIT_VERIFY_CHECK(metrics.max_weight_error <= node_tolerance);
    QMOMKIT_VERIFY_CHECK(metrics.max_abscissa_error <= node_tolerance);
    return metrics;
}

struct RaggedAdaptiveCase
{
    qmomkit::PhaseSpaceHierarchy hierarchy;
    std::vector<std::size_t> extents;
    std::vector<double> weights;
    std::vector<std::vector<double>> abscissas;
    std::vector<qmomkit::ConditioningPath> paths;
};

ValidationMetrics validate_adaptive_case(const std::string& name, const RaggedAdaptiveCase& expected, const double moment_tolerance,
                                         const double node_tolerance)
{
    const std::vector<double> raw_moments = raw_moments_from_terminal_quadrature(expected.abscissas, expected.weights, expected.extents);
    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature(
        {expected.hierarchy, {expected.extents, raw_moments}, {}, {true, 1.0e-10, 1.0e-10, 1.0e-12}, true, qmomkit::CqmomDiagnosticLevel::basic});
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.quadrature.adaptive_branch_counts_enabled());
    QMOMKIT_VERIFY_CHECK(!result.quadrature.has_rectangular_storage());
    QMOMKIT_VERIFY_CHECK(result.quadrature.total_nodes() == expected.weights.size());
    QMOMKIT_VERIFY_CHECK(result.quadrature.terminal_paths() == expected.paths);

    double accumulated_weight_sum = 0.0;
    for (std::size_t terminal = 0; terminal < result.quadrature.total_nodes(); ++terminal)
    {
        const std::size_t node_id = result.quadrature.terminal_node_ids()[terminal];
        const qmomkit::internal::RecursiveQuadrature::TreeNode& node = result.quadrature.tree_nodes()[node_id];
        QMOMKIT_VERIFY_CHECK(node.terminal);
        QMOMKIT_VERIFY_CHECK(node.path == expected.paths[terminal]);
        QMOMKIT_VERIFY_CHECK(std::abs(node.accumulated_weight - expected.weights[terminal]) <= node_tolerance);
        QMOMKIT_VERIFY_CHECK(std::abs(result.quadrature.weights()[terminal] - expected.weights[terminal]) <= node_tolerance);
        accumulated_weight_sum += result.quadrature.weights()[terminal];
        for (std::size_t dimension = 0; dimension < expected.hierarchy.number_of_dimensions(); ++dimension)
        {
            QMOMKIT_VERIFY_CHECK(std::abs(node.abscissas_by_dimension[dimension] - expected.abscissas[dimension][terminal]) <= node_tolerance);
            QMOMKIT_VERIFY_CHECK(std::abs(result.quadrature.abscissas()[dimension][terminal] - expected.abscissas[dimension][terminal]) <=
                                 node_tolerance);
        }
    }
    QMOMKIT_VERIFY_CHECK(std::abs(accumulated_weight_sum - raw_moments.front()) <= node_tolerance);

    bool saw_reduced_branch = false;
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostic : result.diagnostics)
    {
        QMOMKIT_VERIFY_CHECK(diagnostic.conditioning_level < expected.hierarchy.number_of_dimensions());
        QMOMKIT_VERIFY_CHECK(diagnostic.phase_space_dimension == expected.hierarchy.conditioning_order()[diagnostic.conditioning_level]);
        QMOMKIT_VERIFY_CHECK(diagnostic.requested_node_count == expected.hierarchy.node_count(diagnostic.phase_space_dimension));
        QMOMKIT_VERIFY_CHECK(diagnostic.actual_node_count > 0);
        QMOMKIT_VERIFY_CHECK(diagnostic.actual_node_count <= diagnostic.requested_node_count);
        QMOMKIT_VERIFY_CHECK(diagnostic.adaptive_branch_counts_enabled);
        QMOMKIT_VERIFY_CHECK(diagnostic.wheeler_status == qmomkit::QuadratureStatus::success);
        if (diagnostic.adaptive_node_count_reduced)
        {
            saw_reduced_branch = true;
            QMOMKIT_VERIFY_CHECK(diagnostic.actual_node_count < diagnostic.requested_node_count);
        }
        else
        {
            QMOMKIT_VERIFY_CHECK(diagnostic.actual_node_count == diagnostic.requested_node_count);
        }
    }
    QMOMKIT_VERIFY_CHECK(saw_reduced_branch);

    const std::vector<double> reconstructed =
        raw_moments_from_terminal_quadrature(result.quadrature.abscissas(), result.quadrature.weights(), expected.extents);
    QMOMKIT_VERIFY_CHECK(reconstructed.size() == raw_moments.size());

    ValidationMetrics metrics;
    for (std::size_t i = 0; i < raw_moments.size(); ++i)
    {
        const double abs_error = std::abs(reconstructed[i] - raw_moments[i]);
        const double scale = std::max(1.0, std::abs(raw_moments[i]));
        metrics.max_abs_moment_error = std::max(metrics.max_abs_moment_error, abs_error);
        metrics.max_rel_moment_error = std::max(metrics.max_rel_moment_error, abs_error / scale);
    }
    for (std::size_t terminal = 0; terminal < result.quadrature.total_nodes(); ++terminal)
    {
        metrics.max_weight_error = std::max(metrics.max_weight_error, std::abs(result.quadrature.weights()[terminal] - expected.weights[terminal]));
        for (std::size_t dimension = 0; dimension < expected.hierarchy.number_of_dimensions(); ++dimension)
        {
            metrics.max_abscissa_error = std::max(
                metrics.max_abscissa_error, std::abs(result.quadrature.abscissas()[dimension][terminal] - expected.abscissas[dimension][terminal]));
        }
    }

    std::cout << std::setprecision(17) << name << " max_abs_moment_error=" << metrics.max_abs_moment_error
              << " max_rel_moment_error=" << metrics.max_rel_moment_error << " max_weight_error=" << metrics.max_weight_error
              << " max_abscissa_error=" << metrics.max_abscissa_error << '\n';

    QMOMKIT_VERIFY_CHECK(metrics.max_abs_moment_error <= moment_tolerance);
    QMOMKIT_VERIFY_CHECK(metrics.max_rel_moment_error <= moment_tolerance);
    QMOMKIT_VERIFY_CHECK(metrics.max_weight_error <= node_tolerance);
    QMOMKIT_VERIFY_CHECK(metrics.max_abscissa_error <= node_tolerance);
    return metrics;
}

ValidationMetrics validate_moment_closure_case(const std::string& name, const qmomkit::internal::RecursiveQuadrature& expected,
                                               const std::vector<std::size_t>& extents, const double moment_tolerance)
{
    const std::vector<double> raw_moments = raw_moments_from_quadrature(expected, extents);
    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature(
        {expected.hierarchy(), {extents, raw_moments}, {}, {}, false, qmomkit::CqmomDiagnosticLevel::basic});
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(!result.diagnostics.empty());
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostic : result.diagnostics)
    {
        QMOMKIT_VERIFY_CHECK(diagnostic.conditioning_level < expected.number_of_dimensions());
        QMOMKIT_VERIFY_CHECK(diagnostic.phase_space_dimension == expected.hierarchy().conditioning_order()[diagnostic.conditioning_level]);
        QMOMKIT_VERIFY_CHECK(diagnostic.wheeler_status == qmomkit::QuadratureStatus::success);
        QMOMKIT_VERIFY_CHECK(diagnostic.conditional_moment_status == qmomkit::ConditionalMomentStatus::success);
    }

    const std::vector<double> reconstructed = raw_moments_from_quadrature(result.quadrature, extents);
    QMOMKIT_VERIFY_CHECK(reconstructed.size() == raw_moments.size());

    ValidationMetrics metrics;
    for (std::size_t i = 0; i < raw_moments.size(); ++i)
    {
        const double abs_error = std::abs(reconstructed[i] - raw_moments[i]);
        const double scale = std::max(1.0, std::abs(raw_moments[i]));
        metrics.max_abs_moment_error = std::max(metrics.max_abs_moment_error, abs_error);
        metrics.max_rel_moment_error = std::max(metrics.max_rel_moment_error, abs_error / scale);
    }

    std::cout << std::setprecision(17) << name << " max_abs_moment_error=" << metrics.max_abs_moment_error
              << " max_rel_moment_error=" << metrics.max_rel_moment_error << '\n';

    QMOMKIT_VERIFY_CHECK(metrics.max_abs_moment_error <= moment_tolerance);
    QMOMKIT_VERIFY_CHECK(metrics.max_rel_moment_error <= moment_tolerance);
    return metrics;
}

std::vector<std::size_t> moment_extents_for_node_counts(const std::vector<std::size_t>& node_counts)
{
    std::vector<std::size_t> extents;
    extents.reserve(node_counts.size());
    for (const std::size_t node_count : node_counts)
    {
        extents.push_back(2 * node_count);
    }
    return extents;
}

double path_signature(const std::vector<std::size_t>& path)
{
    double signature = 0.0;
    for (std::size_t level = 0; level < path.size(); ++level)
    {
        signature += static_cast<double>((level + 1) * (path[level] + 1));
    }
    return signature;
}

std::vector<double> sweep_local_abscissas(const std::size_t dimension, const std::size_t level, const std::size_t node_count,
                                          const std::vector<std::size_t>& path)
{
    const double signature = path_signature(path);
    const double center = -0.55 + 0.18 * static_cast<double>(dimension) - 0.04 * static_cast<double>(level) + 0.015 * signature;
    std::vector<double> values;
    values.reserve(node_count);
    for (std::size_t node = 0; node < node_count; ++node)
    {
        values.push_back(center + 0.55 * static_cast<double>(node) + 0.03 * static_cast<double>(node * node));
    }
    return values;
}

std::vector<double> sweep_local_weights(const std::size_t dimension, const std::size_t level, const std::size_t node_count,
                                        const std::vector<std::size_t>& path)
{
    const double signature = path_signature(path);
    std::vector<double> weights;
    weights.reserve(node_count);
    double total = 0.0;
    for (std::size_t node = 0; node < node_count; ++node)
    {
        const double raw_weight =
            0.35 + 0.11 * static_cast<double>(node + 1) + 0.025 * static_cast<double>((dimension + node + level) % 3) + 0.006 * signature;
        weights.push_back(raw_weight);
        total += raw_weight;
    }
    for (double& weight : weights)
    {
        weight /= total;
    }
    return weights;
}

void populate_dimension_sweep_recursive(qmomkit::internal::RecursiveQuadrature& quadrature, const std::size_t level, std::vector<std::size_t>& path,
                                        qmomkit::MultiIndex& node_index, std::vector<double>& abscissas, const double accumulated_weight)
{
    const qmomkit::PhaseSpaceHierarchy& hierarchy = quadrature.hierarchy();
    if (level == hierarchy.number_of_dimensions())
    {
        quadrature.set_weight(node_index, accumulated_weight);
        for (std::size_t dimension = 0; dimension < hierarchy.number_of_dimensions(); ++dimension)
        {
            quadrature.set_abscissa(dimension, node_index, abscissas[dimension]);
        }
        return;
    }

    const std::size_t dimension = hierarchy.conditioning_order()[level];
    const std::size_t node_count = hierarchy.node_count(dimension);
    const std::vector<double> local_abscissas = sweep_local_abscissas(dimension, level, node_count, path);
    const std::vector<double> local_weights = sweep_local_weights(dimension, level, node_count, path);
    for (std::size_t node = 0; node < node_count; ++node)
    {
        path.push_back(node);
        node_index[dimension] = node;
        abscissas[dimension] = local_abscissas[node];
        populate_dimension_sweep_recursive(quadrature, level + 1, path, node_index, abscissas, accumulated_weight * local_weights[node]);
        path.pop_back();
    }
}

qmomkit::internal::RecursiveQuadrature make_dimension_sweep_case(const std::vector<std::size_t>& node_counts,
                                                                 const std::vector<std::size_t>& conditioning_order)
{
    qmomkit::internal::RecursiveQuadrature quadrature(qmomkit::PhaseSpaceHierarchy(node_counts, conditioning_order));
    std::vector<std::size_t> path;
    qmomkit::MultiIndex node_index(node_counts.size(), 0);
    std::vector<double> abscissas(node_counts.size(), 0.0);
    populate_dimension_sweep_recursive(quadrature, 0, path, node_index, abscissas, 1.0);
    return quadrature;
}

qmomkit::internal::RecursiveQuadrature make_1d_case()
{
    qmomkit::internal::RecursiveQuadrature quadrature(qmomkit::PhaseSpaceHierarchy({2}));
    quadrature.set_weight({0}, 0.35);
    quadrature.set_abscissa(0, {0}, -1.25);
    quadrature.set_weight({1}, 0.65);
    quadrature.set_abscissa(0, {1}, 2.0);
    return quadrature;
}

qmomkit::internal::RecursiveQuadrature make_2d_identity_case()
{
    qmomkit::internal::RecursiveQuadrature quadrature(qmomkit::PhaseSpaceHierarchy({2, 2}));
    const std::vector<double> x{-1.0, 2.0};
    const std::vector<double> wx{0.25, 0.75};
    const std::vector<std::vector<double>> y{{-3.0, 1.0}, {-2.0, 4.0}};
    const std::vector<std::vector<double>> wy{{0.2, 0.8}, {0.6, 0.4}};
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 2; ++j)
        {
            quadrature.set_weight({i, j}, wx[i] * wy[i][j]);
            quadrature.set_abscissa(0, {i, j}, x[i]);
            quadrature.set_abscissa(1, {i, j}, y[i][j]);
        }
    }
    return quadrature;
}

qmomkit::internal::RecursiveQuadrature make_2d_reversed_case()
{
    qmomkit::internal::RecursiveQuadrature quadrature(qmomkit::PhaseSpaceHierarchy({2, 2}, {1, 0}));
    const std::vector<double> y{-0.5, 1.25};
    const std::vector<double> wy{0.2, 0.8};
    const std::vector<std::vector<double>> x{{-4.0, -1.0}, {0.25, 2.0}};
    const std::vector<std::vector<double>> wx{{0.35, 0.65}, {0.7, 0.3}};
    for (std::size_t y_node = 0; y_node < 2; ++y_node)
    {
        for (std::size_t x_node = 0; x_node < 2; ++x_node)
        {
            const qmomkit::MultiIndex node_index{x_node, y_node};
            quadrature.set_weight(node_index, wy[y_node] * wx[y_node][x_node]);
            quadrature.set_abscissa(0, node_index, x[y_node][x_node]);
            quadrature.set_abscissa(1, node_index, y[y_node]);
        }
    }
    return quadrature;
}

qmomkit::internal::RecursiveQuadrature make_3d_non_identity_case()
{
    qmomkit::internal::RecursiveQuadrature quadrature(qmomkit::PhaseSpaceHierarchy({2, 3, 2}, {2, 0, 1}));
    const std::vector<double> z{-1.5, 0.75};
    const std::vector<double> wz{0.45, 0.55};
    const std::vector<std::vector<double>> x{{-3.0, 0.5}, {-2.0, 4.0}};
    const std::vector<std::vector<double>> wx{{0.25, 0.75}, {0.4, 0.6}};
    for (std::size_t z_node = 0; z_node < 2; ++z_node)
    {
        for (std::size_t x_node = 0; x_node < 2; ++x_node)
        {
            const std::vector<double> y{-2.0 + static_cast<double>(z_node), 0.25 + 0.5 * static_cast<double>(x_node),
                                        2.5 + static_cast<double>(z_node) + static_cast<double>(x_node)};
            const std::vector<double> wy{0.2 + 0.05 * static_cast<double>(x_node), 0.3, 0.5 - 0.05 * static_cast<double>(x_node)};
            for (std::size_t y_node = 0; y_node < 3; ++y_node)
            {
                const qmomkit::MultiIndex node_index{x_node, y_node, z_node};
                quadrature.set_weight(node_index, wz[z_node] * wx[z_node][x_node] * wy[y_node]);
                quadrature.set_abscissa(0, node_index, x[z_node][x_node]);
                quadrature.set_abscissa(1, node_index, y[y_node]);
                quadrature.set_abscissa(2, node_index, z[z_node]);
            }
        }
    }
    return quadrature;
}

qmomkit::internal::RecursiveQuadrature make_non_uniform_case()
{
    qmomkit::internal::RecursiveQuadrature quadrature(qmomkit::PhaseSpaceHierarchy({2, 3}));
    const std::vector<double> x{-2.0, 0.5};
    const std::vector<double> wx{0.3, 0.7};
    const std::vector<std::vector<double>> y{{-2.0, 0.0, 3.0}, {-1.0, 1.0, 4.0}};
    const std::vector<std::vector<double>> wy{{0.2, 0.5, 0.3}, {0.25, 0.25, 0.5}};
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 3; ++j)
        {
            quadrature.set_weight({i, j}, wx[i] * wy[i][j]);
            quadrature.set_abscissa(0, {i, j}, x[i]);
            quadrature.set_abscissa(1, {i, j}, y[i][j]);
        }
    }
    return quadrature;
}

qmomkit::internal::RecursiveQuadrature make_4d_adversarial_non_identity_case()
{
    qmomkit::internal::RecursiveQuadrature quadrature(qmomkit::PhaseSpaceHierarchy({2, 3, 2, 2}, {2, 0, 3, 1}));
    const std::vector<double> z{-1.1, 0.9};
    const std::vector<double> wz{0.91, 0.09};
    const std::vector<std::vector<double>> x{{-0.015, 0.035}, {-2.2, 1.7}};
    const std::vector<std::vector<double>> wx{{0.12, 0.88}, {0.68, 0.32}};

    for (std::size_t z_node = 0; z_node < 2; ++z_node)
    {
        for (std::size_t x_node = 0; x_node < 2; ++x_node)
        {
            const std::vector<double> u{-1.0 + 0.30 * static_cast<double>(z_node) + 0.20 * static_cast<double>(x_node),
                                        0.80 + 0.10 * static_cast<double>(z_node) + 0.25 * static_cast<double>(x_node)};
            const std::vector<double> wu{0.18 + 0.07 * static_cast<double>(z_node) + 0.04 * static_cast<double>(x_node),
                                         0.82 - 0.07 * static_cast<double>(z_node) - 0.04 * static_cast<double>(x_node)};
            for (std::size_t u_node = 0; u_node < 2; ++u_node)
            {
                const double first_y_weight = 0.06 + 0.01 * static_cast<double>(z_node) + 0.015 * static_cast<double>(u_node);
                const double second_y_weight = 0.24 + 0.04 * static_cast<double>(x_node) - 0.02 * static_cast<double>(u_node);
                const std::vector<double> wy{first_y_weight, second_y_weight, 1.0 - first_y_weight - second_y_weight};
                const std::vector<double> y{
                    -2.5 + 0.40 * static_cast<double>(z_node) - 0.25 * static_cast<double>(x_node) + 0.15 * static_cast<double>(u_node),
                    -0.1 + 0.20 * static_cast<double>(z_node) + 0.10 * static_cast<double>(x_node) - 0.05 * static_cast<double>(u_node),
                    1.6 + 0.30 * static_cast<double>(z_node) + 0.35 * static_cast<double>(x_node) + 0.20 * static_cast<double>(u_node)};
                for (std::size_t y_node = 0; y_node < 3; ++y_node)
                {
                    const qmomkit::MultiIndex node_index{x_node, y_node, z_node, u_node};
                    quadrature.set_weight(node_index, wz[z_node] * wx[z_node][x_node] * wu[u_node] * wy[y_node]);
                    quadrature.set_abscissa(0, node_index, x[z_node][x_node]);
                    quadrature.set_abscissa(1, node_index, y[y_node]);
                    quadrature.set_abscissa(2, node_index, z[z_node]);
                    quadrature.set_abscissa(3, node_index, u[u_node]);
                }
            }
        }
    }
    return quadrature;
}

void assert_4d_adversarial_phase_order_storage(const qmomkit::internal::RecursiveQuadrature& quadrature)
{
    // The final entry is valid for conditioning level 3 (phase-space dimension 1),
    // but would be out of range for phase-space dimension 3 under identity-order mapping.
    const qmomkit::ConditioningPath path({1, 0, 1, 2});
    const qmomkit::MultiIndex node_index = qmomkit::recursive_path_to_node_index(quadrature.hierarchy(), path);
    QMOMKIT_VERIFY_CHECK(node_index == qmomkit::MultiIndex({0, 2, 1, 1}));

    const std::vector<double> node = quadrature.node(node_index);
    QMOMKIT_VERIFY_CHECK(node.size() == 4);
    QMOMKIT_VERIFY_CHECK(std::abs(node[0] - (-2.2)) <= 1.0e-8);
    QMOMKIT_VERIFY_CHECK(std::abs(node[1] - 2.1) <= 1.0e-8);
    QMOMKIT_VERIFY_CHECK(std::abs(node[2] - 0.9) <= 1.0e-8);
    QMOMKIT_VERIFY_CHECK(std::abs(node[3] - 0.9) <= 1.0e-8);
}

RaggedAdaptiveCase make_4d_adversarial_adaptive_case()
{
    RaggedAdaptiveCase test_case{qmomkit::PhaseSpaceHierarchy({2, 2, 2, 2}, {2, 0, 3, 1}), {4, 4, 4, 4}, {}, std::vector<std::vector<double>>(4), {}};

    const auto add_terminal =
        [&test_case](const qmomkit::ConditioningPath& path, const double weight, const double x, const double y, const double z, const double u)
    {
        test_case.paths.push_back(path);
        test_case.weights.push_back(weight);
        test_case.abscissas[0].push_back(x);
        test_case.abscissas[1].push_back(y);
        test_case.abscissas[2].push_back(z);
        test_case.abscissas[3].push_back(u);
    };

    add_terminal(qmomkit::ConditioningPath({0, 0, 0, 0}), 0.63 * 0.20 * 0.70 * 0.25, -0.8, -2.0, -1.25, -1.1);
    add_terminal(qmomkit::ConditioningPath({0, 0, 0, 1}), 0.63 * 0.20 * 0.70 * 0.75, -0.8, 0.5, -1.25, -1.1);
    add_terminal(qmomkit::ConditioningPath({0, 0, 1, 0}), 0.63 * 0.20 * 0.30, -0.8, 1.1, -1.25, 0.6);
    add_terminal(qmomkit::ConditioningPath({0, 1, 0, 0}), 0.63 * 0.80 * 0.60, 1.4, -0.4, -1.25, 1.2);
    add_terminal(qmomkit::ConditioningPath({0, 1, 0, 1}), 0.63 * 0.80 * 0.40, 1.4, 2.2, -1.25, 1.2);
    add_terminal(qmomkit::ConditioningPath({1, 0, 0, 0}), 0.37 * 0.45, 0.35, -1.5, 0.75, -0.2);
    add_terminal(qmomkit::ConditioningPath({1, 0, 1, 0}), 0.37 * 0.55 * 0.33, 0.35, 0.2, 0.75, 1.5);
    add_terminal(qmomkit::ConditioningPath({1, 0, 1, 1}), 0.37 * 0.55 * 0.67, 0.35, 1.8, 0.75, 1.5);
    return test_case;
}

} // namespace

int main(const int argc, char** argv)
{
    const std::string selected_case = argc == 2 ? argv[1] : "";
    if (argc > 2)
    {
        std::cerr << "Usage: recursive_cqmom [case-name]\n";
        return 2;
    }
    const auto selected = [&selected_case](const std::string& name) { return selected_case.empty() || selected_case == name; };
    bool matched = selected_case.empty();
    if (selected("1d"))
    {
        matched = true;
        (void)validate_case("1d", make_1d_case(), {4}, 1.0e-12, 1.0e-12);
    }
    if (selected("2d_identity"))
    {
        matched = true;
        (void)validate_case("2d_identity", make_2d_identity_case(), {4, 4}, 1.0e-10, 1.0e-10);
    }
    if (selected("2d_reversed"))
    {
        matched = true;
        (void)validate_case("2d_reversed", make_2d_reversed_case(), {4, 4}, 1.0e-10, 1.0e-10);
    }
    if (selected("3d_non_identity"))
    {
        matched = true;
        (void)validate_case("3d_non_identity", make_3d_non_identity_case(), {4, 6, 4}, 1.0e-7, 1.0e-8);
    }
    if (selected("non_uniform_fixed_counts"))
    {
        matched = true;
        (void)validate_case("non_uniform_fixed_counts", make_non_uniform_case(), {4, 6}, 1.0e-8, 1.0e-8);
    }
    const qmomkit::internal::RecursiveQuadrature adversarial_4d = make_4d_adversarial_non_identity_case();
    if (selected("4d_adversarial_non_identity"))
    {
        matched = true;
        (void)validate_case("4d_adversarial_non_identity", adversarial_4d, {4, 6, 4, 4}, 1.0e-9, 1.0e-8);
        const std::vector<double> adversarial_moments = raw_moments_from_quadrature(adversarial_4d, {4, 6, 4, 4});
        const qmomkit::internal::RecursiveQuadratureResult adversarial_result = qmomkit::internal::compute_recursive_quadrature(
            {adversarial_4d.hierarchy(), {{4, 6, 4, 4}, adversarial_moments}, {}, {}, false, qmomkit::CqmomDiagnosticLevel::basic});
        QMOMKIT_VERIFY_CHECK(adversarial_result.ok());
        assert_4d_adversarial_phase_order_storage(adversarial_result.quadrature);
    }
    if (selected("4d_adversarial_adaptive_non_identity"))
    {
        matched = true;
        (void)validate_adaptive_case("4d_adversarial_adaptive_non_identity", make_4d_adversarial_adaptive_case(), 1.0e-9, 1.0e-8);
    }

    const std::vector<std::vector<std::size_t>> sweep_node_counts{{2, 3}, {2, 2, 3}, {2, 3, 2, 2}, {2, 2, 3, 2, 2}, {2, 2, 2, 3, 2, 2}};
    const std::vector<std::vector<std::size_t>> sweep_orders{{1, 0}, {2, 0, 1}, {2, 0, 3, 1}, {4, 2, 0, 3, 1}, {5, 2, 0, 4, 1, 3}};
    for (std::size_t case_index = 0; case_index < sweep_node_counts.size(); ++case_index)
    {
        const std::string name = "dimension_sweep_D" + std::to_string(sweep_node_counts[case_index].size());
        if (selected(name))
        {
            matched = true;
            const qmomkit::internal::RecursiveQuadrature sweep_case =
                make_dimension_sweep_case(sweep_node_counts[case_index], sweep_orders[case_index]);
            (void)validate_moment_closure_case(name, sweep_case, moment_extents_for_node_counts(sweep_node_counts[case_index]), 1.0e-12);
        }
    }
    if (!matched)
    {
        std::cerr << "Unknown recursive CQMOM case '" << selected_case << "'.\n";
        return 2;
    }
    return 0;
}
