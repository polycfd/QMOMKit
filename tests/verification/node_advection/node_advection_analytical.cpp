// Moving quadrature positions along a velocity field, against trajectories
// that can be written down in closed form.
//
// What is verified
//   That advancing the positions reproduces the exact trajectory, and that
//   advection changes positions only - weights, hierarchy and terminal paths
//   come through untouched, which is what makes the operation a transport step
//   rather than a new inversion.
//
// How the reference is constructed
//   The velocity fields are chosen so the trajectory is solvable on paper, and
//   the exact position is evaluated directly from that solution:
//
//     constant_velocity  dx/dt = 0.6          x(t) = x0 + 0.6*t
//     linear_growth      dx/dt = 0.55*x       x(t) = x0*exp(0.55*t)
//     affine_evolution   dx/dt = -0.4*x + 0.7 x(t) = (x0 + b/a)*exp(a*t) - b/a
//                                             with a = -0.4 and b = 0.7
//     linear_system_2d   dx/dt = -0.3*x, dy/dt = 0.8*x + 0.2*y, solved in
//                        exponentials of -0.3*t and 0.2*t
//
//   Nothing in the reference comes from an integrator.
//
// Error metric and pass criteria
//   Maximum absolute deviation of each position from its exact trajectory, at
//   5.0e-8. That is looser than the inversion studies because the positions are
//   produced by a time integrator asked for a relative accuracy of 1.0e-10, and
//   its error accumulates over the interval; the bound sits above what the
//   controller was asked to deliver rather than at machine precision. Weights
//   are held to the same 5.0e-8 and the terminal node count must not change.
//   Moments rebuilt from the advanced positions are compared with those of the
//   exact positions, unit-floored relative at 5.0e-8 and absolute at 1.0e-5.
//
// Reproduce
//   ctest -R verification_node_advection
//
// Build modes
//   Both. Failures are reported through an exit code, not assert().

#include "qmomkit/configuration.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/node_advection.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace
{

constexpr double tolerance = 5.0e-8;

struct ErrorMetrics
{
    double node = 0.0;
    double weight = 0.0;
    double moment_absolute = 0.0;
    double moment_relative = 0.0;
};

struct TimePlan
{
    double t0 = 0.0;
    double dt = 0.0;
    std::size_t steps = 0;
};

struct CaseDefinition
{
    std::string name;
    bool cqmom = false;
    TimePlan one_step;
    TimePlan multi_step;
};

struct VerificationResult
{
    bool passed = true;
    std::string message;
    ErrorMetrics errors;
    double loose_error = 0.0;
    double tight_error = 0.0;
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
    for (std::size_t power = 0; power < exponent; ++power)
    {
        value *= base;
    }
    return value;
}

std::vector<std::pair<double, double>> sort_rule(std::vector<std::pair<double, double>> rule)
{
    std::sort(rule.begin(), rule.end(), [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });
    return rule;
}

bool coordinate_less(const std::vector<double>& lhs, const std::vector<double>& rhs)
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

qmomkit::NodeAdvectionOptions accurate_options()
{
    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    options.adaptive_options.relative_tolerance = 1.0e-10;
    options.adaptive_options.absolute_tolerance = 1.0e-12;
    options.adaptive_options.initial_step = 0.01;
    options.adaptive_options.max_step = 0.05;
    return options;
}

const CaseDefinition* find_case(const std::string& name)
{
    static const std::vector<CaseDefinition> cases{{"constant_velocity", false, {0.2, 0.15, 1}, {0.2, 0.12, 7}},
                                                   {"linear_growth", false, {-0.1, 0.2, 1}, {-0.1, 0.1, 8}},
                                                   {"affine_evolution", false, {0.4, 0.2, 1}, {0.4, 0.15, 6}},
                                                   {"linear_system_2d", true, {0.3, 0.1, 1}, {0.3, 0.08, 8}}};
    for (const CaseDefinition& test_case : cases)
    {
        if (test_case.name == name)
        {
            return &test_case;
        }
    }
    return nullptr;
}

std::vector<double> exact_qmom_nodes(const std::string& name, const std::vector<double>& initial_nodes, const double elapsed)
{
    std::vector<double> exact(initial_nodes.size(), 0.0);
    if (name == "constant_velocity")
    {
        constexpr double velocity = 0.6;
        for (std::size_t node = 0; node < initial_nodes.size(); ++node)
        {
            exact[node] = initial_nodes[node] + velocity * elapsed;
        }
    }
    else if (name == "linear_growth")
    {
        constexpr double rate = 0.55;
        const double factor = std::exp(rate * elapsed);
        for (std::size_t node = 0; node < initial_nodes.size(); ++node)
        {
            exact[node] = initial_nodes[node] * factor;
        }
    }
    else if (name == "affine_evolution")
    {
        constexpr double rate = -0.4;
        constexpr double shift = 0.7;
        const double factor = std::exp(rate * elapsed);
        for (std::size_t node = 0; node < initial_nodes.size(); ++node)
        {
            exact[node] = (initial_nodes[node] + shift / rate) * factor - shift / rate;
        }
    }
    return exact;
}

std::vector<double> qmom_velocity(const std::string& name, const double, const std::vector<double>& xi)
{
    if (name == "constant_velocity")
    {
        return {0.6};
    }
    if (name == "linear_growth")
    {
        return {0.55 * xi[0]};
    }
    if (name == "affine_evolution")
    {
        return {-0.4 * xi[0] + 0.7};
    }
    return {};
}

std::vector<double> exact_qmom_moments(const std::vector<double>& nodes, const std::vector<double>& weights, const std::size_t count)
{
    std::vector<double> moments(count, 0.0);
    for (std::size_t order = 0; order < count; ++order)
    {
        for (std::size_t node = 0; node < nodes.size(); ++node)
        {
            moments[order] += weights[node] * integer_power(nodes[node], order);
        }
    }
    return moments;
}

std::vector<std::pair<double, double>> expected_initial_qmom_rule(const std::string& name)
{
    if (name == "constant_velocity")
    {
        return {{-0.8, 0.2}, {0.4, 0.35}, {1.7, 0.45}};
    }
    if (name == "linear_growth")
    {
        return {{-1.2, 0.15}, {0.3, 0.55}, {1.4, 0.30}};
    }
    if (name == "affine_evolution")
    {
        return {{0.2, 0.25}, {1.6, 0.45}, {2.8, 0.30}};
    }
    return {};
}

void compare_qmom_state(const qmomkit::Quadrature1D& quadrature, const std::vector<double>& exact_nodes, const std::vector<double>& initial_weights,
                        const std::size_t moment_count, VerificationResult& result)
{
    std::vector<std::pair<double, double>> actual;
    std::vector<std::pair<double, double>> expected;
    for (std::size_t node = 0; node < exact_nodes.size(); ++node)
    {
        actual.emplace_back(quadrature.abscissas[node], quadrature.weights[node]);
        expected.emplace_back(exact_nodes[node], initial_weights[node]);
    }
    actual = sort_rule(std::move(actual));
    expected = sort_rule(std::move(expected));
    for (std::size_t node = 0; node < expected.size(); ++node)
    {
        update_error(result.errors.node, actual[node].first, expected[node].first);
        update_error(result.errors.weight, actual[node].second, expected[node].second);
        if (!(actual[node].second > 0.0))
        {
            append_failure(result, "QMOM weight is not positive.");
        }
    }

    const std::vector<double> reconstructed = qmomkit::reconstruct_raw_moments(quadrature, moment_count);
    const std::vector<double> exact_moments = exact_qmom_moments(exact_nodes, initial_weights, moment_count);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        update_error(result.errors.moment_absolute, reconstructed[order], exact_moments[order]);
        update_relative_error(result.errors.moment_relative, reconstructed[order], exact_moments[order]);
    }
    if (std::abs(reconstructed.front() - exact_moments.front()) > tolerance)
    {
        append_failure(result, "QMOM m0 was not preserved.");
    }
}

void verify_initial_qmom_state(const CaseDefinition& definition, const qmomkit::QmomCase& qmom_case, const qmomkit::Quadrature1D& quadrature,
                               VerificationResult& result)
{
    const std::vector<std::pair<double, double>> expected_rule = expected_initial_qmom_rule(definition.name);
    if (expected_rule.empty())
    {
        append_failure(result, "No independent initial QMOM rule is defined for this case.");
        return;
    }

    std::vector<double> expected_nodes;
    std::vector<double> expected_weights;
    expected_nodes.reserve(expected_rule.size());
    expected_weights.reserve(expected_rule.size());
    for (const auto& [node, weight] : expected_rule)
    {
        expected_nodes.push_back(node);
        expected_weights.push_back(weight);
    }

    compare_qmom_state(quadrature, expected_nodes, expected_weights, qmom_case.moments.size(), result);
    const std::vector<double> expected_moments = exact_qmom_moments(expected_nodes, expected_weights, qmom_case.moments.size());
    for (std::size_t order = 0; order < qmom_case.moments.size(); ++order)
    {
        update_error(result.errors.moment_absolute, qmom_case.moments[order], expected_moments[order]);
        update_relative_error(result.errors.moment_relative, qmom_case.moments[order], expected_moments[order]);
    }
}

void verify_qmom_plan(const CaseDefinition& definition, const qmomkit::Quadrature1D& initial, const TimePlan& plan,
                      const qmomkit::NodeAdvectionOptions& options, VerificationResult& result)
{
    qmomkit::Quadrature1D current = initial;
    const std::vector<double> initial_nodes = initial.abscissas;
    const std::vector<double> initial_weights = initial.weights;
    double time = plan.t0;
    for (std::size_t step = 0; step < plan.steps; ++step)
    {
        const qmomkit::QmomNodeAdvectionResult advanced =
            qmomkit::advance_nodes(current, time, plan.dt,
                                   qmomkit::make_phase_space_velocity([&definition](const double t, const std::vector<double>& xi)
                                                                      { return qmom_velocity(definition.name, t, xi); }),
                                   options);
        if (!advanced.ok())
        {
            append_failure(result, "QMOM node advection failed: " + advanced.message);
            return;
        }
        if (advanced.diagnostics.adaptive_status != qmomkit::AdaptiveRkStatus::success || advanced.diagnostics.step_count == 0)
        {
            append_failure(result, "QMOM RK23 diagnostics did not report successful integration.");
        }
        update_error(result.errors.node, advanced.diagnostics.final_time, time + plan.dt);
        for (std::size_t node = 0; node < initial_weights.size(); ++node)
        {
            update_error(result.errors.weight, advanced.quadrature.weights[node], initial_weights[node]);
        }
        current = advanced.quadrature;
        time += plan.dt;
    }
    compare_qmom_state(current, exact_qmom_nodes(definition.name, initial_nodes, time - plan.t0), initial_weights, 6, result);
}

void verify_linear_growth_convergence(const qmomkit::Quadrature1D& initial, VerificationResult& result)
{
    const TimePlan plan{0.0, 0.9, 1};
    qmomkit::NodeAdvectionOptions loose;
    loose.scheme = qmomkit::NodeAdvectionScheme::rk23;
    loose.adaptive_options.relative_tolerance = 1.0e-3;
    loose.adaptive_options.absolute_tolerance = 1.0e-6;
    loose.adaptive_options.initial_step = 0.1;
    loose.adaptive_options.max_step = 0.9;
    qmomkit::NodeAdvectionOptions tight = accurate_options();
    tight.adaptive_options.max_step = 0.9;

    auto run = [&](const qmomkit::NodeAdvectionOptions& options)
    {
        return qmomkit::advance_nodes(
            initial, plan.t0, plan.dt,
            qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>{0.55 * xi[0]}; }), options);
    };
    const auto loose_result = run(loose);
    const auto tight_result = run(tight);
    if (!loose_result.ok() || !tight_result.ok())
    {
        append_failure(result, "Linear-growth convergence run failed.");
        return;
    }
    const std::vector<double> exact = exact_qmom_nodes("linear_growth", initial.abscissas, plan.dt);
    for (std::size_t node = 0; node < exact.size(); ++node)
    {
        result.loose_error = std::max(result.loose_error, std::abs(loose_result.quadrature.abscissas[node] - exact[node]));
        result.tight_error = std::max(result.tight_error, std::abs(tight_result.quadrature.abscissas[node] - exact[node]));
    }
    if (!(result.tight_error < result.loose_error))
    {
        append_failure(result, "Tighter RK23 settings did not reduce the linear-growth final node error.");
    }
}

VerificationResult verify_qmom_case(const CaseDefinition& definition, const std::filesystem::path& configuration_path)
{
    VerificationResult result;
    const auto configured = qmomkit::config::load_qmom(configuration_path);
    if (!configured.ok())
    {
        append_failure(result, "Failed to load QMOM configuration.");
        return result;
    }
    const qmomkit::QmomResult qmom = qmomkit::compute_qmom(configured.value());
    if (!qmom.ok())
    {
        append_failure(result, "Initial QMOM inversion failed: " + qmom.message);
        return result;
    }
    const qmomkit::Quadrature1D initial{qmom.abscissas, qmom.weights, qmomkit::QuadratureStatus::success, qmom.message};
    verify_initial_qmom_state(definition, configured.value(), initial, result);
    verify_qmom_plan(definition, initial, definition.one_step, accurate_options(), result);
    verify_qmom_plan(definition, initial, definition.multi_step, accurate_options(), result);
    if (definition.name == "linear_growth")
    {
        verify_linear_growth_convergence(initial, result);
    }
    return result;
}

struct TerminalNode
{
    std::vector<double> coordinates;
    double weight = 0.0;
};

std::vector<TerminalNode> expected_initial_linear_system_nodes()
{
    return {{{-0.7, -0.4}, 0.45}, {{0.2, 1.1}, 0.33}, {{1.4, 1.1}, 0.22}};
}

std::vector<TerminalNode> terminal_nodes(const qmomkit::CqmomQuadrature& quadrature)
{
    std::vector<TerminalNode> nodes;
    if (quadrature.has_rectangular_storage())
    {
        for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(quadrature.hierarchy()))
        {
            std::vector<double> coordinates(quadrature.number_of_dimensions(), 0.0);
            for (std::size_t dimension = 0; dimension < coordinates.size(); ++dimension)
            {
                coordinates[dimension] = quadrature.abscissa(dimension, node_index);
            }
            nodes.push_back({std::move(coordinates), quadrature.weight(node_index)});
        }
    }
    else
    {
        for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
        {
            std::vector<double> coordinates(quadrature.number_of_dimensions(), 0.0);
            for (std::size_t dimension = 0; dimension < coordinates.size(); ++dimension)
            {
                coordinates[dimension] = quadrature.abscissas()[dimension][node];
            }
            nodes.push_back({std::move(coordinates), quadrature.weights()[node]});
        }
    }
    return nodes;
}

std::vector<TerminalNode> exact_linear_system_nodes(const std::vector<TerminalNode>& initial, const double elapsed)
{
    constexpr double a = -0.3;
    constexpr double b = 0.8;
    constexpr double c = 0.2;
    const double exp_a = std::exp(a * elapsed);
    const double exp_c = std::exp(c * elapsed);
    std::vector<TerminalNode> exact;
    for (const TerminalNode& node : initial)
    {
        const double x0 = node.coordinates[0];
        const double y0 = node.coordinates[1];
        const double x = x0 * exp_a;
        const double y = y0 * exp_c + b * x0 * (exp_a - exp_c) / (a - c);
        exact.push_back({{x, y}, node.weight});
    }
    return exact;
}

qmomkit::RawMomentTensor exact_cqmom_moments(const std::vector<TerminalNode>& nodes, const std::vector<std::size_t>& extents)
{
    qmomkit::RawMomentTensor moments{extents, std::vector<double>(qmomkit::product_extents(extents), 0.0)};
    for (std::size_t flat = 0; flat < moments.values.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, extents);
        for (const TerminalNode& node : nodes)
        {
            double contribution = node.weight;
            for (std::size_t dimension = 0; dimension < extents.size(); ++dimension)
            {
                contribution *= integer_power(node.coordinates[dimension], index[dimension]);
            }
            moments.values[flat] += contribution;
        }
    }
    return moments;
}

std::vector<TerminalNode> sort_nodes(std::vector<TerminalNode> nodes)
{
    std::sort(nodes.begin(), nodes.end(),
              [](const TerminalNode& lhs, const TerminalNode& rhs) { return coordinate_less(lhs.coordinates, rhs.coordinates); });
    return nodes;
}

void compare_cqmom_state(const qmomkit::CqmomQuadrature& quadrature, const std::vector<TerminalNode>& exact_nodes,
                         const std::vector<std::size_t>& extents, VerificationResult& result)
{
    std::vector<TerminalNode> actual = sort_nodes(terminal_nodes(quadrature));
    std::vector<TerminalNode> expected = sort_nodes(exact_nodes);
    if (actual.size() != expected.size())
    {
        append_failure(result, "CQMOM terminal node count changed.");
        return;
    }
    for (std::size_t node = 0; node < expected.size(); ++node)
    {
        for (std::size_t dimension = 0; dimension < expected[node].coordinates.size(); ++dimension)
        {
            update_error(result.errors.node, actual[node].coordinates[dimension], expected[node].coordinates[dimension]);
        }
        update_error(result.errors.weight, actual[node].weight, expected[node].weight);
        if (!(actual[node].weight > 0.0))
        {
            append_failure(result, "CQMOM terminal weight is not positive.");
        }
    }

    const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(quadrature, extents);
    const qmomkit::RawMomentTensor exact = exact_cqmom_moments(exact_nodes, extents);
    for (std::size_t index = 0; index < exact.values.size(); ++index)
    {
        update_error(result.errors.moment_absolute, reconstructed.values[index], exact.values[index]);
        update_relative_error(result.errors.moment_relative, reconstructed.values[index], exact.values[index]);
    }
    const std::size_t m00 = qmomkit::flatten_multi_index({0, 0}, extents);
    if (std::abs(reconstructed.values[m00] - exact.values[m00]) > tolerance)
    {
        append_failure(result, "CQMOM M[0,0] was not preserved.");
    }
    for (const qmomkit::MultiIndex& cross_index : {qmomkit::MultiIndex{1, 1}, qmomkit::MultiIndex{2, 1}, qmomkit::MultiIndex{1, 2}})
    {
        const std::size_t flat = qmomkit::flatten_multi_index(cross_index, extents);
        if (std::abs(reconstructed.values[flat] - exact.values[flat]) > tolerance)
        {
            append_failure(result, "CQMOM selected cross-moment check failed.");
        }
    }
}

void verify_initial_cqmom_state(const qmomkit::CqmomCase& cqmom_case, const qmomkit::CqmomQuadrature& quadrature, VerificationResult& result)
{
    const std::vector<TerminalNode> expected_nodes = expected_initial_linear_system_nodes();
    const std::vector<std::size_t> extents = cqmom_case.moments.layout.hierarchy.moment_extents();
    compare_cqmom_state(quadrature, expected_nodes, extents, result);

    // The case holds the moments its hierarchy reconstructs from, so the
    // analytically summed tensor is reduced the same way before comparison.
    const qmomkit::CqmomMomentData expected_moments =
        qmomkit::extract_cqmom_moments(exact_cqmom_moments(expected_nodes, extents), cqmom_case.moments.layout);
    for (std::size_t index = 0; index < expected_moments.values.size(); ++index)
    {
        update_error(result.errors.moment_absolute, cqmom_case.moments.values[index], expected_moments.values[index]);
        update_relative_error(result.errors.moment_relative, cqmom_case.moments.values[index], expected_moments.values[index]);
    }
}

void verify_cqmom_plan(const qmomkit::CqmomQuadrature& initial, const TimePlan& plan, const std::vector<std::size_t>& extents,
                       VerificationResult& result)
{
    qmomkit::CqmomQuadrature current = initial;
    const std::vector<TerminalNode> initial_nodes = terminal_nodes(initial);
    double time = plan.t0;
    for (std::size_t step = 0; step < plan.steps; ++step)
    {
        const qmomkit::CqmomNodeAdvectionResult advanced =
            qmomkit::advance_nodes(current, time, plan.dt,
                                   qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi)
                                                                      { return std::vector<double>{-0.3 * xi[0], 0.8 * xi[0] + 0.2 * xi[1]}; }),
                                   accurate_options());
        if (!advanced.ok())
        {
            append_failure(result, "CQMOM node advection failed: " + advanced.message);
            return;
        }
        if (advanced.diagnostics.adaptive_status != qmomkit::AdaptiveRkStatus::success || advanced.diagnostics.step_count == 0)
        {
            append_failure(result, "CQMOM RK23 diagnostics did not report successful integration.");
        }
        update_error(result.errors.node, advanced.diagnostics.final_time, time + plan.dt);
        const std::vector<TerminalNode> current_nodes = terminal_nodes(advanced.quadrature);
        for (std::size_t node = 0; node < current_nodes.size(); ++node)
        {
            update_error(result.errors.weight, current_nodes[node].weight, initial_nodes[node].weight);
        }
        current = advanced.quadrature;
        time += plan.dt;
    }
    compare_cqmom_state(current, exact_linear_system_nodes(initial_nodes, time - plan.t0), extents, result);
}

VerificationResult verify_cqmom_case(const CaseDefinition& definition, const std::filesystem::path& configuration_path)
{
    VerificationResult result;
    const auto configured = qmomkit::config::load_cqmom(configuration_path);
    if (!configured.ok())
    {
        append_failure(result, "Failed to load CQMOM configuration.");
        return result;
    }
    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(configured.value());
    if (!cqmom.ok())
    {
        append_failure(result, "Initial CQMOM inversion failed: " + cqmom.message);
        return result;
    }
    if (cqmom.quadrature.number_of_dimensions() != 2 || cqmom.quadrature.hierarchy().conditioning_order() != std::vector<std::size_t>{1, 0})
    {
        append_failure(result, "2D linear-system case did not use the expected non-identity conditioning order.");
    }
    if (cqmom.quadrature.has_rectangular_storage() || !cqmom.quadrature.adaptive_branch_counts_enabled() || cqmom.quadrature.weights().size() != 3)
    {
        append_failure(result, "2D linear-system case did not produce the expected adaptive ragged CQMOM quadrature.");
    }
    verify_initial_cqmom_state(configured.value(), cqmom.quadrature, result);
    verify_cqmom_plan(cqmom.quadrature, definition.one_step, configured.value().moments.layout.hierarchy.moment_extents(), result);
    verify_cqmom_plan(cqmom.quadrature, definition.multi_step, configured.value().moments.layout.hierarchy.moment_extents(), result);
    return result;
}

} // namespace

int main(const int argc, char** argv)
{
    if (argc != 3)
    {
        std::cerr << "Usage: node_advection_analytical <configuration.qmom> <case_identifier>\n";
        return 2;
    }
    const std::filesystem::path configuration_path = argv[1];
    const std::string case_identifier = argv[2];
    const CaseDefinition* const definition = find_case(case_identifier);
    if (definition == nullptr)
    {
        std::cerr << "Unknown node-advection analytical case '" << case_identifier << "'.\n";
        return 2;
    }

    const VerificationResult result =
        definition->cqmom ? verify_cqmom_case(*definition, configuration_path) : verify_qmom_case(*definition, configuration_path);
    const bool passed = result.passed && result.errors.node <= tolerance && result.errors.weight <= tolerance &&
                        result.errors.moment_relative <= tolerance && result.errors.moment_absolute <= 1.0e-5;

    // Silent on success. These numbers exist to diagnose a failure, and CTest
    // surfaces stderr through --output-on-failure.
    if (!(passed))
    {
        std::cerr << "QMOMKit node-advection analytical verification\n";
        std::cerr << "Configuration: " << configuration_path << '\n';
        std::cerr << "Tolerance: " << std::scientific << tolerance << " (nodes/weights and unit-floored relative moments)\n\n";
        std::cerr << std::left << std::setw(28) << "case" << std::right << std::setw(14) << "node abs" << std::setw(14) << "weight abs"
                  << std::setw(14) << "moment abs" << std::setw(14) << "moment rel" << std::setw(14) << "loose err" << std::setw(14) << "tight err"
                  << std::setw(9) << "status" << '\n';
        std::cerr << std::string(121, '-') << '\n';
        std::cerr << std::left << std::setw(28) << definition->name << std::right << std::scientific << std::setprecision(4) << std::setw(14)
                  << result.errors.node << std::setw(14) << result.errors.weight << std::setw(14) << result.errors.moment_absolute << std::setw(14)
                  << result.errors.moment_relative << std::setw(14) << result.loose_error << std::setw(14) << result.tight_error << std::setw(9)
                  << (passed ? "PASS" : "FAIL") << '\n';

        if (!result.passed)
        {
            std::cerr << result.message << '\n';
        }
    }

    return passed ? 0 : 1;
}
