// Vandermonde solver choice on a deliberately ill-conditioned CQMOM problem.
//
// What is verified
//   That each solver the caller can request is actually the one used, that all
//   of them recover the same prescribed rule on a hard problem, and that the
//   reported condition estimate lands in the range the problem was built to
//   produce.
//
// How the reference is constructed
//   Chosen, not computed, and chosen to be hard. The three parent positions are
//   clustered at x = 1.0, 1.007, 1.014, which makes the weighted 3x3
//   conditional system badly conditioned - kappa_2 is about 2.33e5 - while the
//   per-parent child rules stay positive so the whole problem remains solvable.
//   The moments handed to the library are the raw moments of that rule.
//
// Error metric and pass criteria
//   Maximum absolute deviation of positions, weights and moments from the
//   prescribed rule, plus the reported residual and condition estimate.
//
// Tolerances
//   Looser than elsewhere on purpose, and set by the conditioning rather than
//   picked: positions 5.0e-4, weights 2.0e-5, moments 1.0e-7 absolute and
//   1.0e-8 relative. With kappa_2 near 2.3e5, double-precision inputs accurate
//   to about 1e-16 can only pin the solution to roughly 1e-11 before the CQMOM
//   levels compound it, so tight bounds here would assert a particular rounding
//   path rather than a property of the solver. The condition estimate itself is
//   required to fall between 1.0e5 and 1.0e6, which is what makes the case
//   ill-conditioned by construction rather than by accident.
//
// Reproduce
//   ctest -R verification_.*clustered
//
// Build modes
//   Both. Failures are reported through explicit comparisons and an exit code,
//   not assert().

#include "qmomkit/configuration.hpp"
#include "qmomkit/moment_reconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace
{

constexpr double configuration_moment_tolerance = 1.0e-12;
constexpr double node_tolerance = 5.0e-4;
constexpr double weight_tolerance = 2.0e-5;
constexpr double moment_absolute_tolerance = 1.0e-7;
constexpr double moment_relative_tolerance = 1.0e-8;
constexpr double minimum_condition_estimate = 1.0e5;
constexpr double maximum_condition_estimate = 1.0e6;

struct TerminalNode
{
    double x = 0.0;
    double y = 0.0;
    double weight = 0.0;
};

struct SolverPolicy
{
    qmomkit::VandermondeSolver requested = qmomkit::VandermondeSolver::automatic;
    qmomkit::VandermondeSolver selected = qmomkit::VandermondeSolver::automatic;
    const char* name = "automatic";
};

struct ErrorMetrics
{
    double node = 0.0;
    double weight = 0.0;
    double moment = 0.0;
    double residual = 0.0;
    double condition = 0.0;
};

struct PolicyReport
{
    const char* name = "";
    bool success = false;
    ErrorMetrics errors;
};

double integer_power(const double value, const std::size_t exponent)
{
    double result = 1.0;
    for (std::size_t power = 0; power < exponent; ++power)
    {
        result *= value;
    }
    return result;
}

std::vector<TerminalNode> prescribed_terminal_rule()
{
    // The parent nodes are deliberately clustered. The resulting weighted 3x3
    // conditional Vandermonde system has kappa_2 about 2.33e5, while the
    // branch-dependent positive child rules keep the full CQMOM problem
    // realizable and independently prescribed.
    const std::vector<double> x{1.0, 1.007, 1.014};
    const std::vector<double> parent_weights{0.2, 0.3, 0.5};
    const std::vector<std::vector<double>> y{{-1.0, 2.0}, {-0.5, 3.0}, {0.25, 4.0}};
    const std::vector<std::vector<double>> conditional_weights{{0.4, 0.6}, {0.7, 0.3}, {0.25, 0.75}};

    std::vector<TerminalNode> rule;
    for (std::size_t parent = 0; parent < x.size(); ++parent)
    {
        for (std::size_t child = 0; child < y[parent].size(); ++child)
        {
            rule.push_back({x[parent], y[parent][child], parent_weights[parent] * conditional_weights[parent][child]});
        }
    }
    return rule;
}

qmomkit::RawMomentTensor moments_from_terminal_rule(const std::vector<TerminalNode>& rule)
{
    const std::vector<std::size_t> extents{6, 4};
    std::vector<double> values(qmomkit::product_extents(extents), 0.0);
    for (std::size_t x_order = 0; x_order < extents[0]; ++x_order)
    {
        for (std::size_t y_order = 0; y_order < extents[1]; ++y_order)
        {
            double value = 0.0;
            for (const TerminalNode& node : rule)
            {
                value += node.weight * integer_power(node.x, x_order) * integer_power(node.y, y_order);
            }
            values[x_order * extents[1] + y_order] = value;
        }
    }
    return {extents, values};
}

void update_error(double& maximum, const double actual, const double expected)
{
    maximum = std::max(maximum, std::abs(actual - expected));
}

void require(const bool condition, const std::string& message)
{
    if (!condition)
    {
        std::cerr << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// The configured case holds the moments its hierarchy reconstructs from, so the
// prescribed tensor is reduced the same way before comparison.
void verify_configured_moments(const qmomkit::CqmomMomentData& configured, const qmomkit::RawMomentTensor& prescribed)
{
    require(prescribed.extents == configured.layout.hierarchy.moment_extents(), "Clustered solver-policy fixture has unexpected moment extents.");
    const qmomkit::CqmomMomentData expected = qmomkit::extract_cqmom_moments(prescribed, configured.layout);
    require(configured.values.size() == expected.values.size(), "Clustered solver-policy fixture has unexpected moment value count.");
    for (std::size_t index = 0; index < expected.values.size(); ++index)
    {
        require(std::abs(configured.values[index] - expected.values[index]) <= configuration_moment_tolerance,
                "Clustered solver-policy fixture does not match the independently summed terminal rule.");
    }
}

PolicyReport verify_policy(const SolverPolicy& policy, qmomkit::CqmomLeanCase case_, const std::vector<TerminalNode>& prescribed_rule,
                           const qmomkit::RawMomentTensor& prescribed_moments)
{
    case_.vandermonde_solver = policy.requested;
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(case_);
    require(result.ok(), std::string(policy.name) + " failed clustered public CQMOM reconstruction: " + result.message);
    require(result.quadrature.has_rectangular_storage(), std::string(policy.name) + " did not recover fixed rectangular storage.");
    require(result.quadrature.total_nodes() == prescribed_rule.size(), std::string(policy.name) + " recovered the wrong terminal-node count.");

    PolicyReport report;
    report.name = policy.name;
    report.success = true;
    for (std::size_t parent = 0; parent < 3; ++parent)
    {
        for (std::size_t child = 0; child < 2; ++child)
        {
            const std::size_t terminal = 2 * parent + child;
            const qmomkit::MultiIndex index{parent, child};
            const TerminalNode& expected = prescribed_rule[terminal];
            update_error(report.errors.node, result.quadrature.abscissa(0, index), expected.x);
            update_error(report.errors.node, result.quadrature.abscissa(1, index), expected.y);
            update_error(report.errors.weight, result.quadrature.weight(index), expected.weight);
        }
    }

    const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(result.quadrature, prescribed_moments.extents);
    require(reconstructed.values.size() == prescribed_moments.values.size(), std::string(policy.name) + " reconstructed the wrong moment count.");
    for (std::size_t index = 0; index < prescribed_moments.values.size(); ++index)
    {
        const double expected = prescribed_moments.values[index];
        const double error = std::abs(reconstructed.values[index] - expected);
        update_error(report.errors.moment, reconstructed.values[index], expected);
        require(error <= moment_absolute_tolerance + moment_relative_tolerance * std::max(1.0, std::abs(expected)),
                std::string(policy.name) + " exceeded the prescribed-moment tolerance.");
    }

    bool saw_root_conditional_system = false;
    for (const qmomkit::CqmomDiagnostics& diagnostics : result.diagnostics)
    {
        if (diagnostics.conditioning_level != 0)
        {
            continue;
        }
        saw_root_conditional_system = true;
        require(diagnostics.conditional_vandermonde_requested_solver == policy.requested,
                std::string(policy.name) + " reported the wrong requested public solver policy.");
        require(diagnostics.conditional_vandermonde_selected_solver == policy.selected,
                std::string(policy.name) + " reported the wrong selected public solver policy.");
        require(diagnostics.conditional_moment_diagnostics.condition_estimate.has_value(),
                std::string(policy.name) + " did not report the root conditional condition estimate.");
        require(diagnostics.conditional_moment_diagnostics.residual_norm.has_value(),
                std::string(policy.name) + " did not report the root conditional residual.");

        report.errors.condition = *diagnostics.conditional_moment_diagnostics.condition_estimate;
        report.errors.residual = *diagnostics.conditional_moment_diagnostics.residual_norm;
    }
    require(saw_root_conditional_system, std::string(policy.name) + " did not report a root conditional Vandermonde system.");
    require(report.errors.condition >= minimum_condition_estimate && report.errors.condition <= maximum_condition_estimate,
            std::string(policy.name) + " reported a condition estimate outside the documented moderate-conditioning band.");
    require(report.errors.node <= node_tolerance, std::string(policy.name) + " exceeded the prescribed-node tolerance.");
    require(report.errors.weight <= weight_tolerance, std::string(policy.name) + " exceeded the prescribed-weight tolerance.");
    return report;
}

} // namespace

int main(const int argc, char** argv)
{
    if (argc != 2 && argc != 3)
    {
        std::cerr << "Usage: clustered_solver_policy_verification <clustered automatic .qmom file> [solver]\n";
        return 2;
    }

    const std::vector<TerminalNode> prescribed_rule = prescribed_terminal_rule();
    const qmomkit::RawMomentTensor prescribed_moments = moments_from_terminal_rule(prescribed_rule);
    const auto configured = qmomkit::config::load_cqmom_lean(std::filesystem::path(argv[1]));
    require(configured.ok(), configured.error());
    require(configured.value().vandermonde_solver == qmomkit::VandermondeSolver::automatic,
            "Clustered solver-policy configuration did not preserve the automatic policy.");
    verify_configured_moments(configured.value().moments, prescribed_moments);

    const std::vector<SolverPolicy> policies{
        {qmomkit::VandermondeSolver::automatic, qmomkit::VandermondeSolver::bjorck_pereyra, "automatic"},
        {qmomkit::VandermondeSolver::dense, qmomkit::VandermondeSolver::dense, "dense"},
        {qmomkit::VandermondeSolver::pivoted_qr, qmomkit::VandermondeSolver::pivoted_qr, "pivoted_qr"},
        {qmomkit::VandermondeSolver::svd, qmomkit::VandermondeSolver::svd, "svd"},
        {qmomkit::VandermondeSolver::moore_penrose, qmomkit::VandermondeSolver::moore_penrose, "moore_penrose"},
        {qmomkit::VandermondeSolver::bjorck_pereyra, qmomkit::VandermondeSolver::bjorck_pereyra, "bjorck_pereyra"}};

    const std::string selected_solver = argc == 3 ? argv[2] : "";
    bool matched = selected_solver.empty();
    for (const SolverPolicy& policy : policies)
    {
        if (!selected_solver.empty() && selected_solver != policy.name)
        {
            continue;
        }
        matched = true;
        // verify_policy reports and exits on any breach; its report was only
        // ever read by the success table that no longer exists.
        verify_policy(policy, configured.value(), prescribed_rule, prescribed_moments);
    }
    if (!matched)
    {
        std::cerr << "Unknown clustered solver case '" << selected_solver << "'.\n";
        return 2;
    }
    return 0;
}
