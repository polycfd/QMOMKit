// The two public CQMOM entry points agreeing with each other and with a known
// answer.
//
// What is verified
//   That the dense and lean representations of the same problem produce the
//   same quadrature, that the lean case carries exactly the moments its
//   hierarchy requires and no others, that the diagnostics describe the same
//   reconstruction in both, and that rebuilt moments return the configured
//   values. The point is that choosing a storage form is a data-transfer
//   decision, not a numerical one.
//
// How the reference is constructed
//   Chosen, not computed. A rule of weighted points is written down and its raw
//   moments are evaluated from it; both workflows are then measured against
//   that same independent answer rather than against each other. Comparing the
//   two workflows alone would pass even if both were wrong in the same way.
//
// Error metric and pass criteria
//   Maximum absolute deviation of positions, weights and rebuilt moments from
//   the prescribed rule, at 1.0e-12, applied separately to each workflow. The
//   lean representation is additionally required to hold exactly the required
//   index set, which is a structural check rather than a numerical one.
//
// Reproduce
//   ctest -R verification_cqmom_public_workflows
//
// Build modes
//   Both. Failures are reported through an exit code, not assert().

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

constexpr double tolerance = 1.0e-12;

struct ErrorMetrics
{
    double node = 0.0;
    double weight = 0.0;
    double moment = 0.0;
    double workflow_node = 0.0;
    double workflow_weight = 0.0;
    double workflow_moment = 0.0;
};

struct TerminalNode
{
    qmomkit::ConditioningPath path;
    double x = 0.0;
    double y = 0.0;
    double weight = 0.0;
};

struct VerificationResult
{
    bool passed = true;
    std::string message;
    ErrorMetrics errors;
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

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::CqmomLeanCase>& configured)
{
    std::cerr << configured.error() << '\n';
}

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::CqmomCase>& configured)
{
    std::cerr << configured.error() << '\n';
}

std::vector<TerminalNode> terminal_nodes(const qmomkit::CqmomQuadrature& quadrature)
{
    std::vector<TerminalNode> nodes;
    if (quadrature.has_rectangular_storage())
    {
        for (const qmomkit::MultiIndex& index : qmomkit::enumerate_quadrature_node_multi_indices(quadrature.hierarchy()))
        {
            nodes.push_back(
                {qmomkit::ConditioningPath(index), quadrature.abscissa(0, index), quadrature.abscissa(1, index), quadrature.weight(index)});
        }
    }
    else
    {
        for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
        {
            nodes.push_back({quadrature.terminal_paths()[node].to_conditioning_path(), quadrature.abscissas()[0][node],
                             quadrature.abscissas()[1][node], quadrature.weights()[node]});
        }
    }
    std::sort(nodes.begin(), nodes.end(), [](const TerminalNode& lhs, const TerminalNode& rhs) { return lhs.path < rhs.path; });
    return nodes;
}

std::vector<TerminalNode> independent_product_reference()
{
    // This reference is the independent analytical quadrature behind
    // independent_product.qmom: x has nodes {1, 3} with weights {0.4, 0.6},
    // y has nodes {-1, 2} with local weights {0.25, 0.75}, and independence
    // makes every terminal weight the product of the marginal weights.
    return {{qmomkit::ConditioningPath({0, 0}), 1.0, -1.0, 0.4 * 0.25},
            {qmomkit::ConditioningPath({0, 1}), 1.0, 2.0, 0.4 * 0.75},
            {qmomkit::ConditioningPath({1, 0}), 3.0, -1.0, 0.6 * 0.25},
            {qmomkit::ConditioningPath({1, 1}), 3.0, 2.0, 0.6 * 0.75}};
}

std::vector<qmomkit::MultiIndex> expected_required_indices()
{
    // A fixed 2x2 CQMOM hierarchy with conditioning_order = {0, 1} requires
    // only the recursive moment subset below. The explicit lean workflow stores
    // this selected representation directly instead of requiring a full 4x4
    // dense tensor from the caller.
    return {{0, 0}, {1, 0}, {2, 0}, {3, 0}, {0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}, {1, 3}};
}

double analytical_product_moment(const qmomkit::MultiIndex& index)
{
    double moment = 0.0;
    for (const TerminalNode& terminal : independent_product_reference())
    {
        moment += terminal.weight * std::pow(terminal.x, static_cast<double>(index[0])) * std::pow(terminal.y, static_cast<double>(index[1]));
    }
    return moment;
}

void verify_diagnostics(const std::string& workflow_name, const qmomkit::CqmomResult& result, VerificationResult& verification)
{
    // For this fixed 2x2 workflow, diagnostics should describe one root
    // inversion and one conditional inversion for each x parent branch. A
    // failure here usually means a public workflow did not preserve the same
    // reconstruction policy as the dense analytical CQMOM path.
    if (result.diagnostics.size() != 3)
    {
        append_failure(verification, workflow_name + " diagnostics did not contain the expected three branch records.");
    }

    bool saw_root = false;
    std::map<std::size_t, bool> saw_child_branch;
    for (const qmomkit::CqmomDiagnostics& diagnostics : result.diagnostics)
    {
        if (diagnostics.qmom_status != qmomkit::QmomStatus::success ||
            diagnostics.conditional_moment_status != qmomkit::ConditionalMomentStatus::success)
        {
            append_failure(verification, workflow_name + " diagnostics contain a non-success status.");
        }
        if (diagnostics.requested_node_count != 2 || diagnostics.actual_node_count != 2)
        {
            append_failure(verification, workflow_name + " diagnostics did not report fixed two-node branches.");
        }
        if (diagnostics.adaptive_branch_counts_enabled || diagnostics.adaptive_node_count_reduced)
        {
            append_failure(verification, workflow_name + " diagnostics unexpectedly reported adaptive branch behavior.");
        }
        if (diagnostics.conditioning_level == 0 && diagnostics.parent_path.empty() && diagnostics.phase_space_dimension == 0)
        {
            saw_root = true;
        }
        if (diagnostics.conditioning_level == 1 && diagnostics.parent_path.size() == 1 && diagnostics.phase_space_dimension == 1)
        {
            saw_child_branch[diagnostics.parent_path.node_index_at_level(0)] = true;
        }
    }
    if (!saw_root || !saw_child_branch[0] || !saw_child_branch[1])
    {
        append_failure(verification, workflow_name + " diagnostics did not report the expected root and child branches.");
    }
}

void verify_required_representation(const qmomkit::CqmomLeanCase& case_, VerificationResult& verification)
{
    // Expert callers choose CqmomLeanCase when they already work with the
    // minimal recursive moment subset. This check verifies the public loader
    // provides exactly that selected representation and that the selected
    // values are sufficient to recover the same analytical product quadrature.
    const std::vector<qmomkit::MultiIndex> expected_indices = expected_required_indices();
    const qmomkit::CqmomMomentData& moments = case_.moments;
    if (moments.layout.required_indices.to_multi_indices() != expected_indices)
    {
        append_failure(verification, "explicit lean workflow required-index order did not match the documented fixed 2x2 CQMOM layout.");
    }
    if (moments.values.size() != expected_indices.size())
    {
        append_failure(verification, "explicit lean workflow did not provide one value per required index.");
    }
    if (moments.values.size() >= qmomkit::product_extents(moments.layout.hierarchy.moment_extents()))
    {
        append_failure(verification, "the required moment set was not smaller than the full moment tensor the hierarchy spans.");
    }

    const std::size_t entries = std::min(moments.values.size(), expected_indices.size());
    for (std::size_t entry = 0; entry < entries; ++entry)
    {
        const qmomkit::MultiIndex& index = expected_indices[entry];
        const double expected = analytical_product_moment(index);
        update_error(verification.errors.moment, moments.values[entry], expected);
        if (!moments.contains_required_index(index))
        {
            append_failure(verification, "explicit lean workflow does not contain an expected required moment index.");
            continue;
        }
        update_error(verification.errors.moment, moments.moment(index), expected);
    }
}

void verify_reference_quadrature(const std::string& workflow_name, const qmomkit::CqmomQuadrature& quadrature, VerificationResult& verification)
{
    const std::vector<TerminalNode> expected = independent_product_reference();
    const std::vector<TerminalNode> actual = terminal_nodes(quadrature);
    if (actual.size() != expected.size())
    {
        append_failure(verification, workflow_name + " terminal node count differs from the analytical reference.");
        return;
    }

    for (std::size_t node = 0; node < expected.size(); ++node)
    {
        if (actual[node].path != expected[node].path)
        {
            append_failure(verification, workflow_name + " terminal path ordering differs from the analytical reference.");
            continue;
        }
        update_error(verification.errors.node, actual[node].x, expected[node].x);
        update_error(verification.errors.node, actual[node].y, expected[node].y);
        update_error(verification.errors.weight, actual[node].weight, expected[node].weight);
    }
}

void verify_reconstructed_moments(const std::string& workflow_name, const qmomkit::CqmomQuadrature& quadrature,
                                  const qmomkit::CqmomMomentData& configured_moments, VerificationResult& verification)
{
    // The moment tensor in independent_product.qmom is generated from the same
    // independent analytical quadrature documented above. Matching every
    // configured raw moment verifies that the public workflow returned a
    // mathematically equivalent CQMOM rule, not only visually similar nodes.
    const qmomkit::CqmomMomentData reconstructed = qmomkit::reconstruct_cqmom_moments(quadrature, configured_moments.layout);
    if (reconstructed.values.size() != configured_moments.values.size())
    {
        append_failure(verification, workflow_name + " reconstructed moment set has the wrong size.");
        return;
    }
    for (std::size_t index = 0; index < configured_moments.values.size(); ++index)
    {
        update_error(verification.errors.moment, reconstructed.values[index], configured_moments.values[index]);
    }
}

void compare_workflows(const std::string& workflow_name, const qmomkit::CqmomQuadrature& actual, const qmomkit::CqmomQuadrature& reference,
                       const qmomkit::CqmomMomentLayout& layout, VerificationResult& verification)
{
    // The explicit case and the configured CqmomCase exist for different user
    // needs: an explicit moment set the caller owns, and the recommended
    // configured state ownership. For the same mathematical problem they should
    // produce indistinguishable public quadrature and reconstructed moments.
    const std::vector<TerminalNode> actual_nodes = terminal_nodes(actual);
    const std::vector<TerminalNode> reference_nodes = terminal_nodes(reference);
    if (actual_nodes.size() != reference_nodes.size())
    {
        append_failure(verification, workflow_name + " terminal node count differs from the explicit workflow.");
        return;
    }
    for (std::size_t node = 0; node < actual_nodes.size(); ++node)
    {
        if (actual_nodes[node].path != reference_nodes[node].path)
        {
            append_failure(verification, workflow_name + " terminal path differs from the explicit workflow.");
            continue;
        }
        update_error(verification.errors.workflow_node, actual_nodes[node].x, reference_nodes[node].x);
        update_error(verification.errors.workflow_node, actual_nodes[node].y, reference_nodes[node].y);
        update_error(verification.errors.workflow_weight, actual_nodes[node].weight, reference_nodes[node].weight);
    }

    const qmomkit::CqmomMomentData actual_moments = qmomkit::reconstruct_cqmom_moments(actual, layout);
    const qmomkit::CqmomMomentData reference_moments = qmomkit::reconstruct_cqmom_moments(reference, layout);
    if (actual_moments.values.size() != reference_moments.values.size())
    {
        append_failure(verification, workflow_name + " reconstructed moment count differs from the explicit workflow.");
        return;
    }
    for (std::size_t index = 0; index < actual_moments.values.size(); ++index)
    {
        update_error(verification.errors.workflow_moment, actual_moments.values[index], reference_moments.values[index]);
    }
}

} // namespace

int main(const int argc, char** argv)
{
    // Verify that dense, explicit lean, and default CQMOM workflows produce
    // equivalent quadrature for the same independent product distribution. The
    // dense tensor is convenient as configured input, while the analytical
    // terminal rule and reconstructed moments are the actual reference.
    if (argc != 2)
    {
        std::cerr << "Usage: cqmom_public_workflows <independent_product.qmom>\n";
        return 2;
    }

    const std::filesystem::path configuration_path = argv[1];

    const auto explicit_case = qmomkit::config::load_cqmom_lean(configuration_path);
    if (!explicit_case.ok())
    {
        print_configuration_errors(explicit_case);
        return 1;
    }
    const auto default_case = qmomkit::config::load_cqmom(configuration_path);
    if (!default_case.ok())
    {
        print_configuration_errors(default_case);
        return 1;
    }

    const qmomkit::CqmomResult explicit_result = qmomkit::compute_cqmom(explicit_case.value());
    const qmomkit::CqmomResult default_workflow = qmomkit::compute_cqmom(default_case.value());

    VerificationResult verification;
    if (!explicit_result.ok())
    {
        append_failure(verification, "explicit workflow compute_cqmom failed: " + explicit_result.message);
    }
    if (!default_workflow.ok())
    {
        append_failure(verification, "default CqmomCase workflow compute_cqmom failed: " + default_workflow.message);
    }
    if (!verification.passed)
    {
        std::cerr << verification.message << '\n';
        return 1;
    }

    if (explicit_case.value().moments.layout.hierarchy.node_counts() != std::vector<std::size_t>({2, 2}) ||
        explicit_case.value().moments.layout.hierarchy.conditioning_order() != std::vector<std::size_t>({0, 1}))
    {
        append_failure(verification, "explicit workflow did not load the expected fixed 2x2 hierarchy with conditioning_order = {0, 1}.");
    }
    if (default_case.value().hierarchy().node_counts() != explicit_case.value().moments.layout.hierarchy.node_counts() ||
        default_case.value().hierarchy().conditioning_order() != explicit_case.value().moments.layout.hierarchy.conditioning_order())
    {
        append_failure(verification, "default CqmomCase workflow did not preserve the explicit workflow hierarchy.");
    }

    verify_required_representation(explicit_case.value(), verification);
    verify_diagnostics("explicit workflow", explicit_result, verification);
    verify_diagnostics("default CqmomCase workflow", default_workflow, verification);
    verify_reference_quadrature("explicit workflow", explicit_result.quadrature, verification);
    verify_reference_quadrature("default CqmomCase workflow", default_workflow.quadrature, verification);
    verify_reconstructed_moments("explicit workflow", explicit_result.quadrature, explicit_case.value().moments, verification);
    verify_reconstructed_moments("default CqmomCase workflow", default_workflow.quadrature, explicit_case.value().moments, verification);
    compare_workflows("default CqmomCase workflow", default_workflow.quadrature, explicit_result.quadrature, explicit_case.value().moments.layout,
                      verification);

    const bool passed = verification.passed && verification.errors.node <= tolerance && verification.errors.weight <= tolerance &&
                        verification.errors.moment <= tolerance && verification.errors.workflow_node <= tolerance &&
                        verification.errors.workflow_weight <= tolerance && verification.errors.workflow_moment <= tolerance;

    // Silent on success. These numbers exist to diagnose a failure, and CTest
    // surfaces stderr through --output-on-failure.
    if (!(passed))
    {
        std::cerr << "QMOMKit CQMOM public workflow verification\n";
        std::cerr << "Configuration: " << configuration_path << '\n';
        std::cerr << "Reference: independent 2D product quadrature with known terminal weights and coordinates\n";
        std::cerr << "Tolerance: " << std::scientific << tolerance << " absolute\n\n";
        std::cerr << std::left << std::setw(28) << "check" << std::right << std::setw(14) << "max error" << '\n';
        std::cerr << std::string(42, '-') << '\n';
        std::cerr << std::left << std::setw(28) << "reference nodes" << std::right << std::scientific << std::setprecision(4) << std::setw(14)
                  << verification.errors.node << '\n';
        std::cerr << std::left << std::setw(28) << "reference weights" << std::right << std::setw(14) << verification.errors.weight << '\n';
        std::cerr << std::left << std::setw(28) << "reference moments" << std::right << std::setw(14) << verification.errors.moment << '\n';
        std::cerr << std::left << std::setw(28) << "workflow nodes" << std::right << std::setw(14) << verification.errors.workflow_node << '\n';
        std::cerr << std::left << std::setw(28) << "workflow weights" << std::right << std::setw(14) << verification.errors.workflow_weight << '\n';
        std::cerr << std::left << std::setw(28) << "workflow moments" << std::right << std::setw(14) << verification.errors.workflow_moment << '\n';
        std::cerr << "\nStatus: " << (passed ? "PASS" : "FAIL") << '\n';

        if (!verification.passed)
        {
            std::cerr << verification.message << '\n';
        }
    }

    return passed ? 0 : 1;
}
