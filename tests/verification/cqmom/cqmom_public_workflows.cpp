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

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::CqmomDenseCase>& configured)
{
    std::cerr << configured.error() << '\n';
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
            nodes.push_back(
                {quadrature.terminal_paths()[node], quadrature.abscissas()[0][node], quadrature.abscissas()[1][node], quadrature.weights()[node]});
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

void verify_lean_required_representation(const qmomkit::CqmomLeanCase& lean_case, const qmomkit::CqmomDenseCase& dense_case,
                                         VerificationResult& verification)
{
    // Expert callers choose CqmomLeanCase when they already work with the
    // minimal recursive moment subset. This check verifies the public loader
    // provides exactly that selected representation and that the selected
    // values are sufficient to recover the same analytical product quadrature.
    const std::vector<qmomkit::MultiIndex> expected_indices = expected_required_indices();
    const qmomkit::CqmomMomentData& moments = lean_case.moments;
    if (moments.layout.required_indices.to_multi_indices() != expected_indices)
    {
        append_failure(verification, "explicit lean workflow required-index order did not match the documented fixed 2x2 CQMOM layout.");
    }
    if (moments.values.size() != expected_indices.size())
    {
        append_failure(verification, "explicit lean workflow did not provide one value per required index.");
    }
    if (moments.values.size() >= dense_case.moments.values.size())
    {
        append_failure(verification, "explicit lean workflow did not use a smaller selected moment representation than the dense tensor.");
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

std::vector<double> reconstructed_values(const qmomkit::CqmomQuadrature& quadrature, const std::vector<std::size_t>& extents)
{
    return qmomkit::reconstruct_raw_moments(quadrature, extents).values;
}

void verify_reconstructed_moments(const std::string& workflow_name, const qmomkit::CqmomQuadrature& quadrature,
                                  const qmomkit::RawMomentTensor& dense_moments, VerificationResult& verification)
{
    // The moment tensor in independent_product.qmom is generated from the same
    // independent analytical quadrature documented above. Matching every
    // configured raw moment verifies that the public workflow returned a
    // mathematically equivalent CQMOM rule, not only visually similar nodes.
    const std::vector<double> reconstructed = reconstructed_values(quadrature, dense_moments.extents);
    if (reconstructed.size() != dense_moments.values.size())
    {
        append_failure(verification, workflow_name + " reconstructed moment tensor has the wrong size.");
        return;
    }
    for (std::size_t index = 0; index < dense_moments.values.size(); ++index)
    {
        update_error(verification.errors.moment, reconstructed[index], dense_moments.values[index]);
    }
}

void compare_workflows(const std::string& workflow_name, const qmomkit::CqmomQuadrature& actual, const qmomkit::CqmomQuadrature& dense_reference,
                       const std::vector<std::size_t>& extents, VerificationResult& verification)
{
    // Dense, default, and lean workflows exist for different user needs:
    // explicit full tensors, recommended configured state ownership, and
    // explicit minimal required-moment storage. For the same mathematical
    // problem they should produce indistinguishable public quadrature and
    // reconstructed moments.
    const std::vector<TerminalNode> actual_nodes = terminal_nodes(actual);
    const std::vector<TerminalNode> reference_nodes = terminal_nodes(dense_reference);
    if (actual_nodes.size() != reference_nodes.size())
    {
        append_failure(verification, workflow_name + " terminal node count differs from dense workflow.");
        return;
    }
    for (std::size_t node = 0; node < actual_nodes.size(); ++node)
    {
        if (actual_nodes[node].path != reference_nodes[node].path)
        {
            append_failure(verification, workflow_name + " terminal path differs from dense workflow.");
            continue;
        }
        update_error(verification.errors.workflow_node, actual_nodes[node].x, reference_nodes[node].x);
        update_error(verification.errors.workflow_node, actual_nodes[node].y, reference_nodes[node].y);
        update_error(verification.errors.workflow_weight, actual_nodes[node].weight, reference_nodes[node].weight);
    }

    const std::vector<double> actual_moments = reconstructed_values(actual, extents);
    const std::vector<double> reference_moments = reconstructed_values(dense_reference, extents);
    if (actual_moments.size() != reference_moments.size())
    {
        append_failure(verification, workflow_name + " reconstructed moment count differs from dense workflow.");
        return;
    }
    for (std::size_t index = 0; index < actual_moments.size(); ++index)
    {
        update_error(verification.errors.workflow_moment, actual_moments[index], reference_moments[index]);
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

    const auto dense_case = qmomkit::config::load_cqmom_dense(configuration_path);
    if (!dense_case.ok())
    {
        print_configuration_errors(dense_case);
        return 1;
    }
    const auto lean_case = qmomkit::config::load_cqmom_lean(configuration_path);
    if (!lean_case.ok())
    {
        print_configuration_errors(lean_case);
        return 1;
    }
    const auto default_case = qmomkit::config::load_cqmom(configuration_path);
    if (!default_case.ok())
    {
        print_configuration_errors(default_case);
        return 1;
    }

    const qmomkit::CqmomResult dense = qmomkit::compute_cqmom(dense_case.value());
    const qmomkit::CqmomResult lean = qmomkit::compute_cqmom(lean_case.value());
    const qmomkit::CqmomResult default_workflow = qmomkit::compute_cqmom(default_case.value());

    VerificationResult verification;
    if (!dense.ok())
    {
        append_failure(verification, "dense workflow compute_cqmom failed: " + dense.message);
    }
    if (!lean.ok())
    {
        append_failure(verification, "lean workflow compute_cqmom failed: " + lean.message);
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

    if (dense_case.value().hierarchy.node_counts() != std::vector<std::size_t>({2, 2}) ||
        dense_case.value().hierarchy.conditioning_order() != std::vector<std::size_t>({0, 1}))
    {
        append_failure(verification, "dense workflow did not load the expected fixed 2x2 hierarchy with conditioning_order = {0, 1}.");
    }
    if (default_case.value().hierarchy().node_counts() != dense_case.value().hierarchy.node_counts() ||
        default_case.value().hierarchy().conditioning_order() != dense_case.value().hierarchy.conditioning_order())
    {
        append_failure(verification, "default CqmomCase workflow did not preserve the dense workflow hierarchy.");
    }
    if (lean_case.value().moments.layout.hierarchy.node_counts() != dense_case.value().hierarchy.node_counts() ||
        lean_case.value().moments.layout.hierarchy.conditioning_order() != dense_case.value().hierarchy.conditioning_order())
    {
        append_failure(verification, "explicit lean workflow did not preserve the dense workflow hierarchy.");
    }

    verify_lean_required_representation(lean_case.value(), dense_case.value(), verification);
    verify_diagnostics("dense workflow", dense, verification);
    verify_diagnostics("explicit lean workflow", lean, verification);
    verify_diagnostics("default CqmomCase workflow", default_workflow, verification);
    verify_reference_quadrature("dense workflow", dense.quadrature, verification);
    verify_reference_quadrature("explicit lean workflow", lean.quadrature, verification);
    verify_reference_quadrature("default CqmomCase workflow", default_workflow.quadrature, verification);
    verify_reconstructed_moments("dense workflow", dense.quadrature, dense_case.value().moments, verification);
    verify_reconstructed_moments("explicit lean workflow", lean.quadrature, dense_case.value().moments, verification);
    verify_reconstructed_moments("default CqmomCase workflow", default_workflow.quadrature, dense_case.value().moments, verification);
    compare_workflows("explicit lean workflow", lean.quadrature, dense.quadrature, dense_case.value().moments.extents, verification);
    compare_workflows("default CqmomCase workflow", default_workflow.quadrature, dense.quadrature, dense_case.value().moments.extents, verification);

    const bool passed = verification.passed && verification.errors.node <= tolerance && verification.errors.weight <= tolerance &&
                        verification.errors.moment <= tolerance && verification.errors.workflow_node <= tolerance &&
                        verification.errors.workflow_weight <= tolerance && verification.errors.workflow_moment <= tolerance;

    std::cout << "QMOMKit CQMOM public workflow verification\n";
    std::cout << "Configuration: " << configuration_path << '\n';
    std::cout << "Reference: independent 2D product quadrature with known terminal weights and coordinates\n";
    std::cout << "Tolerance: " << std::scientific << tolerance << " absolute\n\n";
    std::cout << std::left << std::setw(28) << "check" << std::right << std::setw(14) << "max error" << '\n';
    std::cout << std::string(42, '-') << '\n';
    std::cout << std::left << std::setw(28) << "reference nodes" << std::right << std::scientific << std::setprecision(4) << std::setw(14)
              << verification.errors.node << '\n';
    std::cout << std::left << std::setw(28) << "reference weights" << std::right << std::setw(14) << verification.errors.weight << '\n';
    std::cout << std::left << std::setw(28) << "reference moments" << std::right << std::setw(14) << verification.errors.moment << '\n';
    std::cout << std::left << std::setw(28) << "workflow nodes" << std::right << std::setw(14) << verification.errors.workflow_node << '\n';
    std::cout << std::left << std::setw(28) << "workflow weights" << std::right << std::setw(14) << verification.errors.workflow_weight << '\n';
    std::cout << std::left << std::setw(28) << "workflow moments" << std::right << std::setw(14) << verification.errors.workflow_moment << '\n';
    std::cout << "\nStatus: " << (passed ? "PASS" : "FAIL") << '\n';

    if (!verification.passed)
    {
        std::cerr << verification.message << '\n';
    }
    return passed ? 0 : 1;
}
