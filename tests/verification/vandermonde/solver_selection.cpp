#include "../verification_helpers.hpp"
#include "internal/recursive_quadrature.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/moment_reconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace
{

constexpr double reference_tolerance = 1.0e-12;
constexpr double workflow_tolerance = 1.0e-12;

struct TerminalNode
{
    qmomkit::ConditioningPath path;
    double x = 0.0;
    double y = 0.0;
    double weight = 0.0;
};

struct SolverCase
{
    std::string name;
    std::filesystem::path filename;
    qmomkit::VandermondeSolver requested_solver = qmomkit::VandermondeSolver::automatic;
    qmomkit::VandermondeSolver expected_selected_solver = qmomkit::VandermondeSolver::automatic;
    qmomkit::internal::VandermondeAlgorithm expected_algorithm = qmomkit::internal::VandermondeAlgorithm::dense;
};

struct ErrorMetrics
{
    double reference_node = 0.0;
    double reference_weight = 0.0;
    double reference_moment = 0.0;
    double dense_bp_node = 0.0;
    double dense_bp_weight = 0.0;
    double dense_bp_moment = 0.0;
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

const char* solver_name(const qmomkit::VandermondeSolver solver)
{
    switch (solver)
    {
    case qmomkit::VandermondeSolver::automatic:
        return "automatic";
    case qmomkit::VandermondeSolver::dense:
        return "dense";
    case qmomkit::VandermondeSolver::pivoted_qr:
        return "pivoted_qr";
    case qmomkit::VandermondeSolver::svd:
        return "svd";
    case qmomkit::VandermondeSolver::moore_penrose:
        return "moore_penrose";
    case qmomkit::VandermondeSolver::bjorck_pereyra:
        return "bjorck_pereyra";
    }
    return "unknown";
}

const char* algorithm_name(const qmomkit::internal::VandermondeAlgorithm algorithm)
{
    switch (algorithm)
    {
    case qmomkit::internal::VandermondeAlgorithm::dense:
        return "dense";
    case qmomkit::internal::VandermondeAlgorithm::pivoted_qr:
        return "pivoted_qr";
    case qmomkit::internal::VandermondeAlgorithm::svd:
        return "svd";
    case qmomkit::internal::VandermondeAlgorithm::moore_penrose:
        return "moore_penrose";
    case qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra:
        return "bjorck_pereyra";
    }
    return "unknown";
}

qmomkit::WheelerOptions wheeler_options_from_qmom(const qmomkit::QmomOptions& options)
{
    return {options.adaptive, options.rmin, options.eabs, options.cutoff};
}

std::vector<TerminalNode> analytical_reference()
{
    // The selected Vandermonde system is deliberately small and
    // well-conditioned: a fixed 2x2 independent product CQMOM case. The x
    // marginal has separated nodes {1, 3}; the y conditional rule is identical
    // on every parent branch. This makes the expected terminal quadrature
    // independent of the solver implementation while still requiring a
    // branch-local weighted Vandermonde solve.
    return {{qmomkit::ConditioningPath({0, 0}), 1.0, -1.0, 0.4 * 0.25},
            {qmomkit::ConditioningPath({0, 1}), 1.0, 2.0, 0.4 * 0.75},
            {qmomkit::ConditioningPath({1, 0}), 3.0, -1.0, 0.6 * 0.25},
            {qmomkit::ConditioningPath({1, 1}), 3.0, 2.0, 0.6 * 0.75}};
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

void verify_against_reference(const SolverCase& solver_case, const qmomkit::CqmomResult& result, const qmomkit::RawMomentTensor& moments,
                              VerificationResult& verification)
{
    // Each solver policy is first judged against the independent product
    // quadrature, not against another solver. Reconstructing the raw tensor
    // then verifies the public mathematical contract of the returned rule.
    const std::vector<TerminalNode> actual = terminal_nodes(result.quadrature);
    const std::vector<TerminalNode> expected = analytical_reference();
    if (actual.size() != expected.size())
    {
        append_failure(verification, solver_case.name + " recovered the wrong number of terminal nodes.");
        return;
    }
    for (std::size_t node = 0; node < expected.size(); ++node)
    {
        if (actual[node].path != expected[node].path)
        {
            append_failure(verification, solver_case.name + " recovered an unexpected terminal path.");
            continue;
        }
        update_error(verification.errors.reference_node, actual[node].x, expected[node].x);
        update_error(verification.errors.reference_node, actual[node].y, expected[node].y);
        update_error(verification.errors.reference_weight, actual[node].weight, expected[node].weight);
    }

    const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(result.quadrature, moments.extents);
    if (reconstructed.values.size() != moments.values.size())
    {
        append_failure(verification, solver_case.name + " reconstructed the wrong number of raw moments.");
        return;
    }
    for (std::size_t index = 0; index < moments.values.size(); ++index)
    {
        update_error(verification.errors.reference_moment, reconstructed.values[index], moments.values[index]);
    }
}

void verify_solver_diagnostics(const SolverCase& solver_case, const qmomkit::CqmomDenseCase& dense_case, VerificationResult& verification)
{
    // Public CqmomResult diagnostics currently report CQMOM branch status but
    // not the Vandermonde requested/selected/algorithm fields. The recursive
    // diagnostics are inspected here only for that diagnostic contract. The
    // configuration, computation, quadrature, and moment checks above remain on
    // public user-facing APIs.
    const qmomkit::internal::RecursiveQuadratureRequest request{dense_case.hierarchy,
                                                                dense_case.moments,
                                                                dense_case.conditional_options,
                                                                wheeler_options_from_qmom(dense_case.qmom_options),
                                                                dense_case.adaptive_branch_counts,
                                                                dense_case.diagnostic_level,
                                                                dense_case.vandermonde_solver};
    const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature(request);
    if (!result.ok())
    {
        append_failure(verification, solver_case.name + " internal diagnostic reconstruction failed: " + result.message);
        return;
    }

    bool saw_conditional_solve = false;
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostics : result.diagnostics)
    {
        if (diagnostics.conditioning_level + 1 >= dense_case.hierarchy.number_of_dimensions())
        {
            continue;
        }

        saw_conditional_solve = true;
        if (diagnostics.conditional_vandermonde_requested_solver != solver_case.requested_solver)
        {
            append_failure(verification, solver_case.name + " diagnostics reported the wrong requested solver.");
        }
        if (diagnostics.conditional_vandermonde_selected_solver != solver_case.expected_selected_solver)
        {
            append_failure(verification, solver_case.name + " diagnostics reported the wrong selected solver policy.");
        }
        if (diagnostics.conditional_vandermonde_algorithm_used != solver_case.expected_algorithm)
        {
            append_failure(verification, solver_case.name + " diagnostics reported the wrong algorithm actually used.");
        }
    }
    if (!saw_conditional_solve)
    {
        append_failure(verification, solver_case.name + " did not produce a conditional Vandermonde diagnostic branch.");
    }
}

bool reports_selected_solver(const qmomkit::CqmomResult& result, const qmomkit::PhaseSpaceHierarchy& hierarchy,
                             const qmomkit::VandermondeSolver expected_solver)
{
    bool saw_conditional_solve = false;
    for (const qmomkit::CqmomDiagnostics& diagnostics : result.diagnostics)
    {
        if (diagnostics.conditioning_level + 1 >= hierarchy.number_of_dimensions())
        {
            continue;
        }

        saw_conditional_solve = true;
        if (diagnostics.conditional_vandermonde_selected_solver != expected_solver)
        {
            return false;
        }
    }
    return saw_conditional_solve;
}

void compare_dense_and_bp(const qmomkit::CqmomResult& dense, const qmomkit::CqmomResult& bp, const std::vector<std::size_t>& extents,
                          VerificationResult& verification)
{
    // This is an additional consistency check. It is not the only reference:
    // each solver has already been compared against the independent analytical
    // quadrature and reconstructed raw moments.
    const std::vector<TerminalNode> dense_nodes = terminal_nodes(dense.quadrature);
    const std::vector<TerminalNode> bp_nodes = terminal_nodes(bp.quadrature);
    if (dense_nodes.size() != bp_nodes.size())
    {
        append_failure(verification, "dense and Björck-Pereyra results have different terminal node counts.");
        return;
    }
    for (std::size_t node = 0; node < dense_nodes.size(); ++node)
    {
        if (dense_nodes[node].path != bp_nodes[node].path)
        {
            append_failure(verification, "dense and Björck-Pereyra results have different terminal paths.");
            continue;
        }
        update_error(verification.errors.dense_bp_node, bp_nodes[node].x, dense_nodes[node].x);
        update_error(verification.errors.dense_bp_node, bp_nodes[node].y, dense_nodes[node].y);
        update_error(verification.errors.dense_bp_weight, bp_nodes[node].weight, dense_nodes[node].weight);
    }

    const qmomkit::RawMomentTensor dense_reconstructed = qmomkit::reconstruct_raw_moments(dense.quadrature, extents);
    const qmomkit::RawMomentTensor bp_reconstructed = qmomkit::reconstruct_raw_moments(bp.quadrature, extents);
    if (dense_reconstructed.values.size() != bp_reconstructed.values.size())
    {
        append_failure(verification, "dense and Björck-Pereyra reconstructed moment counts differ.");
        return;
    }
    for (std::size_t index = 0; index < dense_reconstructed.values.size(); ++index)
    {
        update_error(verification.errors.dense_bp_moment, bp_reconstructed.values[index], dense_reconstructed.values[index]);
    }
}

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::CqmomDenseCase>& configured)
{
    std::cerr << configured.error() << '\n';
}

} // namespace

int main(const int argc, char** argv)
{
    // Verify public Vandermonde solver selection policies on the same small
    // analytical CQMOM problem. The case files differ only by requested solver
    // policy, so any numerical difference or diagnostic mismatch is meaningful.
    if (argc != 2 && argc != 3)
    {
        std::cerr << "Usage: vandermonde_solver_selection_verification <verification/vandermonde directory> [solver]\n";
        return 2;
    }

    const std::filesystem::path case_directory = argv[1];
    const std::vector<SolverCase> cases{{"automatic", case_directory / "automatic.qmom", qmomkit::VandermondeSolver::automatic,
                                         qmomkit::VandermondeSolver::bjorck_pereyra, qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra},
                                        {"dense", case_directory / "dense.qmom", qmomkit::VandermondeSolver::dense, qmomkit::VandermondeSolver::dense,
                                         qmomkit::internal::VandermondeAlgorithm::dense},
                                        {"pivoted_qr", case_directory / "pivoted_qr.qmom", qmomkit::VandermondeSolver::pivoted_qr,
                                         qmomkit::VandermondeSolver::pivoted_qr, qmomkit::internal::VandermondeAlgorithm::pivoted_qr},
                                        {"svd", case_directory / "svd.qmom", qmomkit::VandermondeSolver::svd, qmomkit::VandermondeSolver::svd,
                                         qmomkit::internal::VandermondeAlgorithm::svd},
                                        {"moore_penrose", case_directory / "moore_penrose.qmom", qmomkit::VandermondeSolver::moore_penrose,
                                         qmomkit::VandermondeSolver::moore_penrose, qmomkit::internal::VandermondeAlgorithm::moore_penrose},
                                        {"bjorck_pereyra", case_directory / "bjorck_pereyra.qmom", qmomkit::VandermondeSolver::bjorck_pereyra,
                                         qmomkit::VandermondeSolver::bjorck_pereyra, qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra}};

    const std::string selected_solver = argc == 3 ? argv[2] : "";
    bool matched = selected_solver.empty();
    VerificationResult verification;
    std::vector<qmomkit::CqmomResult> computed;
    std::vector<qmomkit::CqmomDenseCase> loaded_cases;
    computed.reserve(cases.size());
    loaded_cases.reserve(cases.size());

    for (const SolverCase& solver_case : cases)
    {
        if (!selected_solver.empty() && selected_solver != solver_case.name)
        {
            continue;
        }
        matched = true;
        const auto configured = qmomkit::config::load_cqmom_dense(solver_case.filename);
        if (!configured.ok())
        {
            print_configuration_errors(configured);
            return 1;
        }
        if (configured.value().vandermonde_solver != solver_case.requested_solver)
        {
            append_failure(verification, solver_case.name + " configuration did not load the requested solver policy.");
        }

        qmomkit::CqmomDenseCase dense_case = configured.value();
        const qmomkit::CqmomResult result = qmomkit::compute_cqmom(dense_case);
        if (!result.ok())
        {
            append_failure(verification, solver_case.name + " public CQMOM computation failed: " + result.message);
            continue;
        }

        verify_against_reference(solver_case, result, dense_case.moments, verification);
        verify_solver_diagnostics(solver_case, dense_case, verification);
        computed.push_back(result);
        loaded_cases.push_back(std::move(dense_case));
    }

    if (computed.size() == cases.size())
    {
        const qmomkit::CqmomResult* dense_result = nullptr;
        const qmomkit::CqmomResult* bp_result = nullptr;
        const qmomkit::CqmomDenseCase* dense_case = nullptr;
        const qmomkit::CqmomDenseCase* bp_case = nullptr;
        for (std::size_t index = 0; index < loaded_cases.size(); ++index)
        {
            if (loaded_cases[index].vandermonde_solver == qmomkit::VandermondeSolver::dense)
            {
                dense_result = &computed[index];
                dense_case = &loaded_cases[index];
            }
            if (loaded_cases[index].vandermonde_solver == qmomkit::VandermondeSolver::bjorck_pereyra)
            {
                bp_result = &computed[index];
                bp_case = &loaded_cases[index];
            }
        }

        if (dense_result == nullptr || bp_result == nullptr || dense_case == nullptr || bp_case == nullptr)
        {
            append_failure(verification, "dense/Björck-Pereyra comparison could not locate both explicit solver results.");
        }
        else if (!reports_selected_solver(*dense_result, dense_case->hierarchy, qmomkit::VandermondeSolver::dense))
        {
            append_failure(verification, "dense/Björck-Pereyra comparison did not receive the explicit dense result.");
        }
        else if (!reports_selected_solver(*bp_result, bp_case->hierarchy, qmomkit::VandermondeSolver::bjorck_pereyra))
        {
            append_failure(verification, "dense/Björck-Pereyra comparison did not receive the explicit Björck-Pereyra result.");
        }
        else
        {
            compare_dense_and_bp(*dense_result, *bp_result, dense_case->moments.extents, verification);
        }
    }

    if (!matched)
    {
        std::cerr << "Unknown Vandermonde solver case '" << selected_solver << "'.\n";
        return 2;
    }

    const bool passed = verification.passed && verification.errors.reference_node <= reference_tolerance &&
                        verification.errors.reference_weight <= reference_tolerance && verification.errors.reference_moment <= reference_tolerance &&
                        verification.errors.dense_bp_node <= workflow_tolerance && verification.errors.dense_bp_weight <= workflow_tolerance &&
                        verification.errors.dense_bp_moment <= workflow_tolerance;

    std::cout << "QMOMKit Vandermonde solver selection verification\n";
    std::cout << "Reference: fixed 2D independent-product CQMOM system\n";
    std::cout << "Tolerance: " << std::scientific << reference_tolerance << " absolute\n\n";
    std::cout << std::left << std::setw(38) << "solver policy" << std::setw(38) << "selected policy" << "algorithm used\n";
    std::cout << std::string(96, '-') << '\n';
    for (const SolverCase& solver_case : cases)
    {
        std::cout << std::left << std::setw(38) << solver_name(solver_case.requested_solver) << std::setw(38)
                  << solver_name(solver_case.expected_selected_solver) << algorithm_name(solver_case.expected_algorithm) << '\n';
    }
    std::cout << '\n';
    std::cout << std::left << std::setw(30) << "reference node error" << std::right << std::scientific << std::setprecision(4) << std::setw(14)
              << verification.errors.reference_node << '\n';
    std::cout << std::left << std::setw(30) << "reference weight error" << std::right << std::setw(14) << verification.errors.reference_weight
              << '\n';
    std::cout << std::left << std::setw(30) << "reference moment error" << std::right << std::setw(14) << verification.errors.reference_moment
              << '\n';
    std::cout << std::left << std::setw(30) << "dense/bp node difference" << std::right << std::setw(14) << verification.errors.dense_bp_node << '\n';
    std::cout << std::left << std::setw(30) << "dense/bp weight difference" << std::right << std::setw(14) << verification.errors.dense_bp_weight
              << '\n';
    std::cout << std::left << std::setw(30) << "dense/bp moment difference" << std::right << std::setw(14) << verification.errors.dense_bp_moment
              << '\n';
    std::cout << "\nStatus: " << (passed ? "PASS" : "FAIL") << '\n';

    if (!verification.passed)
    {
        std::cerr << verification.message << '\n';
    }
    return passed ? 0 : 1;
}
