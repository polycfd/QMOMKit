// Rebuilding raw moments from a CQMOM quadrature, and refusing to do so when
// the request does not match the quadrature.
//
// What is verified
//   That moments rebuilt from a recovered quadrature match the moments it came
//   from, for a dense QMOM rule, for a selected subset of CQMOM indices, and
//   for the full set a hierarchy requires; and that mismatched requests are
//   rejected rather than answered with a partially filled tensor.
//
// How the reference is constructed
//   Chosen, not computed. A rule of weighted points is written down, its raw
//   moments are evaluated from it, and those moments are what the
//   reconstruction must return. Because the quadrature is the starting point
//   rather than something the library derived, a reconstruction that silently
//   dropped or reordered an index would show up as a moment mismatch.
//
// Error metric and pass criteria
//   Maximum absolute deviation of each rebuilt moment from the prescribed
//   value, at 1.0e-12. The failure cases additionally require a non-success
//   status rather than a best-effort answer, because a partially filled tensor
//   that looks plausible is the dangerous outcome here.
//
// Reproduce
//   ctest -R verification_.*moment_reconstruction
//
// Build modes
//   Both. Failures are reported through an exit code, not assert().

#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace
{

constexpr double tolerance = 1.0e-12;

struct VerificationResult
{
    bool passed = true;
    std::string message;
    double maximum_error = 0.0;
};

struct TerminalNode
{
    double x = 0.0;
    double y = 0.0;
    double weight = 0.0;
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

bool contains_text(const std::string& text, const std::string& fragment)
{
    return text.find(fragment) != std::string::npos;
}

double integer_power(const double base, const std::size_t exponent)
{
    double value = 1.0;
    for (std::size_t i = 0; i < exponent; ++i)
    {
        value *= base;
    }
    return value;
}

void update_error(VerificationResult& result, const double actual, const double expected)
{
    result.maximum_error = std::max(result.maximum_error, std::abs(actual - expected));
}

std::vector<TerminalNode> independent_product_reference()
{
    // Independent analytical reference: x has nodes {1, 3} with weights
    // {0.4, 0.6}; y has nodes {-1, 2} with weights {0.25, 0.75}. Terminal
    // weights are products of the marginal weights.
    return {{1.0, -1.0, 0.4 * 0.25}, {1.0, 2.0, 0.4 * 0.75}, {3.0, -1.0, 0.6 * 0.25}, {3.0, 2.0, 0.6 * 0.75}};
}

double analytical_cqmom_moment(const qmomkit::MultiIndexView index)
{
    double moment = 0.0;
    for (const TerminalNode& node : independent_product_reference())
    {
        moment += node.weight * integer_power(node.x, index[0]) * integer_power(node.y, index[1]);
    }
    return moment;
}

double analytical_qmom_moment(const std::size_t order)
{
    // One-dimensional independent reference used only for the public
    // reconstruct_raw_moments(Quadrature1D, count) path.
    const std::vector<double> nodes{-1.0, 0.5, 2.0};
    const std::vector<double> weights{0.25, 0.5, 0.25};
    double moment = 0.0;
    for (std::size_t node = 0; node < nodes.size(); ++node)
    {
        moment += weights[node] * integer_power(nodes[node], order);
    }
    return moment;
}

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::CqmomCase>& configured)
{
    std::cerr << configured.error() << '\n';
}

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::CqmomLeanCase>& configured)
{
    std::cerr << configured.error() << '\n';
}

void verify_qmom_dense_reconstruction(VerificationResult& verification)
{
    // Public 1D reconstruction is verified against an independent three-node
    // discrete quadrature. The failure subcase simulates a caller accidentally
    // passing diagnostic data from an unsuccessful inversion.
    const qmomkit::Quadrature1D quadrature{{-1.0, 0.5, 2.0}, {0.25, 0.5, 0.25}, qmomkit::QuadratureStatus::success, ""};
    const std::vector<double> moments = qmomkit::reconstruct_raw_moments(quadrature, 6);
    if (moments.size() != 6)
    {
        append_failure(verification, "QMOM dense reconstruction returned the wrong number of moments.");
        return;
    }
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        update_error(verification, moments[order], analytical_qmom_moment(order));
    }

    try
    {
        const qmomkit::Quadrature1D failed{{0.0}, {1.0}, qmomkit::QuadratureStatus::not_realizable, "synthetic failed quadrature"};
        (void)qmomkit::reconstruct_raw_moments(failed, 2);
        append_failure(verification, "QMOM reconstruction accepted an unsuccessful quadrature.");
    }
    catch (const std::exception& error)
    {
        if (!contains_text(error.what(), "unsuccessful QMOM quadrature"))
        {
            append_failure(verification, "QMOM reconstruction failure message did not explain the unsuccessful quadrature.");
        }
    }
}

void verify_selected_cqmom_reconstruction(const qmomkit::CqmomQuadrature& quadrature, VerificationResult& verification)
{
    // Selected reconstruction exists for callers that need a sparse, ordered
    // moment list rather than a dense tensor. This request intentionally mixes
    // ordering, duplicates, and moment orders so the result must preserve the
    // user's request exactly.
    const qmomkit::RawMomentRequest request{{{2, 0}, {0, 0}, {1, 1}, {0, 2}, {2, 0}, {3, 1}}};
    const qmomkit::SelectedRawMomentResult selected = qmomkit::reconstruct_raw_moments(quadrature, request);
    if (!selected.ok())
    {
        append_failure(verification, "selected CQMOM reconstruction did not report success: " + selected.message);
        return;
    }
    if (selected.indices != request.indices)
    {
        append_failure(verification, "selected CQMOM reconstruction did not preserve request index order.");
    }
    if (selected.values.size() != request.indices.size())
    {
        append_failure(verification, "selected CQMOM reconstruction returned the wrong number of values.");
        return;
    }
    if (selected.values[0] != selected.values[4])
    {
        append_failure(verification, "selected CQMOM reconstruction did not preserve duplicate requested moments.");
    }

    const qmomkit::RawMomentTensor dense = qmomkit::reconstruct_raw_moments(quadrature, {4, 4});
    for (std::size_t entry = 0; entry < request.indices.size(); ++entry)
    {
        const qmomkit::MultiIndex& index = request.indices[entry];
        const double analytical = analytical_cqmom_moment(index);
        const double from_dense = dense.values[qmomkit::flatten_multi_index(index, dense.extents)];
        update_error(verification, selected.values[entry], analytical);
        update_error(verification, selected.values[entry], from_dense);
    }
}

void verify_cqmom_required_layout_reconstruction(const qmomkit::CqmomQuadrature& quadrature, const qmomkit::CqmomMomentLayout& layout,
                                                 VerificationResult& verification)
{
    // reconstruct_cqmom_moments is the public way to reconstruct exactly the
    // fixed recursive CQMOM moment layout. Values must align with
    // layout.required_indices, not with dense tensor order.
    const qmomkit::CqmomMomentData reconstructed = qmomkit::reconstruct_cqmom_moments(quadrature, layout);
    if (reconstructed.layout.required_indices != layout.required_indices)
    {
        append_failure(verification, "CQMOM required-layout reconstruction did not preserve layout index order.");
    }
    if (reconstructed.values.size() != layout.required_indices.size())
    {
        append_failure(verification, "CQMOM required-layout reconstruction returned the wrong number of values.");
        return;
    }

    const qmomkit::RawMomentTensor dense = qmomkit::reconstruct_raw_moments(quadrature, layout.hierarchy.moment_extents());
    for (std::size_t entry = 0; entry < layout.required_indices.size(); ++entry)
    {
        const qmomkit::MultiIndexView index = layout.required_indices[entry];
        const double analytical = analytical_cqmom_moment(index);
        const double from_dense = dense.values[qmomkit::flatten_multi_index(index, dense.extents)];
        update_error(verification, reconstructed.values[entry], analytical);
        update_error(verification, reconstructed.values[entry], from_dense);
    }
}

void verify_cqmom_reconstruction_failures(const qmomkit::CqmomQuadrature& quadrature, const qmomkit::CqmomMomentLayout& layout,
                                          VerificationResult& verification)
{
    // Failure paths model invalid public reconstruction requests: a selected
    // moment with the wrong phase-space rank and a required-layout request that
    // does not match the quadrature hierarchy. Both should fail with useful
    // messages instead of silently returning misindexed values.
    try
    {
        const qmomkit::RawMomentRequest wrong_rank{{{2, 0, 1}}};
        (void)qmomkit::reconstruct_raw_moments(quadrature, wrong_rank);
        append_failure(verification, "selected CQMOM reconstruction accepted a request with the wrong rank.");
    }
    catch (const std::exception& error)
    {
        if (!contains_text(error.what(), "rank"))
        {
            append_failure(verification, "selected CQMOM reconstruction wrong-rank failure message was not meaningful.");
        }
    }

    try
    {
        qmomkit::CqmomMomentLayout mismatched = layout;
        mismatched.hierarchy = qmomkit::PhaseSpaceHierarchy({2, 1}, {0, 1});
        (void)qmomkit::reconstruct_cqmom_moments(quadrature, mismatched);
        append_failure(verification, "CQMOM required-layout reconstruction accepted a mismatched hierarchy.");
    }
    catch (const std::exception& error)
    {
        if (!contains_text(error.what(), "hierarchy"))
        {
            append_failure(verification, "CQMOM required-layout mismatch failure message was not meaningful.");
        }
    }
}

} // namespace

int main(const int argc, char** argv)
{
    // Reconstruct raw moments from independent 1D and 2D quadratures using the
    // public dense, selected, and CQMOM-required-layout APIs. The reference
    // values are evaluated analytically from the defining terminal rules.
    if (argc != 2 && argc != 3)
    {
        std::cerr << "Usage: moment_reconstruction_verification <independent_product.qmom> [case-name]\n";
        return 2;
    }

    const std::filesystem::path configuration_path = argv[1];
    const auto configured = qmomkit::config::load_cqmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        return 1;
    }
    const auto lean = qmomkit::config::load_cqmom_lean(configuration_path);
    if (!lean.ok())
    {
        print_configuration_errors(lean);
        return 1;
    }

    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(configured.value());
    if (!cqmom.ok())
    {
        std::cerr << "Unable to compute reference CQMOM quadrature: " << cqmom.message << '\n';
        return 1;
    }

    const std::string selected_case = argc == 3 ? argv[2] : "";
    const auto selected = [&selected_case](const std::string& name) { return selected_case.empty() || selected_case == name; };
    bool matched = selected_case.empty();
    VerificationResult verification;
    if (selected("qmom_dense"))
    {
        matched = true;
        verify_qmom_dense_reconstruction(verification);
    }
    if (selected("cqmom_selected"))
    {
        matched = true;
        verify_selected_cqmom_reconstruction(cqmom.quadrature, verification);
    }
    if (selected("cqmom_required_layout"))
    {
        matched = true;
        verify_cqmom_required_layout_reconstruction(cqmom.quadrature, lean.value().moments.layout, verification);
    }
    if (selected("failure_paths"))
    {
        matched = true;
        verify_cqmom_reconstruction_failures(cqmom.quadrature, lean.value().moments.layout, verification);
    }
    if (!matched)
    {
        std::cerr << "Unknown moment-reconstruction case '" << selected_case << "'.\n";
        return 2;
    }

    if (verification.maximum_error > tolerance)
    {
        append_failure(verification, "raw moment reconstruction exceeded the analytical tolerance.");
    }
    // Silent on success; CTest surfaces stderr with --output-on-failure.
    if (!verification.passed)
    {
        std::cerr << "Configuration: " << configuration_path << '\n';
        std::cerr << "Maximum error: " << verification.maximum_error << " against tolerance " << tolerance << '\n';
        std::cerr << verification.message << '\n';
        return 1;
    }
    return 0;
}
