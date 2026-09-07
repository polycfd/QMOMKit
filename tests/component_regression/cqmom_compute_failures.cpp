#include "../verification/verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{

bool contains_text(const std::string& text, const std::string& fragment)
{
    return text.find(fragment) != std::string::npos;
}

std::vector<double> indexed_values(const std::vector<std::size_t>& extents)
{
    std::vector<double> values(qmomkit::product_extents(extents), 0.0);
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        values[i] = 1.0 + static_cast<double>(i);
    }
    return values;
}

std::vector<double> raw_moments_from_terminals(const std::vector<std::size_t>& extents, const std::vector<double>& weights,
                                               const std::vector<std::vector<double>>& abscissas)
{
    // Build a realizable dense tensor directly from terminal nodes and weights
    // when the failure being tested is policy validation rather than moment
    // mathematics.
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        double moment = 0.0;
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            double contribution = weights[node];
            for (std::size_t dimension = 0; dimension < moment_index.size(); ++dimension)
            {
                contribution *= std::pow(abscissas[dimension][node], static_cast<double>(moment_index[dimension]));
            }
            moment += contribution;
        }
        moments[flat_moment] = moment;
    }
    return moments;
}

std::string nonrealizable_cqmom_configuration()
{
    // Configuration-level non-realizability: all public CQMOM workflows parse
    // successfully, but the root QMOM inversion sees M[0,0] < 0 and must report
    // qmom_failed with branch diagnostics instead of invalid_input.
    return R"([phase_space]
names = x, y
node_counts = 2, 2
conditioning_order = 0, 1
adaptive_branch_counts = false

[qmom]
adaptive = false

[moments]
extents = 4, 4
M[0,0] = -1.0
M[0,1] = 1.25
M[0,2] = 3.25
M[0,3] = 5.75
M[1,0] = 2.2
M[1,1] = 2.75
M[1,2] = 7.15
M[1,3] = 12.65
M[2,0] = 5.8
M[2,1] = 7.25
M[2,2] = 18.85
M[2,3] = 33.35
M[3,0] = 16.6
M[3,1] = 20.75
M[3,2] = 53.95
M[3,3] = 95.45
)";
}

std::string bounded_abscissa_violation_configuration()
{
    // This contains no realizability policy. Before the separation cleanup,
    // the bounded candidate was rejected by a Hausdorff check before Wheeler
    // could return its nodes. The current path must reach Wheeler and reject
    // only the returned node above the configured quadrature interval.
    return R"([phase_space]
names = x, y
node_counts = 2, 1
conditioning_order = 0, 1
adaptive_branch_counts = false
support_lower = 0.0, none
support_upper = 1.0, none

[qmom]
adaptive = false

[moments]
extents = 4, 2
M[0,0] = 1.0
M[0,1] = 0.0
M[1,0] = 0.75
M[1,1] = 0.0
M[2,0] = 0.8125
M[2,1] = 0.0
M[3,0] = 0.984375
M[3,1] = 0.0
)";
}

void verify_nonrealizable_result(const std::string& workflow_name, const qmomkit::CqmomResult& result)
{
    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::CqmomStatus::qmom_failed);
    QMOMKIT_VERIFY_CHECK(!result.message.empty());
    QMOMKIT_VERIFY_CHECK(contains_text(result.message, "m0 < 0"));
    QMOMKIT_VERIFY_CHECK(!result.diagnostics.empty());

    const qmomkit::CqmomDiagnostics& diagnostics = result.diagnostics.front();
    QMOMKIT_VERIFY_CHECK(diagnostics.conditioning_level == 0);
    QMOMKIT_VERIFY_CHECK(diagnostics.phase_space_dimension == 0);
    QMOMKIT_VERIFY_CHECK(diagnostics.qmom_status == qmomkit::QmomStatus::not_realizable);
    QMOMKIT_VERIFY_CHECK(!diagnostics.message.empty());
    QMOMKIT_VERIFY_CHECK(contains_text(diagnostics.message, "m0 < 0"));

    std::cout << std::left << std::setw(34) << workflow_name << std::setw(22) << "qmom_failed" << "root m0 < 0\n";
}

void verify_nonrealizable_public_workflows()
{
    // The same non-realizable mathematical case is run through default, dense,
    // and lean public workflows to ensure status mapping is workflow-invariant.
    const std::string configuration = nonrealizable_cqmom_configuration();
    const auto default_case = qmomkit::config::parse_cqmom(configuration, "nonrealizable_cqmom");
    const auto dense_case = qmomkit::config::parse_cqmom_dense(configuration, "nonrealizable_cqmom");
    const auto lean_case = qmomkit::config::parse_cqmom_lean(configuration, "nonrealizable_cqmom");

    QMOMKIT_VERIFY_CHECK(default_case.ok());
    QMOMKIT_VERIFY_CHECK(dense_case.ok());
    QMOMKIT_VERIFY_CHECK(lean_case.ok());

    verify_nonrealizable_result("default CqmomCase", qmomkit::compute_cqmom(default_case.value()));
    verify_nonrealizable_result("dense CqmomDenseCase", qmomkit::compute_cqmom(dense_case.value()));
    verify_nonrealizable_result("lean CqmomLeanCase", qmomkit::compute_cqmom(lean_case.value()));
}

void verify_bounded_abscissa_violation_reaches_wheeler()
{
    const auto configured = qmomkit::config::parse_cqmom(bounded_abscissa_violation_configuration(), "bounded_abscissa_violation");
    QMOMKIT_VERIFY_CHECK(configured.ok());

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(configured.value());
    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::CqmomStatus::qmom_failed);
    QMOMKIT_VERIFY_CHECK(result.failed_phase_space_dimension == 0);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.size() == 1);

    const qmomkit::CqmomDiagnostics& diagnostic = result.diagnostics.front();
    QMOMKIT_VERIFY_CHECK(diagnostic.phase_space_dimension == 0);
    QMOMKIT_VERIFY_CHECK(diagnostic.bounded_support_checked);
    QMOMKIT_VERIFY_CHECK(diagnostic.bounded_candidate_node_count == 2);
    QMOMKIT_VERIFY_CHECK(diagnostic.bounded_returned_nodes.size() == 2);
    QMOMKIT_VERIFY_CHECK(diagnostic.bounded_returned_weights.size() == 2);
    QMOMKIT_VERIFY_NEAR(diagnostic.bounded_returned_nodes[0], 0.25, 1.0e-12);
    QMOMKIT_VERIFY_NEAR(diagnostic.bounded_returned_nodes[1], 1.25, 1.0e-12);
    QMOMKIT_VERIFY_NEAR(diagnostic.bounded_returned_weights[0], 0.5, 1.0e-12);
    QMOMKIT_VERIFY_NEAR(diagnostic.bounded_returned_weights[1], 0.5, 1.0e-12);
    QMOMKIT_VERIFY_CHECK(contains_text(result.message, "returned nodes outside the configured interval"));
    QMOMKIT_VERIFY_CHECK(!contains_text(result.message, "Bounded-support realizability check failed"));
    QMOMKIT_VERIFY_CHECK(!contains_text(result.message, "Support: bounded_interval"));
    QMOMKIT_VERIFY_CHECK(!contains_text(result.message, "L_lower"));
    QMOMKIT_VERIFY_CHECK(!contains_text(diagnostic.message, "Bounded-support realizability check failed"));
    QMOMKIT_VERIFY_CHECK(!contains_text(diagnostic.message, "Support: bounded_interval"));
    QMOMKIT_VERIFY_CHECK(!contains_text(diagnostic.message, "L_lower"));

    std::cout << std::left << std::setw(34) << "bounded abscissa violation" << std::setw(22) << "qmom_failed"
              << "returned node above interval\n";
}

void verify_dense_invalid_result(const std::string& case_name, const qmomkit::CqmomResult& result, const std::string& expected_message)
{
    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::CqmomStatus::invalid_input);
    QMOMKIT_VERIFY_CHECK(!result.message.empty());
    QMOMKIT_VERIFY_CHECK(contains_text(result.message, expected_message));
    QMOMKIT_VERIFY_CHECK(result.diagnostics.empty());

    std::cout << std::left << std::setw(34) << case_name << std::setw(22) << "invalid_input" << expected_message << '\n';
}

void verify_dense_invalid_inputs(const std::string& selected_case = {})
{
    // Dense invalid-input cases represent structural caller mistakes in a full
    // raw tensor: rank mismatch, wrong value count, non-finite data, or missing
    // moments required by the configured hierarchy.
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});

    if (selected_case.empty() || selected_case == "dense_extent_rank")
    {
        verify_dense_invalid_result(
            "dense extent rank mismatch",
            qmomkit::compute_cqmom(qmomkit::CqmomDenseCase{hierarchy, qmomkit::RawMomentTensor{{4}, std::vector<double>(4, 1.0)}, {}, {}}),
            "extents must match");
    }

    if (selected_case.empty() || selected_case == "dense_value_count")
    {
        verify_dense_invalid_result(
            "dense value count mismatch",
            qmomkit::compute_cqmom(qmomkit::CqmomDenseCase{hierarchy, qmomkit::RawMomentTensor{{4, 4}, std::vector<double>(15, 1.0)}, {}, {}}),
            "data size");
    }

    if (selected_case.empty() || selected_case == "dense_nonfinite_moment")
    {
        qmomkit::RawMomentTensor nonfinite{{4, 4}, std::vector<double>(16, 1.0)};
        nonfinite.values[0] = std::numeric_limits<double>::infinity();
        verify_dense_invalid_result("dense non-finite moment", qmomkit::compute_cqmom(qmomkit::CqmomDenseCase{hierarchy, nonfinite, {}, {}}),
                                    "finite");
    }

    if (selected_case.empty() || selected_case == "dense_missing_required_moment")
    {
        verify_dense_invalid_result(
            "dense missing required moment",
            qmomkit::compute_cqmom(qmomkit::CqmomDenseCase{hierarchy, qmomkit::RawMomentTensor{{3, 4}, indexed_values({3, 4})}, {}, {}}),
            "at least 2*N");
    }
}

void verify_lean_invalid_result(const std::string& case_name, const qmomkit::CqmomResult& result, const std::string& expected_message)
{
    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::CqmomStatus::invalid_input);
    QMOMKIT_VERIFY_CHECK(!result.message.empty());
    QMOMKIT_VERIFY_CHECK(contains_text(result.message, expected_message));

    std::cout << std::left << std::setw(34) << case_name << std::setw(22) << "invalid_input" << expected_message << '\n';
}

void verify_lean_invalid_inputs(const std::string& selected_case = {})
{
    // Lean invalid-input cases exercise the selected required-moment API:
    // missing values, non-finite required moments, and a layout whose hierarchy
    // no longer matches the stored data.
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);

    if (selected_case.empty() || selected_case == "lean_missing_required_value")
    {
        qmomkit::CqmomMomentData missing_values = qmomkit::extract_cqmom_moments({extents, indexed_values(extents)}, layout);
        missing_values.values.pop_back();
        verify_lean_invalid_result("lean missing required value", qmomkit::compute_cqmom(qmomkit::CqmomLeanCase{missing_values}),
                                   "values must match");
    }

    if (selected_case.empty() || selected_case == "lean_nonfinite_value")
    {
        qmomkit::CqmomMomentData nonfinite_values = qmomkit::extract_cqmom_moments({extents, indexed_values(extents)}, layout);
        nonfinite_values.values[0] = std::numeric_limits<double>::quiet_NaN();
        verify_lean_invalid_result("lean non-finite value", qmomkit::compute_cqmom(qmomkit::CqmomLeanCase{nonfinite_values}), "finite");
    }

    if (selected_case.empty() || selected_case == "lean_incompatible_layout")
    {
        const qmomkit::PhaseSpaceHierarchy incompatible_hierarchy({2, 2, 2}, {2, 0, 1});
        const qmomkit::CqmomMomentLayout incompatible_layout = qmomkit::create_cqmom_moment_layout(incompatible_hierarchy);
        qmomkit::CqmomMomentData incompatible = qmomkit::extract_cqmom_moments({{4, 4, 4}, indexed_values({4, 4, 4})}, incompatible_layout);
        incompatible.layout.hierarchy = qmomkit::PhaseSpaceHierarchy({2, 2, 2});
        verify_lean_invalid_result("lean incompatible layout", qmomkit::compute_cqmom(qmomkit::CqmomLeanCase{incompatible}), "layout");
    }
}

void verify_adaptive_policy_rejection()
{
    // Adaptive Wheeler can reduce the root node count. When adaptive branch
    // counts are disabled, that reduction would make a fixed rectangular CQMOM
    // hierarchy ambiguous, so the public API rejects the policy combination
    // before running a reconstruction.
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.1, 0.2, 0.3, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 1.0, 1.0}, {-2.0, 2.0, -2.0, 2.0}};
    const std::vector<double> moments = raw_moments_from_terminals(extents, weights, abscissas);
    const qmomkit::QmomOptions adaptive_qmom{true, 1.0e-10, 1.0e-10, 1.0e-12};

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom({hierarchy, {extents, moments}, {}, adaptive_qmom, false});

    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::CqmomStatus::invalid_input);
    QMOMKIT_VERIFY_CHECK(result.failed_phase_space_dimension == 0);
    QMOMKIT_VERIFY_CHECK(result.diagnostics.empty());
    QMOMKIT_VERIFY_CHECK(contains_text(result.message, "adaptive_branch_counts"));

    std::cout << std::left << std::setw(34) << "adaptive policy mismatch" << std::setw(22) << "invalid_input"
              << "adaptive_branch_counts\n";
}

} // namespace

int main(const int argc, char** argv)
{
    // Direct CQMOM compute failures distinguish invalid public requests from
    // mathematically non-realizable moment sets, and confirm that dense, lean,
    // and default workflows expose the same failure meaning where applicable.
    std::cout << "QMOMKit direct CQMOM compute failure verification\n";
    std::cout << std::left << std::setw(34) << "case" << std::setw(22) << "expected status" << "expected reason\n";
    std::cout << std::string(82, '-') << '\n';

    const std::string selected_case = argc == 2 ? argv[1] : "";
    if (argc > 2)
    {
        std::cerr << "Usage: cqmom_compute_failures [case-name]\n";
        return 2;
    }
    if (selected_case.empty())
    {
        verify_nonrealizable_public_workflows();
        verify_bounded_abscissa_violation_reaches_wheeler();
        verify_dense_invalid_inputs();
        verify_lean_invalid_inputs();
        verify_adaptive_policy_rejection();
    }
    else if (selected_case == "nonrealizable_public_workflow")
        verify_nonrealizable_public_workflows();
    else if (selected_case == "bounded_abscissa_violation")
        verify_bounded_abscissa_violation_reaches_wheeler();
    else if (selected_case.rfind("dense_", 0) == 0)
        verify_dense_invalid_inputs(selected_case);
    else if (selected_case.rfind("lean_", 0) == 0)
        verify_lean_invalid_inputs(selected_case);
    else if (selected_case == "adaptive_policy")
        verify_adaptive_policy_rejection();
    else
    {
        std::cerr << "Unknown CQMOM failure case '" << selected_case << "'.\n";
        return 2;
    }

    return 0;
}
