#include "../verification/verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
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

std::string nonrealizable_cqmom_configuration()
{
    // Configuration-level non-realizability: all public CQMOM workflows parse
    // successfully, but the root QMOM inversion sees M[0,0] < 0 and must report
    // qmom_failed with branch diagnostics instead of invalid_input.
    return R"([phase_space]
names = x, y
node_counts = 2, 2
conditioning_order = 0, 1

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
    // This contains no realizability policy, so no Hausdorff check may reject
    // the bounded candidate: reconstruction must reach Wheeler and reject only
    // the returned node above the configured quadrature interval.
    return R"([phase_space]
names = x, y
node_counts = 2, 1
conditioning_order = 0, 1
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
    // The same non-realizable mathematical case is run through both public
    // workflows to ensure status mapping is workflow-invariant.
    const std::string configuration = nonrealizable_cqmom_configuration();
    const auto cqmom_prepared_case_result = qmomkit::config::parse_cqmom_prepared(configuration, "nonrealizable_cqmom");
    const auto cqmom_case_result = qmomkit::config::parse_cqmom(configuration, "nonrealizable_cqmom");

    QMOMKIT_VERIFY_CHECK(cqmom_prepared_case_result.ok());
    QMOMKIT_VERIFY_CHECK(cqmom_case_result.ok());

    verify_nonrealizable_result("default CqmomPreparedCase", qmomkit::compute_cqmom(cqmom_prepared_case_result.value()));
    verify_nonrealizable_result("CqmomCase", qmomkit::compute_cqmom(cqmom_case_result.value()));
}

void verify_bounded_abscissa_violation_reaches_wheeler()
{
    const auto configured = qmomkit::config::parse_cqmom_prepared(bounded_abscissa_violation_configuration(), "bounded_abscissa_violation");
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

template <typename Exception, typename Callable>
void verify_extraction_throws(const std::string& case_name, const Callable& extraction, const std::string& expected_message)
{
    bool threw_expected_type = false;
    try
    {
        extraction();
    }
    catch (const Exception& error)
    {
        threw_expected_type = true;
        QMOMKIT_VERIFY_CHECK(contains_text(error.what(), expected_message));
    }
    catch (const std::exception& error)
    {
        // Reported rather than allowed to escape, so that the wrong exception
        // type fails this case with a readable message instead of aborting.
        std::cerr << case_name << " threw an unexpected exception type: " << error.what() << '\n';
        QMOMKIT_VERIFY_CHECK(false);
    }
    QMOMKIT_VERIFY_CHECK(threw_expected_type);

    std::cout << std::left << std::setw(34) << case_name << std::setw(22) << "throws" << expected_message << '\n';
}

void verify_extraction_rejects_structural_mistakes(const std::string& selected_case = {})
{
    // Structural caller mistakes in a full raw tensor - rank mismatch, wrong
    // value count, or a tensor too small for the moments the hierarchy needs -
    // are rejected when the required moments are extracted, before any case is
    // built. Extraction reports them as exceptions rather than statuses because
    // there is no result object to carry a status at that point.
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);

    if (selected_case.empty() || selected_case == "extraction_extent_rank")
    {
        verify_extraction_throws<std::invalid_argument>(
            "extraction extent rank mismatch",
            [&] { return qmomkit::extract_cqmom_moments(qmomkit::RawMomentTensor{{4}, std::vector<double>(4, 1.0)}, layout); },
            "The rank of the full tensor of raw moments must match the CQMOM moment layout dimension count.");
    }

    if (selected_case.empty() || selected_case == "extraction_value_count")
    {
        verify_extraction_throws<std::invalid_argument>(
            "extraction value count mismatch",
            [&] { return qmomkit::extract_cqmom_moments(qmomkit::RawMomentTensor{{4, 4}, std::vector<double>(15, 1.0)}, layout); },
            "The value count of the full tensor of raw moments does not match its extents.");
    }

    if (selected_case.empty() || selected_case == "extraction_missing_required_moment")
    {
        verify_extraction_throws<std::out_of_range>(
            "extraction missing required moment",
            [&] { return qmomkit::extract_cqmom_moments(qmomkit::RawMomentTensor{{3, 4}, indexed_values({3, 4})}, layout); },
            "The full tensor of raw moments does not contain a required CQMOM moment.");
    }
}

void verify_moment_data_invalid_result(const std::string& case_name, const qmomkit::CqmomResult& result, const std::string& expected_message)
{
    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::CqmomStatus::invalid_input);
    QMOMKIT_VERIFY_CHECK(!result.message.empty());
    QMOMKIT_VERIFY_CHECK(contains_text(result.message, expected_message));

    std::cout << std::left << std::setw(34) << case_name << std::setw(22) << "invalid_input" << expected_message << '\n';
}

void verify_moment_data_invalid_inputs(const std::string& selected_case = {})
{
    // Invalid moment-data cases exercise the required-moment API:
    // missing values, non-finite required moments, and a layout whose hierarchy
    // no longer matches the stored data.
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);

    if (selected_case.empty() || selected_case == "moment_data_missing_required_value")
    {
        qmomkit::CqmomMomentData missing_values = qmomkit::extract_cqmom_moments({extents, indexed_values(extents)}, layout);
        missing_values.values.pop_back();
        verify_moment_data_invalid_result("moment data missing value", qmomkit::compute_cqmom(qmomkit::CqmomCase{missing_values}),
                                          "values must match");
    }

    if (selected_case.empty() || selected_case == "moment_data_nonfinite_value")
    {
        qmomkit::CqmomMomentData nonfinite_values = qmomkit::extract_cqmom_moments({extents, indexed_values(extents)}, layout);
        nonfinite_values.values[0] = std::numeric_limits<double>::quiet_NaN();
        verify_moment_data_invalid_result("moment data non-finite value", qmomkit::compute_cqmom(qmomkit::CqmomCase{nonfinite_values}), "finite");
    }

    if (selected_case.empty() || selected_case == "moment_data_incompatible_layout")
    {
        const qmomkit::PhaseSpaceHierarchy incompatible_hierarchy({2, 2, 2}, {2, 0, 1});
        const qmomkit::CqmomMomentLayout incompatible_layout = qmomkit::create_cqmom_moment_layout(incompatible_hierarchy);
        qmomkit::CqmomMomentData incompatible = qmomkit::extract_cqmom_moments({{4, 4, 4}, indexed_values({4, 4, 4})}, incompatible_layout);
        incompatible.layout.hierarchy = qmomkit::PhaseSpaceHierarchy({2, 2, 2});
        verify_moment_data_invalid_result("moment data incompatible layout", qmomkit::compute_cqmom(qmomkit::CqmomCase{incompatible}), "layout");
    }
}

} // namespace

int main(const int argc, char** argv)
{
    // Direct CQMOM compute failures distinguish invalid public requests from
    // mathematically non-realizable moment sets, and confirm that
    // CqmomPreparedCase and CqmomCase report the same failure where applicable.
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
        verify_extraction_rejects_structural_mistakes();
        verify_moment_data_invalid_inputs();
    }
    else if (selected_case == "nonrealizable_public_workflow")
        verify_nonrealizable_public_workflows();
    else if (selected_case == "bounded_abscissa_violation")
        verify_bounded_abscissa_violation_reaches_wheeler();
    else if (selected_case.rfind("extraction_", 0) == 0)
        verify_extraction_rejects_structural_mistakes(selected_case);
    else if (selected_case.rfind("moment_data_", 0) == 0)
        verify_moment_data_invalid_inputs(selected_case);
    else
    {
        std::cerr << "Unknown CQMOM failure case '" << selected_case << "'.\n";
        return 2;
    }

    return 0;
}
