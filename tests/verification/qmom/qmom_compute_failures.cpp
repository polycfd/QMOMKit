// Failure reporting for moment sets that cannot produce a quadrature.
//
// What is verified
//   That each invalid input produces a specific non-success status and a
//   message naming the cause, rather than a silent fallback or a quadrature
//   built from unusable data.
//
// How the reference is constructed
//   Chosen by construction: each case is a moment set that violates a stated
//   requirement - too few moments for the requested node count (also with
//   abscissa bounds and adaptive inversion), a negative zeroth moment, a
//   negative variance, a zero node count. The expected status
//   and a fragment of the expected message are written next to each case.
//   There is no numerical reference here because no quadrature should exist.
//
// Pass criteria
//   The returned status must equal the expected status, and the returned
//   message must contain the expected fragment. A match on status alone would
//   not show that the diagnostic points at the right cause.
//
// Reproduce
//   ctest -R verification_qmom_compute_failure
//
// Build modes
//   Both. Checks go through QMOMKIT_VERIFY_*, not assert().

#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

namespace
{

struct FailureCase
{
    std::string name;
    qmomkit::MomentSet moments;
    std::size_t node_count = 0;
    qmomkit::QmomStatus expected_status = qmomkit::QmomStatus::invalid_input;
    std::string expected_message_fragment;
    std::optional<qmomkit::AbscissaBounds> bounds;
    bool adaptive = false;
};

bool contains_text(const std::string& text, const std::string& fragment)
{
    return text.find(fragment) != std::string::npos;
}

void verify_failure_case(const FailureCase& failure_case)
{
    // Direct public compute_qmom failure cases model common caller mistakes or
    // impossible moment sets: invalid request shape, negative number density,
    // and violated positivity. The expected behavior is a clear non-success status
    // and message, with no successful quadrature implied by diagnostic data.
    // The requests are fixed-count unless a case says otherwise: adaptive
    // inversion would recover a reduced rule from the negative variance
    // instead of reporting it.
    const qmomkit::QmomOptions options{failure_case.adaptive};
    const qmomkit::QmomResult result = failure_case.bounds.has_value()
                                           ? qmomkit::compute_qmom(failure_case.moments, failure_case.node_count, options, *failure_case.bounds)
                                           : qmomkit::compute_qmom(failure_case.moments, failure_case.node_count, options);

    QMOMKIT_VERIFY_CHECK(!result.ok());
    QMOMKIT_VERIFY_CHECK(result.status == failure_case.expected_status);
    QMOMKIT_VERIFY_CHECK(!result.message.empty());
    QMOMKIT_VERIFY_CHECK(contains_text(result.message, failure_case.expected_message_fragment));

    if (result.status == qmomkit::QmomStatus::invalid_input)
    {
        // Structurally invalid public requests should not expose fallback
        // quadrature data. Non-realizable moment sets may carry diagnostic
        // fallback data, but they must still report non-success.
        QMOMKIT_VERIFY_CHECK(result.abscissas.empty());
        QMOMKIT_VERIFY_CHECK(result.weights.empty());
    }
    else
    {
        QMOMKIT_VERIFY_CHECK(result.status != qmomkit::QmomStatus::success);
    }
}

} // namespace

int main(const int argc, char** argv)
{
    // These cases are intentionally small so the mathematical reason for each
    // failure is visible in the moment sequence itself.
    const std::vector<FailureCase> cases{
        {"negative zeroth moment", qmomkit::MomentSet{-1.0, 0.0, 1.0, 0.0}, 2, qmomkit::QmomStatus::not_realizable, "m0 < 0"},
        {"negative variance", qmomkit::MomentSet{1.0, 0.0, -1.0, 0.0}, 2, qmomkit::QmomStatus::not_realizable, "negative recurrence"},
        {"insufficient moments", qmomkit::MomentSet{1.0, 0.0, 1.0}, 2, qmomkit::QmomStatus::invalid_input, "fewer than 2*n"},
        // Bounds and adaptive inversion do not answer the request with the
        // one-node rule the three moments support.
        {"insufficient moments bounded adaptive", qmomkit::MomentSet{1.0, 0.0, 1.0}, 2, qmomkit::QmomStatus::invalid_input, "fewer than 2*n",
         qmomkit::AbscissaBounds{-10.0, 10.0}, true},
        {"zero node count", qmomkit::MomentSet{1.0, 0.0}, 0, qmomkit::QmomStatus::invalid_input, "positive"},
    };

    const std::string selected_case = argc == 2 ? argv[1] : "";
    if (argc > 2)
    {
        std::cerr << "Usage: qmom_compute_failures [case-name]\n";
        return 2;
    }

    bool matched = selected_case.empty();
    for (const FailureCase& failure_case : cases)
    {
        if (!selected_case.empty() && failure_case.name != selected_case)
        {
            continue;
        }
        matched = true;
        verify_failure_case(failure_case);
    }

    if (!matched)
    {
        std::cerr << "Unknown QMOM failure case '" << selected_case << "'.\n";
        return 2;
    }

    return 0;
}
