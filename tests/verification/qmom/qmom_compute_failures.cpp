#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"

#include <cstddef>
#include <iomanip>
#include <iostream>
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
};

bool contains_text(const std::string& text, const std::string& fragment)
{
    return text.find(fragment) != std::string::npos;
}

void verify_failure_case(const FailureCase& failure_case)
{
    // Direct public compute_qmom failure cases model common caller mistakes or
    // impossible moment sets: invalid request shape, negative total mass, and
    // violated positivity. The expected behavior is a clear non-success status
    // and message, with no successful quadrature implied by diagnostic data.
    const qmomkit::QmomResult result = qmomkit::compute_qmom(failure_case.moments, failure_case.node_count);

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
        {"zero node count", qmomkit::MomentSet{1.0, 0.0}, 0, qmomkit::QmomStatus::invalid_input, "positive"},
    };

    std::cout << "QMOMKit direct QMOM compute failure verification\n";
    std::cout << std::left << std::setw(28) << "case" << std::setw(20) << "expected status" << "message fragment\n";
    std::cout << std::string(72, '-') << '\n';

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
        std::cout << std::left << std::setw(28) << failure_case.name << std::setw(20);
        switch (failure_case.expected_status)
        {
        case qmomkit::QmomStatus::invalid_input:
            std::cout << "invalid_input";
            break;
        case qmomkit::QmomStatus::not_realizable:
            std::cout << "not_realizable";
            break;
        case qmomkit::QmomStatus::success:
            std::cout << "success";
            break;
        }
        std::cout << failure_case.expected_message_fragment << '\n';
    }

    if (!matched)
    {
        std::cerr << "Unknown QMOM failure case '" << selected_case << "'.\n";
        return 2;
    }

    return 0;
}
