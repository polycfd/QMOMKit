#include "qmomkit/realizability.hpp"

#include <iostream>

namespace
{

const char* status_name(const qmomkit::RealizabilityStatus status)
{
    switch (status)
    {
    case qmomkit::RealizabilityStatus::strictly_realizable:
        return "strictly_realizable";
    case qmomkit::RealizabilityStatus::boundary:
        return "boundary";
    case qmomkit::RealizabilityStatus::not_realizable:
        return "not_realizable";
    case qmomkit::RealizabilityStatus::invalid_input:
        return "invalid_input";
    case qmomkit::RealizabilityStatus::numerically_indeterminate:
        return "numerically_indeterminate";
    }
    return "unknown";
}

} // namespace

int main()
{
    // Raw moments m0..m3 for an exponential distribution with density 1.
    const qmomkit::MomentSet moments{1.0, 1.0, 2.0, 6.0};

    const qmomkit::RealizabilityResult result = qmomkit::check_realizability(moments, qmomkit::MomentSupport::nonnegative);

    std::cout << "Explicit realizability check\n";
    std::cout << "support = nonnegative\n";
    std::cout << "status = " << status_name(result.status) << '\n';
    std::cout << "message = " << result.message << '\n';

    return result.realizable() ? 0 : 1;
}
