#include "qmomkit/realizability.hpp"

#include <iostream>

int main()
{
    // Raw moments m0..m3 for an exponential distribution with density 1.
    const qmomkit::MomentSet moments{1.0, 1.0, 2.0, 6.0};

    const qmomkit::RealizabilityResult result = qmomkit::check_realizability(moments, qmomkit::MomentSupport::nonnegative);

    std::cout << "Explicit realizability check\n";
    std::cout << "support = nonnegative\n";
    std::cout << "status = ";
    switch (result.status)
    {
    case qmomkit::RealizabilityStatus::strictly_realizable:
        std::cout << "strictly_realizable";
        break;
    case qmomkit::RealizabilityStatus::boundary:
        std::cout << "boundary";
        break;
    case qmomkit::RealizabilityStatus::not_realizable:
        std::cout << "not_realizable";
        break;
    case qmomkit::RealizabilityStatus::invalid_input:
        std::cout << "invalid_input";
        break;
    case qmomkit::RealizabilityStatus::numerically_indeterminate:
        std::cout << "numerically_indeterminate";
        break;
    }
    std::cout << '\n';
    std::cout << "message = " << result.message << '\n';

    return result.realizable() ? 0 : 1;
}
