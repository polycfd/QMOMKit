#include "test_case_dispatch.hpp"

#include <cassert>
#include <qmomkit/abscissa_bounds.hpp>

namespace
{

void test_default_construction_is_zero_initialized()
{
    // AbscissaBounds has no invariant tying lower/upper together, so a
    // default-constructed instance is not expected to represent a valid
    // interval; this only guards against reading uninitialized memory.
    const qmomkit::AbscissaBounds bounds;
    assert(bounds.lower == 0.0);
    assert(bounds.upper == 0.0);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"default_construction_is_zero_initialized", test_default_construction_is_zero_initialized},
                                             });
}
