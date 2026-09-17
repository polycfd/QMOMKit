#pragma once

#include <initializer_list>
#include <iostream>
#include <string_view>

// Selector dispatch for QMOMKitMonteCarlo test executables, so each logical
// case can be registered with CTest on its own instead of being hidden inside a
// whole-file run. The companion package is a standalone consumer and must not
// reach into the root test tree, so this mirrors the root helper rather than
// including it.
//
// Running the executable with no argument runs every case, which keeps a direct
// invocation useful while developing.

namespace qmomkit::monte_carlo::test
{

struct TestCase
{
    std::string_view selector;
    void (*run)();
};

inline int dispatch_test_case(const int argc, char** argv, const std::initializer_list<TestCase> cases)
{
    if (argc == 1)
    {
        for (const TestCase& test_case : cases)
        {
            test_case.run();
        }
        return 0;
    }

    if (argc == 2)
    {
        const std::string_view selector(argv[1]);
        for (const TestCase& test_case : cases)
        {
            if (test_case.selector == selector)
            {
                test_case.run();
                return 0;
            }
        }
    }

    std::cerr << "Unknown test selector. Available selectors:\n";
    for (const TestCase& test_case : cases)
    {
        std::cerr << "  " << test_case.selector << '\n';
    }
    return 2;
}

} // namespace qmomkit::monte_carlo::test
