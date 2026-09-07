#pragma once

#include <initializer_list>
#include <iostream>
#include <string_view>

namespace qmomkit::test
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

} // namespace qmomkit::test
