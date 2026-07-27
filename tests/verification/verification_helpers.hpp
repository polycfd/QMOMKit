#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace qmomkit::verification
{

[[noreturn]] inline void fail_check(const char* expression, const char* file, const int line, const std::string& detail = {})
{
    std::cerr << "Verification failure at " << file << ':' << line << '\n';
    std::cerr << "  expression: " << expression << '\n';
    if (!detail.empty())
    {
        std::cerr << "  detail: " << detail << '\n';
    }
    std::exit(EXIT_FAILURE);
}

inline void check(const bool condition, const char* expression, const char* file, const int line)
{
    if (!condition)
    {
        fail_check(expression, file, line);
    }
}

inline void check_near(const double actual, const double expected, const double tolerance, const char* actual_expression,
                       const char* expected_expression, const char* tolerance_expression, const char* file, const int line)
{
    if (!std::isfinite(actual) || !std::isfinite(expected) || !std::isfinite(tolerance) || tolerance < 0.0 || actual > expected + tolerance ||
        actual < expected - tolerance)
    {
        std::ostringstream detail;
        detail << std::setprecision(17) << actual_expression << " = " << actual << ", " << expected_expression << " = " << expected << ", "
               << tolerance_expression << " = " << tolerance;
        fail_check("absolute near comparison", file, line, detail.str());
    }
}

inline void check_absolute_or_relative(const double actual, const double expected, const double absolute_tolerance, const double relative_tolerance,
                                       const char* actual_expression, const char* expected_expression, const char* absolute_tolerance_expression,
                                       const char* relative_tolerance_expression, const char* file, const int line)
{
    const double absolute_error = std::abs(actual - expected);
    const double denominator = std::max(std::abs(expected), 1.0);
    const double relative_error = absolute_error / denominator;
    if (!std::isfinite(actual) || !std::isfinite(expected) || !std::isfinite(absolute_tolerance) || !std::isfinite(relative_tolerance) ||
        absolute_tolerance < 0.0 || relative_tolerance < 0.0 || (absolute_error > absolute_tolerance && relative_error > relative_tolerance))
    {
        std::ostringstream detail;
        detail << std::setprecision(17) << actual_expression << " = " << actual << ", " << expected_expression << " = " << expected << ", "
               << "absolute_error = " << absolute_error << ", relative_error = " << relative_error << ", " << absolute_tolerance_expression << " = "
               << absolute_tolerance << ", " << relative_tolerance_expression << " = " << relative_tolerance;
        fail_check("absolute/relative near comparison", file, line, detail.str());
    }
}

} // namespace qmomkit::verification

#define QMOMKIT_VERIFY_CHECK(expression) ::qmomkit::verification::check(static_cast<bool>(expression), #expression, __FILE__, __LINE__)
#define QMOMKIT_VERIFY_NEAR(actual, expected, tolerance)                                                                                             \
    ::qmomkit::verification::check_near((actual), (expected), (tolerance), #actual, #expected, #tolerance, __FILE__, __LINE__)
#define QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(actual, expected, absolute_tolerance, relative_tolerance)                                                \
    ::qmomkit::verification::check_absolute_or_relative((actual), (expected), (absolute_tolerance), (relative_tolerance), #actual, #expected,        \
                                                        #absolute_tolerance, #relative_tolerance, __FILE__, __LINE__)
