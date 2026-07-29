#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace qmomkit::monte_carlo::test
{

[[noreturn]] inline void fail(const char* file, const int line, const std::string& message)
{
    std::cerr << file << ':' << line << ": " << message << '\n';
    std::exit(EXIT_FAILURE);
}

inline void require(const bool condition, const char* expression, const char* file, const int line, const std::string& context = {})
{
    if (!condition)
    {
        std::ostringstream message;
        message << "check failed: " << expression;
        if (!context.empty())
        {
            message << " [" << context << ']';
        }
        fail(file, line, message.str());
    }
}

inline void require_close(const double actual, const double expected, const double absolute_tolerance, const double relative_tolerance,
                          const char* quantity, const char* file, const int line, const std::string& context = {})
{
    const double absolute_error = std::abs(actual - expected);
    const double relative_scale = std::max({1.0, std::abs(actual), std::abs(expected)});
    const double allowed_error = absolute_tolerance + relative_tolerance * relative_scale;

    if (!std::isfinite(actual) || !std::isfinite(expected) || !std::isfinite(absolute_error) || absolute_error > allowed_error)
    {
        std::ostringstream message;
        message << std::setprecision(17) << "close check failed: " << quantity << " actual=" << actual << " expected=" << expected
                << " absolute_error=" << absolute_error << " allowed_error=" << allowed_error << " absolute_tolerance=" << absolute_tolerance
                << " relative_tolerance=" << relative_tolerance;
        if (!context.empty())
        {
            message << " [" << context << ']';
        }
        fail(file, line, message.str());
    }
}

inline void require_finite(const double value, const char* quantity, const char* file, const int line, const std::string& context = {})
{
    if (!std::isfinite(value))
    {
        std::ostringstream message;
        message << std::setprecision(17) << "finite check failed: " << quantity << " value=" << value;
        if (!context.empty())
        {
            message << " [" << context << ']';
        }
        fail(file, line, message.str());
    }
}

} // namespace qmomkit::monte_carlo::test

#define CHECK(expression) ::qmomkit::monte_carlo::test::require(static_cast<bool>(expression), #expression, __FILE__, __LINE__)

#define CHECK_CONTEXT(expression, context)                                                                                                           \
    ::qmomkit::monte_carlo::test::require(static_cast<bool>(expression), #expression, __FILE__, __LINE__, (context))

#define CHECK_CLOSE(actual, expected, absolute_tolerance)                                                                                            \
    ::qmomkit::monte_carlo::test::require_close((actual), (expected), (absolute_tolerance), 0.0, #actual, __FILE__, __LINE__)

#define CHECK_CLOSE_CONTEXT(actual, expected, absolute_tolerance, relative_tolerance, context)                                                       \
    ::qmomkit::monte_carlo::test::require_close((actual), (expected), (absolute_tolerance), (relative_tolerance), #actual, __FILE__, __LINE__,       \
                                                (context))

#define CHECK_FINITE(value) ::qmomkit::monte_carlo::test::require_finite((value), #value, __FILE__, __LINE__)
