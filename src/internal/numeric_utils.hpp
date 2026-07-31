#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace qmomkit::internal
{

[[nodiscard]] inline bool all_finite(const std::vector<double>& values)
{
    return std::all_of(values.begin(), values.end(), [](const double value) { return std::isfinite(value); });
}

[[nodiscard]] inline std::string format_diagnostic_vector(const std::vector<double>& values)
{
    std::ostringstream output;
    output << std::setprecision(std::numeric_limits<double>::max_digits10) << '[';
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (i > 0)
        {
            output << ", ";
        }
        output << values[i];
    }
    output << ']';
    return output.str();
}

[[nodiscard]] inline std::string phase_space_velocity_wrong_dimension_message(const std::size_t expected_size, const std::size_t returned_size)
{
    std::ostringstream output;
    output << "Phase-space velocity callback returned the wrong dimension:\n"
           << "  expected size: " << expected_size << '\n'
           << "  returned size: " << returned_size;
    return output.str();
}

[[nodiscard]] inline std::string phase_space_velocity_non_finite_message(const double time, const std::size_t dimension, const double value,
                                                                         const std::vector<double>& state, const std::vector<double>& velocity)
{
    std::ostringstream output;
    output << std::setprecision(std::numeric_limits<double>::max_digits10) << "Phase-space velocity callback returned a non-finite value:\n"
           << "  time: " << time << '\n'
           << "  dimension: " << dimension << '\n'
           << "  value: " << value << '\n'
           << "  state: " << format_diagnostic_vector(state) << '\n'
           << "  velocity: " << format_diagnostic_vector(velocity);
    return output.str();
}

inline void validate_phase_space_velocity_result(const std::vector<double>& velocity, const std::size_t expected_size, const double time,
                                                 const std::vector<double>& state)
{
    if (velocity.size() != expected_size)
    {
        throw std::invalid_argument(phase_space_velocity_wrong_dimension_message(expected_size, velocity.size()));
    }
    for (std::size_t dimension = 0; dimension < velocity.size(); ++dimension)
    {
        if (!std::isfinite(velocity[dimension]))
        {
            throw std::invalid_argument(phase_space_velocity_non_finite_message(time, dimension, velocity[dimension], state, velocity));
        }
    }
}

} // namespace qmomkit::internal
