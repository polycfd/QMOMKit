#pragma once

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

// Packed-callback counterpart. Diagnostics stay per node: the offending entry
// is decomposed back into its node and phase-space dimension so the reported
// coordinate and velocity slices match the per-node message exactly. The
// slices are built only on the failure path.
inline void validate_packed_phase_space_velocity_result(const std::vector<double>& velocities, const std::size_t expected_size, const double time,
                                                        const std::vector<double>& coordinates, const std::size_t dimensions)
{
    if (velocities.size() != expected_size)
    {
        throw std::invalid_argument(phase_space_velocity_wrong_dimension_message(expected_size, velocities.size()));
    }
    for (std::size_t entry = 0; entry < velocities.size(); ++entry)
    {
        if (!std::isfinite(velocities[entry]))
        {
            const std::size_t node = entry / dimensions;
            const std::size_t dimension = entry % dimensions;
            const auto node_begin = static_cast<std::ptrdiff_t>(node * dimensions);
            const auto node_end = static_cast<std::ptrdiff_t>((node + 1) * dimensions);
            const std::vector<double> node_coordinate(coordinates.begin() + node_begin, coordinates.begin() + node_end);
            const std::vector<double> node_velocity(velocities.begin() + node_begin, velocities.begin() + node_end);
            throw std::invalid_argument(phase_space_velocity_non_finite_message(time, dimension, velocities[entry], node_coordinate, node_velocity));
        }
    }
}

} // namespace qmomkit::internal
