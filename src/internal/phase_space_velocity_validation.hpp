#pragma once

#include "numeric_utils.hpp"

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

[[nodiscard]] inline std::vector<double> phase_space_node_slice(const std::vector<double>& values, const std::size_t node,
                                                                const std::size_t dimensions)
{
    const auto node_begin = static_cast<std::ptrdiff_t>(node * dimensions);
    return std::vector<double>(values.begin() + node_begin, values.begin() + node_begin + static_cast<std::ptrdiff_t>(dimensions));
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

// The packed callback fills a buffer the library supplies rather than
// returning one, so a callback that writes only part of it would otherwise
// advance the remaining nodes with whatever the library happened to leave
// there. The unwritten-entry marking itself is shared with the other
// fill-style callback boundary in `numeric_utils.hpp`; only the diagnostics
// below are specific to velocities.
#ifndef NDEBUG
[[nodiscard]] inline std::string phase_space_velocity_unfilled_message(const double time, const std::size_t node, const std::size_t dimension,
                                                                       const std::vector<double>& coordinate)
{
    std::ostringstream output;
    output << std::setprecision(std::numeric_limits<double>::max_digits10) << "Phase-space velocity callback did not fill every velocity entry:\n"
           << "  time: " << time << '\n'
           << "  node: " << node << '\n'
           << "  dimension: " << dimension << '\n'
           << "  coordinate: " << format_diagnostic_vector(coordinate);
    return output.str();
}
#endif

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
#ifndef NDEBUG
    // Ahead of the finiteness scan so that an unwritten entry is reported as
    // the missing write it is, rather than as a non-finite value the callback
    // never produced.
    for (std::size_t entry = 0; entry < velocities.size(); ++entry)
    {
        if (is_output_buffer_poison(velocities[entry]))
        {
            const std::size_t node = entry / dimensions;
            throw std::invalid_argument(
                phase_space_velocity_unfilled_message(time, node, entry % dimensions, phase_space_node_slice(coordinates, node, dimensions)));
        }
    }
#endif
    for (std::size_t entry = 0; entry < velocities.size(); ++entry)
    {
        if (!std::isfinite(velocities[entry]))
        {
            const std::size_t node = entry / dimensions;
            const std::vector<double> node_coordinate = phase_space_node_slice(coordinates, node, dimensions);
            const std::vector<double> node_velocity = phase_space_node_slice(velocities, node, dimensions);
            throw std::invalid_argument(
                phase_space_velocity_non_finite_message(time, entry % dimensions, velocities[entry], node_coordinate, node_velocity));
        }
    }
}

} // namespace qmomkit::internal
