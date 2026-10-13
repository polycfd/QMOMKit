#pragma once

#include "qmomkit/distributions.hpp"

#include <iosfwd>
#include <string_view>
#include <vector>

namespace cqmom_5d_non_spherical_particles
{

/**
 * The moments both executables print: the number density, the first raw
 * moments, the second raw moment of each coordinate, and the diameter-shape
 * moment M[1,1,0,0,0] that breakage couples.
 */
[[nodiscard]] qmomkit::RawMomentRequest printed_moments();

/**
 * Prints the initial and final values of `printed_moments()` side by side. Both
 * executables use it, so their summaries line up for comparison.
 */
void print_moments(std::ostream& stream, std::string_view run, const std::vector<double>& initial_values, const std::vector<double>& final_values);

} // namespace cqmom_5d_non_spherical_particles
