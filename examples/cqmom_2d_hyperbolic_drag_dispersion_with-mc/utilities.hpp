#pragma once

#include "qmomkit/distributions.hpp"

#include <iosfwd>
#include <string_view>
#include <vector>

namespace cqmom_2d_hyperbolic_drag_dispersion
{

/**
 * The moments all three executables print: the number density, the first and
 * second raw moments of (u, d), and M[2,1], which the hierarchy holds as a mixed
 * moment of parent order 2 when the hyperbolic velocity is the parent.
 */
[[nodiscard]] qmomkit::RawMomentRequest printed_moments();

/**
 * Prints the initial and final values of `printed_moments()` side by side. All
 * three executables use it, so their summaries line up for comparison.
 */
void print_moments(std::ostream& stream, std::string_view run, const std::vector<double>& initial_values, const std::vector<double>& final_values);

} // namespace cqmom_2d_hyperbolic_drag_dispersion
