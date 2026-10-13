#pragma once

#include "qmomkit/distributions.hpp"

#include <iosfwd>
#include <string_view>
#include <vector>

namespace cqmom_2d_bubbles_cavitation
{

/**
 * The moments both executables print: the number density, the first and the
 * second raw moments of (R, U).
 */
[[nodiscard]] qmomkit::RawMomentRequest printed_moments();

/**
 * Prints the initial and final values of `printed_moments()` side by side. Both
 * executables use it, so their summaries line up for comparison.
 */
void print_moments(std::ostream& stream, std::string_view run, const std::vector<double>& initial_values, const std::vector<double>& final_values);

} // namespace cqmom_2d_bubbles_cavitation
