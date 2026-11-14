#pragma once

#include <iosfwd>
#include <string_view>
#include <vector>

namespace gqmom_1d_droplet_breakup_with_mc
{

/**
 * Prints the initial and final raw moments side by side. Both executables use
 * it, so their summaries line up for comparison.
 */
void print_moments(std::ostream& stream, std::string_view run, const std::vector<double>& initial_moments, const std::vector<double>& final_moments);

} // namespace gqmom_1d_droplet_breakup_with_mc
