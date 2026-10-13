#pragma once

#include <cstddef>
#include <filesystem>
#include <iosfwd>
#include <string_view>
#include <tuple>
#include <vector>

namespace qmom_1d_breakage_advection_with_mc
{

/**
 * Prints the initial and final raw moments side by side. Both executables use
 * it, so their summaries line up for comparison.
 */
void print_moments(std::ostream& stream, std::string_view run, const std::vector<double>& initial_moments, const std::vector<double>& final_moments);

void write_particle_history(const std::filesystem::path& path, const std::vector<std::tuple<double, std::size_t, double>>& history);

} // namespace qmom_1d_breakage_advection_with_mc
