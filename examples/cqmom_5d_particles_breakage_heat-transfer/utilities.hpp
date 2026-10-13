#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/distributions.hpp"

#include <cstddef>
#include <iosfwd>
#include <string_view>
#include <vector>

namespace cqmom_5d_particles_breakage_heat_transfer
{

[[nodiscard]] std::vector<std::size_t> maximum_orders(const std::vector<std::size_t>& extents);

/**
 * The raw moments of the ensemble within the given extents, with each sample
 * carrying `particle_weight`. Throws when the estimate fails.
 */
[[nodiscard]] qmomkit::RawMomentTensor ensemble_moments(const qmomkit::monte_carlo::Ensemble& ensemble, double particle_weight,
                                                        const std::vector<std::size_t>& extents);

/**
 * The moments both executables print: the number density, the first raw
 * moments, the second raw moments of each coordinate, and the
 * volume-temperature moment M[1,1,0,0,0] that size-dependent heating couples.
 */
[[nodiscard]] qmomkit::RawMomentRequest printed_moments();

/**
 * Prints the initial and final values of `printed_moments()` side by side. Both
 * executables use it, so their summaries line up for comparison.
 */
void print_moments(std::ostream& stream, std::string_view run, const std::vector<double>& initial_values, const std::vector<double>& final_values);

} // namespace cqmom_5d_particles_breakage_heat_transfer
