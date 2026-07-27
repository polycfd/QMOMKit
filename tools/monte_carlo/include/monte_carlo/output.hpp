#pragma once

#include "monte_carlo/moments.hpp"

#include <iosfwd>

namespace qmomkit::monte_carlo
{

/**
 * Write one-dimensional moment estimates as a QMOMKit-style text table.
 *
 * Metadata records the sample count, seed, and absolute number density used to
 * scale the sample averages. Data rows use whitespace-delimited columns:
 *
 *     alpha_0 value standard_error
 */
void write_moment_text(std::ostream& output, const OneDimensionalMomentEstimate& estimate, std::size_t sample_count, unsigned long long seed,
                       double number_density);

/**
 * Write tensor moment estimates as a QMOMKit-style text table.
 *
 * The `alpha_*` columns store the multi-index in phase-space dimension order.
 * The value and standard-error columns use the flat order of the supplied
 * `RawMomentTensor`.
 */
void write_moment_text(std::ostream& output, const TensorMomentEstimate& estimate, std::size_t sample_count, unsigned long long seed,
                       double number_density);

} // namespace qmomkit::monte_carlo
