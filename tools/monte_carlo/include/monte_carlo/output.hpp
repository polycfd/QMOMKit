#pragma once

#include "monte_carlo/moments.hpp"
#include "qmomkit/output.hpp"

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

/**
 * Append standard tensor moments and configured generalized moments to a
 * QMOMKit moment-history writer. Generalized values use the same equally
 * weighted sample convention as the standard Monte Carlo estimators.
 */
void write_moment_history_row(qmomkit::MomentWriter& writer, double time, const qmomkit::RawMomentTensor& moments, const Ensemble& ensemble,
                              double number_density, const std::vector<std::vector<double>>& generalized_moment_orders);

/** Append layout-ordered standard moments and configured generalized moments. */
void write_moment_history_row(qmomkit::MomentWriter& writer, double time, const qmomkit::CqmomMomentData& moments, const Ensemble& ensemble,
                              double number_density, const std::vector<std::vector<double>>& generalized_moment_orders);

} // namespace qmomkit::monte_carlo
