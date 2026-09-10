#pragma once

#include "monte_carlo/moments.hpp"
#include "qmomkit/output.hpp"

#include <filesystem>
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
 * Create generalized-output-capable history writers from the shared output
 * policy. These helpers keep generalized column selection out of applications.
 */
[[nodiscard]] qmomkit::MomentWriter make_moment_history_writer(const std::filesystem::path& path, std::size_t moment_count,
                                                               const qmomkit::OutputOptions& output_options);
[[nodiscard]] qmomkit::MomentWriter make_moment_history_writer(const std::filesystem::path& path, std::vector<std::size_t> extents,
                                                               const qmomkit::OutputOptions& output_options);
[[nodiscard]] qmomkit::MomentWriter make_moment_history_writer(const std::filesystem::path& path, std::vector<std::size_t> extents,
                                                               std::vector<qmomkit::MultiIndex> selected_moments,
                                                               const qmomkit::OutputOptions& output_options);
[[nodiscard]] qmomkit::MomentWriter make_moment_history_writer(const std::filesystem::path& path, const qmomkit::CqmomMomentLayout& layout,
                                                               const qmomkit::OutputOptions& output_options);

/** Append standard moments and any configured generalized estimates from the current ensemble. */
void write_moment_history_row(qmomkit::MomentWriter& writer, double time, const qmomkit::MomentSet& moments, const Ensemble& ensemble,
                              double number_density, const qmomkit::OutputOptions& output_options);
void write_moment_history_row(qmomkit::MomentWriter& writer, double time, const qmomkit::RawMomentTensor& moments, const Ensemble& ensemble,
                              double number_density, const qmomkit::OutputOptions& output_options);
void write_moment_history_row(qmomkit::MomentWriter& writer, double time, const qmomkit::CqmomMomentData& moments, const Ensemble& ensemble,
                              double number_density, const qmomkit::OutputOptions& output_options);

} // namespace qmomkit::monte_carlo
