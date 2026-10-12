#pragma once

#include "monte_carlo/moments.hpp"
#include "qmomkit/output.hpp"
#include "qmomkit/raw_moment_tensor.hpp"

#include <cstddef>
#include <filesystem>
#include <vector>

namespace qmomkit::monte_carlo
{

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
