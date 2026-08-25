#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/distributions.hpp"

#include <filesystem>
#include <vector>

namespace cqmom_5d_spherical_breakup_moments
{

[[nodiscard]] std::vector<std::size_t> maximum_orders(const std::vector<std::size_t>& extents);
[[nodiscard]] qmomkit::RawMomentTensor ensemble_moments(const qmomkit::monte_carlo::Ensemble& ensemble, double particle_weight,
                                                        const std::vector<std::size_t>& extents);
[[nodiscard]] double maximum_tensor_difference(const qmomkit::RawMomentTensor& left, const qmomkit::RawMomentTensor& right);
[[nodiscard]] qmomkit::SelectedRawMomentResult select_moments(const qmomkit::RawMomentTensor& moments, const qmomkit::RawMomentRequest& request);
void print_selected_moments(const char* label, const qmomkit::SelectedRawMomentResult& moments);
void ensure_output_directory(const std::filesystem::path& path);

} // namespace cqmom_5d_spherical_breakup_moments
