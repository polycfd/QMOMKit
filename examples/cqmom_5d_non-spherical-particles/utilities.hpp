#pragma once

#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/multi_index.hpp"
#include "qmomkit/raw_moment_tensor.hpp"

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace cqmom_5d_non_spherical_particles
{

void print_configuration_errors(const std::string& error);

void ensure_output_directory(const std::filesystem::path& path);

[[nodiscard]] std::vector<qmomkit::MultiIndex> select_terminal_moment_indices(const qmomkit::MultiIndexTable& required_indices);

void print_result_summary(const qmomkit::SelectedRawMomentResult& initial_moments, const qmomkit::SelectedRawMomentResult& final_moments,
                          const std::filesystem::path& output_directory, bool wrote_moments, bool wrote_quadrature);

void print_monte_carlo_result_summary(const qmomkit::RawMomentTensor& initial_moments, const qmomkit::RawMomentTensor& final_moments,
                                      const std::vector<qmomkit::MultiIndex>& terminal_moment_indices, const std::filesystem::path& output_directory,
                                      std::size_t final_particle_count, bool wrote_moments);

} // namespace cqmom_5d_non_spherical_particles
