#pragma once

#include "qmomkit/cqmom.hpp"
#include "qmomkit/distributions.hpp"
#include "qmomkit/multi_index.hpp"

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace cqmom_reactive_spherical_particles
{

void print_configuration_errors(const std::string& error);

void ensure_output_directory(const std::filesystem::path& path);

[[nodiscard]] std::vector<qmomkit::MultiIndex> select_terminal_moment_indices(const std::vector<qmomkit::MultiIndex>& required_indices);

void print_result_summary(const qmomkit::SelectedRawMomentResult& initial_moments, const qmomkit::SelectedRawMomentResult& final_moments,
                          const qmomkit::CqmomQuadrature& final_quadrature, const std::filesystem::path& output_directory, std::size_t sample_count,
                          unsigned long long initialization_seed, double time_step, double final_time, bool wrote_moments, bool wrote_quadrature);

void print_monte_carlo_result_summary(const qmomkit::RawMomentTensor& initial_moments, const qmomkit::RawMomentTensor& final_moments,
                                      const std::vector<qmomkit::MultiIndex>& terminal_moment_indices, const std::filesystem::path& output_directory,
                                      std::size_t final_particle_count, std::size_t sample_count, unsigned long long initialization_seed,
                                      double time_step, double final_time, bool wrote_moments);

} // namespace cqmom_reactive_spherical_particles
