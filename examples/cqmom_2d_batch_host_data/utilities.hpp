#pragma once

#include "host_data.hpp"
#include "qmomkit/cqmom.hpp"

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace cqmom_2d_batch_host_data
{

void print_configuration_errors(const std::string& error);

void ensure_output_directory(const std::filesystem::path& path);

// Output file of one location, for example output/location_3_moment_history.txt.
[[nodiscard]] std::filesystem::path location_output_file(const std::filesystem::path& directory, std::size_t location, const std::string& name);

// Prints one row per location: the computed number density and mean particle
// velocity, each next to the exact solution of the model.
void print_location_summary(const LocalHostData& host_data, const std::vector<qmomkit::CqmomState>& states, double initial_time,
                            double initial_number_density, double initial_mean_velocity, double time);

} // namespace cqmom_2d_batch_host_data
