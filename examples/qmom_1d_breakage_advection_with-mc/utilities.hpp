#pragma once

#include "qmomkit/configuration.hpp"

#include <cstddef>
#include <filesystem>
#include <iosfwd>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace qmom_1d_breakage_advection_with_mc
{

void print_configuration_errors(const std::string& error);

void print_common_run_summary(std::ostream& stream, std::string_view workflow, const std::filesystem::path& configuration_file,
                              const std::filesystem::path& output_directory, std::size_t sample_count, unsigned long long initialization_seed,
                              double initial_time, double final_time, double time_step, double output_interval);

void print_moments(std::ostream& stream, std::string_view label, const std::vector<double>& moments);

void write_particle_history(const std::filesystem::path& path, const std::vector<std::tuple<double, std::size_t, double>>& history);

} // namespace qmom_1d_breakage_advection_with_mc
