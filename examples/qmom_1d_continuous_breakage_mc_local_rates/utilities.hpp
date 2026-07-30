#pragma once

#include "qmomkit/configuration.hpp"
#include "qmomkit/moment_rates.hpp"
#include "qmomkit/output.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/quadrature.hpp"
#include "qmomkit/runge_kutta.hpp"

#include <cstddef>
#include <filesystem>
#include <iosfwd>
#include <optional>
#include <string_view>
#include <vector>

namespace qmom_1d_continuous_breakage_mc_local_rates
{

struct OutputFiles
{
    std::optional<qmomkit::MomentWriter> moment_history_writer;
    std::optional<qmomkit::QuadratureWriter> quadrature_history_writer;
};

[[nodiscard]] const char* to_string(qmomkit::ContinuousMomentRateStatus status);
[[nodiscard]] const char* to_string(qmomkit::DiscontinuousMomentRateStatus status);
[[nodiscard]] const char* to_string(qmomkit::QmomStatus status);
[[nodiscard]] const char* to_string(qmomkit::SspRk3Status status);

void print_common_run_summary(std::ostream& stream, std::string_view workflow, const std::filesystem::path& configuration_file,
                              const std::filesystem::path& output_directory, std::size_t sample_count, unsigned long long initialization_seed,
                              double initial_time, double final_time, double time_step, double output_interval);

void print_moments(std::ostream& stream, std::string_view label, const std::vector<double>& moments);

[[nodiscard]] OutputFiles open_output_files(const qmomkit::OutputOptions& output_options, double initial_time,
                                            const std::vector<double>& initial_moments, const qmomkit::Quadrature1D& initial_quadrature);

void write_output_state(OutputFiles& output_files, double time, const std::vector<double>& moments, const qmomkit::Quadrature1D& quadrature);

void print_run_report(std::ostream& stream, const std::filesystem::path& configuration_file, const qmomkit::OutputOptions& output_options,
                      std::size_t sample_count, unsigned long long initialization_seed, double initial_time, double final_time, double time_step,
                      double output_interval, const std::vector<double>& initial_moments, const qmomkit::Quadrature1D& initial_quadrature,
                      const std::vector<double>& final_moments, const qmomkit::Quadrature1D& final_quadrature, double final_simulation_time);

} // namespace qmom_1d_continuous_breakage_mc_local_rates
