#include "utilities.hpp"

#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <iostream>

namespace qmom_1d_continuous_breakage_mc_local_rates
{
const char* to_string(const qmomkit::ContinuousMomentRateStatus status)
{
    switch (status)
    {
    case qmomkit::ContinuousMomentRateStatus::success:
        return "success";
    case qmomkit::ContinuousMomentRateStatus::invalid_input:
        return "invalid_input";
    case qmomkit::ContinuousMomentRateStatus::inversion_failed:
        return "inversion_failed";
    case qmomkit::ContinuousMomentRateStatus::velocity_failed:
        return "velocity_failed";
    }
    return "unknown";
}

const char* to_string(const qmomkit::DiscontinuousMomentRateStatus status)
{
    switch (status)
    {
    case qmomkit::DiscontinuousMomentRateStatus::success:
        return "success";
    case qmomkit::DiscontinuousMomentRateStatus::invalid_input:
        return "invalid_input";
    case qmomkit::DiscontinuousMomentRateStatus::source_failed:
        return "source_failed";
    case qmomkit::DiscontinuousMomentRateStatus::invalid_rate:
        return "invalid_rate";
    }
    return "unknown";
}

const char* to_string(const qmomkit::QmomStatus status)
{
    switch (status)
    {
    case qmomkit::QmomStatus::success:
        return "success";
    case qmomkit::QmomStatus::invalid_input:
        return "invalid_input";
    case qmomkit::QmomStatus::not_realizable:
        return "not_realizable";
    case qmomkit::QmomStatus::no_adaptive_quadrature:
        return "no_adaptive_quadrature";
    }
    return "unknown";
}

const char* to_string(const qmomkit::SspRk3Status status)
{
    switch (status)
    {
    case qmomkit::SspRk3Status::success:
        return "success";
    case qmomkit::SspRk3Status::rejected_by_predicate:
        return "rejected_by_predicate";
    case qmomkit::SspRk3Status::minimum_timestep_reached:
        return "minimum_timestep_reached";
    }
    return "unknown";
}

void print_configuration_errors(const std::vector<qmomkit::config::ConfigurationError>& errors)
{
    for (const qmomkit::config::ConfigurationError& error : errors)
    {
        std::cerr << error.source_name << ':' << error.line << ": " << error.message << '\n';
    }
}

void print_common_run_summary(std::ostream& stream, const std::string_view workflow, const std::filesystem::path& configuration_file,
                              const std::filesystem::path& output_directory, const std::size_t sample_count,
                              const unsigned long long initialization_seed, const double initial_time, const double final_time,
                              const double time_step, const double output_interval)
{
    stream << "1D continuous + discontinuous source example with MC initialization\n";
    stream << "Workflow: " << workflow << '\n';
    stream << "Configuration: " << configuration_file << '\n';
    stream << "Output directory: " << output_directory << '\n';
    stream << "Samples: " << sample_count << ", initialization seed: " << initialization_seed << '\n';
    stream << "Continuous model: dd/dt = 1 - d\n";
    stream << "Binary fragmentation: d -> d/2\n";
    stream << "Time: initial = " << initial_time << ", final = " << final_time << ", step = " << time_step
           << ", output interval = " << output_interval << '\n';
}

void print_moments(std::ostream& stream, const std::string_view label, const std::vector<double>& moments)
{
    stream << '\n' << label << ":\n";
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        stream << "M" << order << " = " << moments[order] << '\n';
    }
}

OutputFiles open_output_files(const qmomkit::OutputOptions& output_options, const double initial_time, const std::vector<double>& initial_moments,
                              const qmomkit::Quadrature1D& initial_quadrature)
{
    OutputFiles output_files;
    if (output_options.write_moments || output_options.write_quadrature)
    {
        std::filesystem::create_directories(output_options.directory);
    }
    if (output_options.write_moments)
    {
        qmomkit::MomentWriter initial_moment_writer =
            qmomkit::MomentWriter::for_qmom(output_options.directory / "initial_moments.txt", initial_moments.size());
        initial_moment_writer.write(initial_time, qmomkit::MomentSet(initial_moments));

        output_files.moment_history_writer.emplace(
            qmomkit::MomentWriter::for_qmom(output_options.directory / "moment_history.txt", initial_moments.size()));
        output_files.moment_history_writer->write(initial_time, qmomkit::MomentSet(initial_moments));
    }
    if (output_options.write_quadrature)
    {
        qmomkit::QuadratureWriter initial_quadrature_writer =
            qmomkit::QuadratureWriter::for_qmom(output_options.directory / "initial_quadrature.txt", initial_quadrature.weights.size());
        initial_quadrature_writer.write(initial_time, initial_quadrature);

        output_files.quadrature_history_writer.emplace(
            qmomkit::QuadratureWriter::for_qmom(output_options.directory / "quadrature_history.txt", initial_quadrature.weights.size()));
        output_files.quadrature_history_writer->write(initial_time, initial_quadrature);
    }
    return output_files;
}

void write_output_state(OutputFiles& output_files, const double time, const std::vector<double>& moments, const qmomkit::Quadrature1D& quadrature)
{
    if (output_files.moment_history_writer.has_value())
    {
        output_files.moment_history_writer->write(time, qmomkit::MomentSet(moments));
    }
    if (output_files.quadrature_history_writer.has_value())
    {
        output_files.quadrature_history_writer->write(time, quadrature);
    }
}

void print_run_report(std::ostream& stream, const std::filesystem::path& configuration_file, const qmomkit::OutputOptions& output_options,
                      const std::size_t sample_count, const unsigned long long initialization_seed, const double initial_time,
                      const double final_time, const double time_step, const double output_interval, const std::vector<double>& initial_moments,
                      const qmomkit::Quadrature1D& initial_quadrature, const std::vector<double>& final_moments,
                      const qmomkit::Quadrature1D& final_quadrature, const double final_simulation_time)
{
    stream << std::setprecision(12);

    print_common_run_summary(stream, "QMOM public local moment-rate integration", configuration_file, output_options.directory, sample_count,
                             initialization_seed, initial_time, final_time, time_step, output_interval);
    stream << "Representation: sampled moments retained as the application-owned state vector\n";
    stream << "Continuous update: public moment-rate evaluator inside application-managed SSP-RK3\n";
    stream << "Source update: public discontinuous rate evaluator inside application-managed SSP-RK3\n";

    print_moments(stream, "Initial moments", initial_moments);

    stream << "\nInitial quadrature:\n";
    stream << "node  weight        abscissa\n";
    for (std::size_t node = 0; node < initial_quadrature.weights.size(); ++node)
    {
        stream << node << "     " << initial_quadrature.weights[node] << "   " << initial_quadrature.abscissas[node] << '\n';
    }

    print_moments(stream, "Final moments", final_moments);

    stream << "\nFinal quadrature at t = " << final_simulation_time << ":\n";
    stream << "node  initial_abscissa  final_abscissa\n";
    for (std::size_t node = 0; node < final_quadrature.abscissas.size(); ++node)
    {
        stream << node << "     " << initial_quadrature.abscissas[node] << "         " << final_quadrature.abscissas[node] << '\n';
    }
    if (output_options.write_moments || output_options.write_quadrature)
    {
        stream << "\nCompleted application-managed local-rate QMOM workflow. Outputs written to " << output_options.directory << '\n';
    }
    else
    {
        stream << "\nCompleted application-managed local-rate QMOM workflow. Output disabled by configuration.\n";
    }
}

} // namespace qmom_1d_continuous_breakage_mc_local_rates
