#include "utilities.hpp"

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace qmom_1d_continuous_breakage_mc_initialization
{
namespace
{

std::ofstream open_output_file(const std::filesystem::path& path)
{
    std::ofstream file(path);
    if (!file)
    {
        throw std::runtime_error("Could not write " + path.string());
    }
    file << std::setprecision(17);
    return file;
}

} // namespace

void print_configuration_errors(const std::string& error)
{
    std::cerr << error << '\n';
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

void write_particle_history(const std::filesystem::path& path, const std::vector<std::tuple<double, std::size_t, double>>& history)
{
    std::ofstream file = open_output_file(path);
    file << "# time particle_id diameter\n";
    for (const auto& [time, particle_id, diameter] : history)
    {
        file << std::setprecision(12) << time << std::setprecision(17) << ' ' << particle_id << ' ' << diameter << '\n';
    }
}

} // namespace qmom_1d_continuous_breakage_mc_initialization
