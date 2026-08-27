#include "utilities.hpp"

#include "qmomkit/multi_index.hpp"

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <system_error>
#include <utility>
#include <vector>

namespace paper_3d_droplets_condensation_nucleation
{
namespace
{
constexpr std::size_t max_console_moments = 12;

void print_multi_index(const qmomkit::MultiIndex& order)
{
    std::cout << '[';
    for (std::size_t dimension = 0; dimension < order.size(); ++dimension)
    {
        if (dimension > 0)
        {
            std::cout << ',';
        }
        std::cout << order[dimension];
    }
    std::cout << ']';
}

void print_path(const qmomkit::ConditioningPath& path)
{
    std::cout << '[';
    for (std::size_t level = 0; level < path.size(); ++level)
    {
        if (level > 0)
        {
            std::cout << ',';
        }
        std::cout << path.node_index_at_level(level);
    }
    std::cout << ']';
}

void print_selected_moments(const char* label, const qmomkit::SelectedRawMomentResult& moments)
{
    std::cout << label << '\n';
    const std::size_t printed_count = std::min(max_console_moments, moments.values.size());
    for (std::size_t index = 0; index < printed_count; ++index)
    {
        std::cout << "M";
        print_multi_index(moments.indices[index]);
        std::cout << " = " << moments.values[index] << '\n';
    }
    if (moments.values.size() > printed_count)
    {
        std::cout << "... " << (moments.values.size() - printed_count) << " more moments omitted from the console summary\n";
    }
}

void print_index_vector(const std::vector<std::size_t>& values)
{
    std::cout << '[';
    for (std::size_t index = 0; index < values.size(); ++index)
    {
        if (index > 0)
        {
            std::cout << ',';
        }
        std::cout << values[index];
    }
    std::cout << ']';
}

void print_quadrature(const qmomkit::CqmomQuadrature& quadrature)
{
    std::cout << "node  path       weight          r               w               T\n";
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        std::cout << node << "     ";
        print_path(quadrature.terminal_paths()[node]);
        std::cout << "   " << quadrature.weights()[node] << "   " << quadrature.abscissas()[0][node] << "   " << quadrature.abscissas()[1][node]
                  << "   " << quadrature.abscissas()[2][node] << '\n';
    }
}

void print_run_metadata(const std::size_t sample_count, const unsigned long long initialization_seed, const double time_step, const double final_time)
{
    std::cout << "sample_count = " << sample_count << '\n';
    std::cout << "initialization_seed = " << initialization_seed << '\n';
    std::cout << "time_step = " << time_step << '\n';
    std::cout << "final_time = " << final_time << '\n';
}

} // namespace

void print_configuration_errors(const std::string& error)
{
    std::cerr << error << '\n';
}

void ensure_output_directory(const std::filesystem::path& path)
{
    std::error_code error;
    std::filesystem::create_directories(path, error);
    if (error)
    {
        throw std::runtime_error("Could not create output directory " + path.string() + ": " + error.message());
    }
}

void print_cqmom_run_configuration(const qmomkit::CqmomCase& cqmom_case,
                                   const qmomkit::CqmomPhaseSpaceEvolutionOptions& phase_space_evolution_options)
{
    std::cout << std::setprecision(12);
    std::cout << "phase-space dimensions: 3\n";
    std::cout << "phase-space ordering: (r,w,T)\n";
    std::cout << "conditioning_order = ";
    print_index_vector(cqmom_case.hierarchy().conditioning_order());
    std::cout << '\n';
    std::cout << "node_counts = ";
    print_index_vector(cqmom_case.hierarchy().node_counts());
    std::cout << '\n';
    std::cout << "phase-space evolution: Strang-split advection and boundary contribution\n";
    std::cout << qmomkit::format_cqmom_settings(cqmom_case.moment_storage(), cqmom_case.adaptive_branch_counts(), cqmom_case.vandermonde_solver(),
                                                cqmom_case.diagnostic_level());
    std::cout << qmomkit::format_node_advection_options(phase_space_evolution_options.node_options);
}

void print_monte_carlo_run_configuration(const qmomkit::CqmomCase& cqmom_case, const qmomkit::NodeAdvectionOptions& node_options)
{
    std::cout << std::setprecision(12);
    std::cout << "phase-space dimensions: 3\n";
    std::cout << "phase-space ordering: (r,w,T)\n";
    std::cout << "node_counts = ";
    print_index_vector(cqmom_case.hierarchy().node_counts());
    std::cout << '\n';
    std::cout << "particle evolution: Strang-split advection and population injection\n";
    std::cout << qmomkit::format_node_advection_options(node_options);
}

void print_result_summary(const qmomkit::SelectedRawMomentResult& initial_moments, const qmomkit::SelectedRawMomentResult& final_moments,
                          const qmomkit::CqmomQuadrature& final_quadrature, const std::filesystem::path& output_directory,
                          const std::size_t sample_count, const unsigned long long initialization_seed, const double time_step,
                          const double final_time, const bool wrote_moments, const bool wrote_quadrature)
{
    std::cout << std::setprecision(12);
    std::cout << "3D CQMOM condensing-droplet advection evolution example\n";
    std::cout << "phase-space ordering: (r,w,T)\n";
    print_run_metadata(sample_count, initialization_seed, time_step, final_time);

    std::cout << '\n';
    print_selected_moments("Initial moments:", initial_moments);

    std::cout << '\n';
    print_selected_moments("Final moments:", final_moments);

    std::cout << "\nFinal quadrature:\n";
    print_quadrature(final_quadrature);

    if (wrote_moments || wrote_quadrature)
    {
        std::cout << "\nOutputs written to:\n";
        if (wrote_moments)
        {
            std::cout << (output_directory / "initial_moments.txt").string() << '\n';
            std::cout << (output_directory / "moment_history.txt").string() << '\n';
        }
        if (wrote_quadrature)
        {
            std::cout << (output_directory / "quadrature_history.txt").string() << '\n';
        }
    }
    else
    {
        std::cout << "\nOutput disabled by configuration.\n";
    }
}

void print_monte_carlo_result_summary(const qmomkit::SelectedRawMomentResult& initial_moments, const qmomkit::SelectedRawMomentResult& final_moments,
                                      const std::filesystem::path& output_directory, const std::size_t final_particle_count,
                                      const std::size_t sample_count, const unsigned long long initialization_seed, const double time_step,
                                      const double final_time, const bool wrote_moments)
{
    std::cout << std::setprecision(12);
    std::cout << "3D direct Monte Carlo condensing-droplet advection evolution example\n";
    std::cout << "phase-space ordering: (r,w,T)\n";
    print_run_metadata(sample_count, initialization_seed, time_step, final_time);

    std::cout << '\n';
    print_selected_moments("Initial moments:", initial_moments);

    std::cout << '\n';
    print_selected_moments("Final moments:", final_moments);

    std::cout << "\nFinal ensemble:\n";
    std::cout << "final_sample_count = " << final_particle_count << '\n';

    if (wrote_moments)
    {
        std::cout << "\nOutputs written to:\n";
        std::cout << (output_directory / "moment_history_mc.txt").string() << '\n';
    }
    else
    {
        std::cout << "\nMoment output disabled by configuration.\n";
    }
}

} // namespace paper_3d_droplets_condensation_nucleation
