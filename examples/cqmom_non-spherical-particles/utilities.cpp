#include "utilities.hpp"

#include "qmomkit/multi_index.hpp"

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <system_error>
#include <vector>

namespace cqmom_non_spherical_particles
{
namespace
{
constexpr std::size_t max_console_moments = 12;

std::size_t moment_total_degree(const qmomkit::MultiIndex& order)
{
    std::size_t degree = 0;
    for (const std::size_t exponent : order)
    {
        degree += exponent;
    }
    return degree;
}

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

void print_moment_entry(const qmomkit::MultiIndex& order, const double value)
{
    std::cout << "M";
    print_multi_index(order);
    std::cout << " = " << value << '\n';
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
        print_moment_entry(moments.indices[index], moments.values[index]);
    }
    if (moments.values.size() > printed_count)
    {
        std::cout << "... " << (moments.values.size() - printed_count) << " more moments omitted from the console summary\n";
    }
}

void print_moment_tensor(const qmomkit::RawMomentTensor& moments, const std::vector<qmomkit::MultiIndex>& terminal_moment_indices)
{
    if (moments.values.size() != qmomkit::product_extents(moments.extents))
    {
        throw std::runtime_error("The Monte Carlo moment tensor does not match its configured extents.");
    }
    for (const qmomkit::MultiIndex& order : terminal_moment_indices)
    {
        if (order.size() != moments.extents.size())
        {
            throw std::runtime_error("A terminal moment index does not match the Monte Carlo moment tensor rank.");
        }
        for (std::size_t dimension = 0; dimension < order.size(); ++dimension)
        {
            if (order[dimension] >= moments.extents[dimension])
            {
                throw std::runtime_error("A terminal moment index is outside the Monte Carlo moment tensor extents.");
            }
        }
        const std::size_t flat_index = qmomkit::flatten_multi_index(order, moments.extents);
        print_moment_entry(order, moments.values[flat_index]);
    }
}

void print_quadrature(const qmomkit::CqmomQuadrature& quadrature)
{
    std::cout << "node  path            weight          d               chi             u               v               w\n";
    if (quadrature.number_of_dimensions() != 5)
    {
        throw std::runtime_error("The final quadrature is not ordered as a 5D non-spherical particle state.");
    }
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        std::cout << node << "     ";
        print_path(quadrature.terminal_paths()[node]);
        std::cout << "   " << quadrature.weights()[node] << "   " << quadrature.abscissas()[0][node] << "   " << quadrature.abscissas()[1][node]
                  << "   " << quadrature.abscissas()[2][node] << "   " << quadrature.abscissas()[3][node] << "   " << quadrature.abscissas()[4][node]
                  << '\n';
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

std::vector<qmomkit::MultiIndex> select_terminal_moment_indices(const std::vector<qmomkit::MultiIndex>& required_indices)
{
    std::vector<qmomkit::MultiIndex> selected_indices = required_indices;
    // Stable sorting by total degree gives readable console output while preserving
    // the CQMOM hierarchy order as the tie-break within each degree. File-output
    // ordering is unaffected.
    std::stable_sort(selected_indices.begin(), selected_indices.end(), [](const qmomkit::MultiIndex& lhs, const qmomkit::MultiIndex& rhs)
                     { return moment_total_degree(lhs) < moment_total_degree(rhs); });
    if (selected_indices.size() > max_console_moments)
    {
        selected_indices.resize(max_console_moments);
    }
    return selected_indices;
}

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

void print_result_summary(const qmomkit::SelectedRawMomentResult& initial_moments, const qmomkit::SelectedRawMomentResult& final_moments,
                          const qmomkit::CqmomQuadrature& final_quadrature, const std::filesystem::path& output_directory,
                          const std::size_t sample_count, const unsigned long long initialization_seed, const double time_step,
                          const double final_time, const bool wrote_moments, const bool wrote_quadrature)
{
    std::cout << std::setprecision(12);
    std::cout << "5D CQMOM non-spherical particle example\n";
    std::cout << "phase-space ordering: (d, chi, u, v, w)\n";
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
            std::cout << (output_directory / "initial_quadrature.txt").string() << '\n';
            std::cout << (output_directory / "quadrature_history.txt").string() << '\n';
        }
    }
    else
    {
        std::cout << "\nOutput disabled by configuration.\n";
    }
}

void print_monte_carlo_result_summary(const qmomkit::RawMomentTensor& initial_moments, const qmomkit::RawMomentTensor& final_moments,
                                      const std::vector<qmomkit::MultiIndex>& terminal_moment_indices, const std::filesystem::path& output_directory,
                                      const std::size_t final_particle_count, const std::size_t sample_count,
                                      const unsigned long long initialization_seed, const double time_step, const double final_time,
                                      const bool wrote_moments)
{
    std::cout << std::setprecision(12);
    std::cout << "5D direct Monte Carlo non-spherical particle example\n";
    std::cout << "phase-space ordering: (d, chi, u, v, w)\n";
    print_run_metadata(sample_count, initialization_seed, time_step, final_time);

    std::cout << "\nInitial moments:\n";
    print_moment_tensor(initial_moments, terminal_moment_indices);

    std::cout << "\nFinal moments:\n";
    print_moment_tensor(final_moments, terminal_moment_indices);

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

} // namespace cqmom_non_spherical_particles
