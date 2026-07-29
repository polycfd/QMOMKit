#include "utilities.hpp"

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <system_error>

namespace cqmom_2d_vanishing_variance
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

void print_paths(const std::vector<qmomkit::ConditioningPath>& paths)
{
    for (std::size_t path = 0; path < paths.size(); ++path)
    {
        if (path > 0)
        {
            std::cout << ',';
        }
        print_path(paths[path]);
    }
}

void print_conditioning_order(const std::vector<std::size_t>& conditioning_order)
{
    std::cout << '[';
    for (std::size_t dimension = 0; dimension < conditioning_order.size(); ++dimension)
    {
        if (dimension > 0)
        {
            std::cout << ',';
        }
        std::cout << conditioning_order[dimension];
    }
    std::cout << ']';
}

void print_selected_moments(const char* label, const qmomkit::SelectedRawMomentResult& moments)
{
    std::cout << label << '\n';
    if (!moments.ok())
    {
        std::cout << "  unavailable: " << moments.message << '\n';
        return;
    }

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

void print_quadrature(const qmomkit::CqmomQuadrature& quadrature)
{
    std::cout << "node  path       weight          d               u\n";
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        std::cout << node << "     ";
        print_path(quadrature.terminal_paths()[node]);
        std::cout << "   " << quadrature.weights()[node] << "   " << quadrature.abscissas()[0][node] << "   " << quadrature.abscissas()[1][node]
                  << '\n';
    }
}

double selected_moment(const qmomkit::SelectedRawMomentResult& moments, const qmomkit::MultiIndex& index)
{
    for (std::size_t value = 0; value < moments.values.size(); ++value)
    {
        if (moments.indices[value] == index)
        {
            return moments.values[value];
        }
    }
    throw std::runtime_error("Requested diagnostic moment is unavailable.");
}

void print_variance_history(const std::vector<std::pair<double, qmomkit::CqmomQuadrature>>& quadrature_history)
{
    std::cout << "\nVariance history:\n";
    std::cout << "time  size_variance  velocity_variance  active_paths\n";
    const qmomkit::RawMomentRequest request{{{0, 0}, {1, 0}, {2, 0}, {0, 1}, {0, 2}}};
    for (const auto& [time, quadrature] : quadrature_history)
    {
        const qmomkit::SelectedRawMomentResult moments = qmomkit::reconstruct_raw_moments(quadrature, request);
        if (!moments.ok())
        {
            std::cout << time << "  unavailable: " << moments.message << '\n';
            continue;
        }

        const double mass = selected_moment(moments, {0, 0});
        const double size_mean = selected_moment(moments, {1, 0}) / mass;
        const double velocity_mean = selected_moment(moments, {0, 1}) / mass;
        const double size_variance = selected_moment(moments, {2, 0}) / mass - size_mean * size_mean;
        const double velocity_variance = selected_moment(moments, {0, 2}) / mass - velocity_mean * velocity_mean;
        std::cout << time << "  " << size_variance << "  " << velocity_variance << "  ";
        print_paths(quadrature.terminal_paths());
        std::cout << '\n';
    }
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

void print_result_summary(const qmomkit::SelectedRawMomentResult& initial_moments, const qmomkit::SelectedRawMomentResult& final_moments,
                          const std::vector<std::pair<double, qmomkit::CqmomQuadrature>>& quadrature_history,
                          const qmomkit::CqmomQuadrature& final_quadrature, const std::filesystem::path& output_directory,
                          const std::vector<std::size_t>& conditioning_order, const double time_step, const double final_time,
                          const bool wrote_moments, const bool wrote_quadrature)
{
    std::cout << std::setprecision(12);
    std::cout << "2D CQMOM vanishing velocity-variance example\n";
    std::cout << "phase-space ordering: (d,u)\n";
    std::cout << "conditioning_order = ";
    print_conditioning_order(conditioning_order);
    std::cout << '\n';
    std::cout << "time_step = " << time_step << '\n';
    std::cout << "final_time = " << final_time << '\n';

    std::cout << '\n';
    print_selected_moments("Initial moments:", initial_moments);

    std::cout << '\n';
    print_selected_moments("Final moments:", final_moments);

    print_variance_history(quadrature_history);

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

} // namespace cqmom_2d_vanishing_variance
