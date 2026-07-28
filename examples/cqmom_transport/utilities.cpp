#include "utilities.hpp"

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <system_error>
#include <vector>

namespace cqmom_transport
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

void print_quadrature(const qmomkit::CqmomQuadrature& quadrature)
{
    std::cout << "node  path       weight";
    for (std::size_t dimension = 0; dimension < quadrature.number_of_dimensions(); ++dimension)
    {
        std::cout << "          x" << dimension;
    }
    std::cout << '\n';

    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        std::cout << node << "     ";
        print_path(quadrature.terminal_paths()[node]);
        std::cout << "   " << quadrature.weights()[node];
        for (std::size_t dimension = 0; dimension < quadrature.number_of_dimensions(); ++dimension)
        {
            std::cout << "   " << quadrature.abscissas()[dimension][node];
        }
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
                          const qmomkit::CqmomQuadrature& final_quadrature, const std::filesystem::path& output_directory, const bool wrote_moments,
                          const bool wrote_quadrature)
{
    std::cout << std::setprecision(12);
    std::cout << "CQMOM continuous + discontinuous evolution example\n";
    std::cout << "phase-space dimensions: " << final_quadrature.number_of_dimensions() << '\n';
    std::cout << "breakage coordinate: dimension 0\n";

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

} // namespace cqmom_transport
