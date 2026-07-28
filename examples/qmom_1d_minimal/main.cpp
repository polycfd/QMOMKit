#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

int main()
{
    const std::size_t node_count = 3;
    std::vector<double> moments{1.0, 2.0, 5.025, 14.195, 42.11255, 127.48675};

    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.01;
    double time = initial_time;

    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>& coordinates)
    {
        if (coordinates.size() != 1)
        {
            throw std::invalid_argument("The minimal QMOM example expects one phase-space coordinate.");
        }

        const double diameter = coordinates[0];
        return std::vector<double>{1.0 - diameter};
    };

    const qmomkit::MomentSource source = [](double, const std::vector<double>& current_moments)
    {
        constexpr double fragmentation_rate = 0.5;

        std::vector<double> derivative(current_moments.size(), 0.0);
        for (std::size_t order = 0; order < current_moments.size(); ++order)
        {
            const double daughter_factor = std::pow(2.0, 1.0 - static_cast<double>(order));
            derivative[order] = fragmentation_rate * (daughter_factor - 1.0) * current_moments[order];
        }
        return derivative;
    };

    qmomkit::QmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::strang;

    std::cout << std::setprecision(12);
    std::cout << "Minimal 1D QMOM transport example\n\n";
    std::cout << "Initial moments:\n";
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        std::cout << "M" << order << " = " << moments[order] << '\n';
    }

    const std::filesystem::path output_directory = "output";
    std::filesystem::create_directories(output_directory);
    qmomkit::MomentWriter moment_history_writer = qmomkit::MomentWriter::for_qmom(output_directory / "moment_history.txt", moments.size());
    qmomkit::QuadratureWriter quadrature_history_writer =
        qmomkit::QuadratureWriter::for_qmom(output_directory / "quadrature_history.txt", node_count);

    const qmomkit::QmomResult initial_quadrature = qmomkit::compute_qmom(qmomkit::MomentSet(moments), node_count, options.qmom_options);
    if (!initial_quadrature.ok())
    {
        std::cerr << "Initial QMOM reconstruction failed: " << initial_quadrature.message << '\n';
        return 1;
    }
    moment_history_writer.write(initial_time, qmomkit::MomentSet(moments));
    quadrature_history_writer.write(initial_time, initial_quadrature);
    qmomkit::QmomResult final_quadrature = initial_quadrature;

    while (time < final_time - 1.0e-12)
    {
        const double step = std::min(time_step, final_time - time);
        const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(moments, node_count, time, step, velocity, source, options);
        if (!result.ok())
        {
            std::cerr << "QMOM transport failed: " << result.message << '\n';
            return 1;
        }

        moments = result.moments;
        time = result.final_time;
        final_quadrature = result.quadrature;
        moment_history_writer.write(time, qmomkit::MomentSet(moments));
        quadrature_history_writer.write(time, result.quadrature);
    }

    std::cout << '\n';
    std::cout << "Final moments:\n";
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        std::cout << "M" << order << " = " << moments[order] << '\n';
    }
    std::cout << '\n';
    std::cout << "Final quadrature:\n";
    std::cout << "node  weight        abscissa\n";
    for (std::size_t node = 0; node < final_quadrature.weights.size(); ++node)
    {
        std::cout << node << "     " << final_quadrature.weights[node] << "   " << final_quadrature.abscissas[node] << '\n';
    }
    std::cout << "\nOutputs written to:\n";
    std::cout << (output_directory / "moment_history.txt").string() << '\n';
    std::cout << (output_directory / "quadrature_history.txt").string() << '\n';

    return 0;
}
