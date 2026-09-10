#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 20.0;
    const double time_step = 0.2;

    auto configured = qmomkit::config::load_qmom(configuration_file);
    if (!configured)
    {
        std::cerr << configured.error() << '\n';
        return 1;
    }

    auto node_advection_configured = qmomkit::config::load_node_advection_options(configuration_file);
    if (!node_advection_configured)
    {
        std::cerr << node_advection_configured.error() << '\n';
        return 1;
    }

    auto output_configured = qmomkit::config::load_output_options(configuration_file);
    if (!output_configured)
    {
        std::cerr << output_configured.error() << '\n';
        return 1;
    }

    const qmomkit::QmomCase qmom_case = std::move(configured).value();
    const qmomkit::NodeAdvectionOptions node_advection_options = std::move(node_advection_configured).value();
    const qmomkit::OutputOptions output_options = std::move(output_configured).value();

    const qmomkit::QmomResult initial_reconstruction = qmomkit::compute_qmom(qmom_case);
    if (!initial_reconstruction.ok())
    {
        std::cerr << "Initial QMOM reconstruction failed: " << initial_reconstruction.message << '\n';
        return 1;
    }
    if (initial_reconstruction.weights.size() != qmom_case.node_count)
    {
        std::cerr << "Initial QMOM reconstruction did not retain the configured fixed node count.\n";
        return 1;
    }

    const qmomkit::Quadrature1D initial_quadrature{initial_reconstruction.abscissas, initial_reconstruction.weights};
    qmomkit::Quadrature1D quadrature = initial_quadrature;
    std::vector<double> moments = qmomkit::reconstruct_raw_moments(quadrature, qmom_case.moments.values().size());
    const std::vector<double> initial_moments = moments;
    if (!std::all_of(moments.begin(), moments.end(), [](const double value) { return std::isfinite(value); }))
    {
        std::cerr << "Initial QMOM moments are non-finite.\n";
        return 1;
    }

    const qmomkit::PhaseSpaceVelocity velocity = qmom_1d_aerosol_evaporation::evaporation_velocity();

    std::optional<qmomkit::MomentWriter> moment_writer;
    std::optional<qmomkit::QuadratureWriter> quadrature_writer;
    try
    {
        if (output_options.write_moments || output_options.write_quadrature)
        {
            std::filesystem::create_directories(output_options.directory);
        }
        if (output_options.write_moments)
        {
            moment_writer.emplace(qmomkit::MomentWriter::for_qmom(output_options.directory / "moments.txt", moments.size(), output_options));
            moment_writer->write(initial_time, qmomkit::MomentSet(moments), quadrature);
        }
        if (output_options.write_quadrature)
        {
            quadrature_writer.emplace(qmomkit::QuadratureWriter::for_qmom(output_options.directory / "quadrature.txt", qmom_case.node_count));
            quadrature_writer->write(initial_time, quadrature);
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    double time = initial_time;
    const double output_interval = output_options.interval.value_or(final_time - initial_time);
    double next_output_time = initial_time + output_interval;
    while (time < final_time - 1.0e-12)
    {
        const double step = std::min(time_step, final_time - time);
        const qmomkit::QmomNodeAdvectionResult advection = qmomkit::advance_nodes(quadrature, time, step, velocity, node_advection_options);
        if (!advection.ok())
        {
            std::cerr << "QMOM node advection failed: " << advection.message << '\n';
            return 1;
        }

        quadrature = advection.quadrature;
        qmom_1d_aerosol_evaporation::enforce_residual_diameter(quadrature);
        if (quadrature.weights.size() != qmom_case.node_count)
        {
            std::cerr << "Fixed QMOM node count changed during evaporation.\n";
            return 1;
        }

        moments = qmomkit::reconstruct_raw_moments(quadrature, qmom_case.moments.values().size());
        if (!std::all_of(moments.begin(), moments.end(), [](const double value) { return std::isfinite(value); }))
        {
            std::cerr << "Evaporated QMOM moments are non-finite.\n";
            return 1;
        }
        time += step;
        if (time >= final_time - 1.0e-12)
        {
            time = final_time;
        }

        if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
        {
            try
            {
                if (moment_writer.has_value())
                {
                    moment_writer->write(time, qmomkit::MomentSet(moments), quadrature);
                }
                if (quadrature_writer.has_value())
                {
                    quadrature_writer->write(time, quadrature);
                }
            }
            catch (const std::exception& error)
            {
                std::cerr << error.what() << '\n';
                return 1;
            }
            while (next_output_time <= time + 1.0e-12)
            {
                next_output_time += output_interval;
            }
        }
    }

    std::cout << std::setprecision(12);
    std::cout << "1D QMOM aerosol evaporation\n";
    std::cout << "configured QMOM nodes: " << qmom_case.node_count << '\n';
    std::cout << "initial: time=" << initial_time << " M0=" << initial_moments[0] << " mean_diameter=" << initial_moments[1] / initial_moments[0]
              << " um\n";
    std::cout << "final:   time=" << time << " M0=" << moments[0] << " mean_diameter=" << moments[1] / moments[0] << " um\n";

    std::cout << "\nInitial quadrature:\nnode  weight        diameter (um)\n";
    for (std::size_t node = 0; node < initial_quadrature.weights.size(); ++node)
    {
        std::cout << node << "     " << initial_quadrature.weights[node] << "   " << initial_quadrature.abscissas[node] << '\n';
    }

    std::cout << "\nFinal quadrature:\nnode  weight        diameter (um)\n";
    for (std::size_t node = 0; node < quadrature.weights.size(); ++node)
    {
        std::cout << node << "     " << quadrature.weights[node] << "   " << quadrature.abscissas[node] << '\n';
    }

    if (output_options.write_moments || output_options.write_quadrature)
    {
        std::cout << "\nOutputs written to:\n";
        if (output_options.write_moments)
        {
            std::cout << (output_options.directory / "moments.txt").string() << '\n';
        }
        if (output_options.write_quadrature)
        {
            std::cout << (output_options.directory / "quadrature.txt").string() << '\n';
        }
    }

    return 0;
}
