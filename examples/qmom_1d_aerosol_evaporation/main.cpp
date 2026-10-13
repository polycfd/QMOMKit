#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>
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

    // The case, the node-advection options and the output options come from the
    // same file. A faulty file is reported as a message naming the key and line.
    auto configured_case = qmomkit::config::load_qmom(configuration_file);
    if (!configured_case)
    {
        std::cerr << configured_case.error() << '\n';
        return 1;
    }
    auto configured_options = qmomkit::config::load_node_advection_options(configuration_file);
    if (!configured_options)
    {
        std::cerr << configured_options.error() << '\n';
        return 1;
    }
    auto configured_output = qmomkit::config::load_output_options(configuration_file);
    if (!configured_output)
    {
        std::cerr << configured_output.error() << '\n';
        return 1;
    }
    const qmomkit::QmomCase qmom_case = std::move(configured_case).value();
    const qmomkit::NodeAdvectionOptions options = std::move(configured_options).value();
    const qmomkit::OutputOptions output = std::move(configured_output).value();

    // This example deliberately works one level below `QmomState`. Evaporation is
    // pure node motion, so it inverts the moments once and then advances the
    // quadrature directly, reconstructing moments from it whenever output is due.
    const qmomkit::QmomResult reconstruction = qmomkit::compute_qmom(qmom_case);
    if (!reconstruction.ok())
    {
        std::cerr << "Initial QMOM reconstruction failed: " << reconstruction.message << '\n';
        return 1;
    }
    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = reconstruction.abscissas;
    quadrature.weights = reconstruction.weights;
    const std::size_t moment_count = qmom_case.moments.values().size();
    const std::vector<double> initial_moments = qmomkit::reconstruct_raw_moments(quadrature, moment_count);

    // The only physics is a shrink rate for the droplet diameter.
    const qmomkit::PhaseSpaceVelocity velocity = qmom_1d_aerosol_evaporation::evaporation_velocity();

    try
    {
        // The writers throw when a file cannot be written.
        std::optional<qmomkit::MomentWriter> moment_writer;
        std::optional<qmomkit::QuadratureWriter> quadrature_writer;
        if (output.write_moments || output.write_quadrature)
        {
            std::filesystem::create_directories(output.directory);
        }
        if (output.write_moments)
        {
            moment_writer.emplace(qmomkit::MomentWriter::for_qmom(output.directory / "moments.txt", moment_count, output));
        }
        if (output.write_quadrature)
        {
            quadrature_writer.emplace(qmomkit::QuadratureWriter::for_qmom(output.directory / "quadrature.txt", qmom_case.node_count));
        }
        const auto write_output = [&](const double time)
        {
            if (moment_writer)
            {
                moment_writer->write(time, qmomkit::MomentSet(qmomkit::reconstruct_raw_moments(quadrature, moment_count)), quadrature);
            }
            if (quadrature_writer)
            {
                quadrature_writer->write(time, quadrature);
            }
        };
        write_output(initial_time);

        const double output_interval = output.interval.value_or(final_time - initial_time);
        double next_output_time = initial_time + output_interval;
        double time = initial_time;
        while (time < final_time - 1.0e-12)
        {
            const double step = std::min(time_step, final_time - time);
            // advance_nodes(...) moves the abscissas and leaves the weights
            // untouched, which is what makes this the right operation for
            // evaporation: droplets shrink, but none are created or destroyed.
            const qmomkit::QmomNodeAdvectionResult advection = qmomkit::advance_nodes(quadrature, time, step, velocity, options);
            if (!advection.ok())
            {
                std::cerr << "QMOM node advection failed: " << advection.message << '\n';
                return 1;
            }
            quadrature = advection.quadrature;
            // Application-side physics, not a QMOMKit step: hold shrinking droplets
            // at a residual diameter so abscissas cannot cross zero.
            qmom_1d_aerosol_evaporation::enforce_residual_diameter(quadrature);
            time += step;
            if (time >= final_time - 1.0e-12)
            {
                time = final_time;
            }

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                write_output(time);
                next_output_time += output_interval;
            }
        }

        const std::vector<double> final_moments = qmomkit::reconstruct_raw_moments(quadrature, moment_count);
        std::cout << "Number density M0:         " << initial_moments[0] << " -> " << final_moments[0] << '\n';
        std::cout << "Mean diameter M1/M0 (um):  " << initial_moments[1] / initial_moments[0] << " -> " << final_moments[1] / final_moments[0]
                  << '\n';
        std::cout << "Final quadrature at t = " << time << " s:\n";
        for (std::size_t node = 0; node < quadrature.weights.size(); ++node)
        {
            std::cout << "  weight " << quadrature.weights[node] << " at d = " << quadrature.abscissas[node] << " um\n";
        }
        if (moment_writer || quadrature_writer)
        {
            std::cout << "Output written to " << output.directory.string() << "/\n";
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
