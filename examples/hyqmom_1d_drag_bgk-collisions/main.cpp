#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <utility>
#include <vector>

namespace
{

// The number density, the mean and standard deviation of the velocity, m_{2n+1}
// of the quadrature, and the characteristic speeds of a state.
bool print_state(const double time, const qmomkit::QmomState& state)
{
    const qmomkit::HyqmomCharacteristicSpeeds characteristic = qmomkit::compute_hyqmom_characteristic_speeds(state);
    if (!characteristic.ok())
    {
        std::cerr << "Characteristic speeds failed: " << characteristic.message << '\n';
        return false;
    }
    const std::vector<double>& moments = state.moment_values();
    const double mean = moments[1] / moments[0];
    const double standard_deviation = std::sqrt(moments[2] / moments[0] - mean * mean);
    // HyQMOM with node count n reconstructs n + 1 nodes from m0..m2n; their
    // m_{2n+1} is the moment that HyQMOM provides beyond those of the state.
    const std::size_t next_order = 2 * state.maximum_node_count() - 1;

    std::cout << "t = " << time << " s, " << state.active_node_count() << " nodes\n";
    std::cout << "  number density m0 (1/m^3):     " << moments[0] << '\n';
    std::cout << "  mean velocity (m/s):           " << mean << '\n';
    std::cout << "  standard deviation (m/s):      " << standard_deviation << '\n';
    std::cout << "  moment m" << next_order << " (m^" << next_order - 3 << "/s^" << next_order
              << "):           " << qmomkit::evaluate_generalized_moment(state.quadrature(), static_cast<double>(next_order)) << '\n';
    std::cout << "  characteristic speeds (m/s):  ";
    for (const double speed : characteristic.speeds)
    {
        std::cout << ' ' << speed;
    }
    std::cout << '\n';
    return true;
}

} // namespace

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 2.0;
    const double time_step = 0.01;

    // The case, the evolution options and the output options come from the same
    // file. A faulty file is reported as a message naming the key and line.
    auto configured_case = qmomkit::config::load_qmom(configuration_file);
    if (!configured_case)
    {
        std::cerr << configured_case.error() << '\n';
        return 1;
    }
    auto configured_options = qmomkit::config::load_qmom_phase_space_evolution_options(configuration_file);
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
    const qmomkit::QmomPhaseSpaceEvolutionOptions options = std::move(configured_options).value();
    const qmomkit::OutputOptions output = std::move(configured_output).value();

    // The physics is supplied as two callbacks defined in physics.cpp: the drag,
    // which moves the quadrature nodes, and the BGK collisions, which relax the
    // moments toward those of a normal distribution.
    const qmomkit::PhaseSpaceVelocity velocity = hyqmom_1d_drag_bgk_collisions::drag_velocity();
    const qmomkit::QmomMomentSource source = hyqmom_1d_drag_bgk_collisions::bgk_collision_source();

    try
    {
        // initial_state() throws when the configured moments admit no quadrature,
        // and the writers throw when a file cannot be written.
        qmomkit::QmomState state = qmom_case.initial_state();
        if (!print_state(initial_time, state))
        {
            return 1;
        }

        std::optional<qmomkit::MomentWriter> moment_writer;
        std::optional<qmomkit::QuadratureWriter> quadrature_writer;
        if (output.write_moments || output.write_quadrature)
        {
            std::filesystem::create_directories(output.directory);
        }
        if (output.write_moments)
        {
            moment_writer.emplace(qmomkit::MomentWriter::for_qmom(output.directory / "moments.txt", state, output));
        }
        if (output.write_quadrature)
        {
            quadrature_writer.emplace(qmomkit::QuadratureWriter::for_qmom(output.directory / "quadrature.txt", state.maximum_node_count()));
        }
        const auto write_output = [&](const double time)
        {
            if (moment_writer)
            {
                moment_writer->write(time, state, state.quadrature());
            }
            if (quadrature_writer)
            {
                quadrature_writer->write(time, state.quadrature());
            }
        };
        write_output(initial_time);

        const double output_interval = output.interval.value_or(final_time - initial_time);
        double next_output_time = initial_time + output_interval;
        double time = initial_time;
        while (time < final_time - 1.0e-12)
        {
            // One step of the configured splitting scheme: the drag on the nodes
            // and the BGK collisions on the moments.
            qmomkit::QmomPhaseSpaceEvolutionResult result =
                qmomkit::advance_qmom_phase_space(state, time, std::min(time_step, final_time - time), velocity, source, options);
            if (!result.ok())
            {
                std::cerr << "QMOM phase-space evolution failed: " << result.message << '\n';
                return 1;
            }
            state = std::move(result.state);
            time = result.final_time;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                write_output(time);
                next_output_time += output_interval;
            }
        }

        if (!print_state(time, state))
        {
            return 1;
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
