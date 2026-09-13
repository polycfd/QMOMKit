#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <utility>
#include <vector>

int main()
{
    const std::size_t node_count = 3;
    const std::vector<double> initial_moments{1.0, 2.0, 5.025, 14.195, 42.11255, 127.48675};

    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.01;

    // Configure the options for the run (default options are used here)
    qmomkit::QmomPhaseSpaceEvolutionOptions options;

    // Create the CASE and initialize the STATE
    const qmomkit::QmomCase qmom_case{qmomkit::MomentSet(initial_moments), node_count};
    qmomkit::QmomState state = qmom_case.initial_state();

    // Define the callback for the phase-space advection velocity. The callback
    // receives every quadrature node at once: `coordinates` holds
    // `nodes * dimensions` values packed node-major, and `velocities` must be
    // filled with dxi/dt in the same layout.
    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>& coordinates, std::vector<double>& velocities, std::size_t)
    {
        for (std::size_t node = 0; node < coordinates.size(); ++node)
        {
            velocities[node] = 1.0 - coordinates[node];
        }
    };

    // Define the callback for the moment source
    const auto source = [](double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> derivative(context.size(), 0.0);
        for (std::size_t order = 0; order < context.size(); ++order)
        {
            const double daughter_factor = std::pow(2.0, 1.0 - static_cast<double>(order));
            derivative[order] = 0.5 * (daughter_factor - 1.0) * context.value(order);
        }
        return derivative;
    };

    // Configure the output
    qmomkit::OutputOptions output_options;
    output_options.directory = "output";
    output_options.interval = 0.1;
    const std::size_t precision = 6;

    std::filesystem::create_directories(output_options.directory);

    qmomkit::MomentWriter moment_history_writer =
        qmomkit::MomentWriter::for_qmom(output_options.directory / "moment_history.txt", state, output_options, precision);

    qmomkit::QuadratureWriter quadrature_history_writer =
        qmomkit::QuadratureWriter::for_qmom(output_options.directory / "quadrature_history.txt", state.maximum_node_count(), precision);

    // Write out the initial state
    moment_history_writer.write(initial_time, state, state.quadrature());
    quadrature_history_writer.write(initial_time, state.quadrature());

    double time = initial_time;
    double next_output_time = initial_time + output_options.interval.value();

    // Time loop
    while (time < final_time - 1.0e-12)
    {
        const double step = std::min(time_step, final_time - time);

        // Advance the moments in phase space
        qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_phase_space(state, time, step, velocity, source, options);

        if (!result.ok())
        {
            std::cerr << "QMOM phase-space evolution failed: " << result.message << '\n';
            return 1;
        }

        state = std::move(result.state);
        time = result.final_time;

        // Write out the results
        if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
        {
            moment_history_writer.write(time, state, state.quadrature());
            quadrature_history_writer.write(time, state.quadrature());
            next_output_time = time + output_options.interval.value();
        }
    }

    return 0;
}
