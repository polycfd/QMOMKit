#include "qmomkit/qmomkit.hpp"

#include <iostream>

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
    const qmomkit::QmomCase qmom_case{qmomkit::MomentSet(initial_moments), node_count, options.qmom_options};
    qmomkit::QmomState state = qmom_case.initial_state();

    // Define the callback for the phase-space advection velocity
    const auto velocity = [](double, const std::vector<double>& coordinates) { return std::vector<double>{1.0 - coordinates[0]}; };

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
    output_options.write_moments = true;
    output_options.write_quadrature = true;
    output_options.interval = 0.1;
    const std::size_t precision = 6;

    std::filesystem::create_directories(output_options.directory);

    qmomkit::MomentWriter moment_history_writer =
        qmomkit::MomentWriter::for_qmom(output_options.directory / "moment_history.txt", state, output_options, precision);

    qmomkit::QuadratureWriter quadrature_history_writer =
        qmomkit::QuadratureWriter::for_qmom(output_options.directory / "quadrature_history.txt", qmom_case.node_count, precision);

    // Write out the initial state
    const auto& quadrature = state.quadrature();
    moment_history_writer.write(initial_time, state, quadrature);
    quadrature_history_writer.write(initial_time, quadrature);

    double next_output_time = initial_time + output_options.interval.value();
    double time = initial_time;

    // Time loop
    while (time < final_time - 1.0e-12)
    {
        const double step = std::min(time_step, final_time - time);

        // Advance the moments in phase space
        const qmomkit::QmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_qmom_phase_space(state, qmom_case.node_count, time, step, velocity, source, options);

        if (!result.ok())
        {
            std::cerr << "QMOM phase-space evolution failed: " << result.message << '\n';
            return 1;
        }

        state = result.state;
        time = result.final_time;

        // Write out the results
        if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
        {
            const auto& quadrature = state.quadrature();
            moment_history_writer.write(time, state, quadrature);
            quadrature_history_writer.write(time, quadrature);
            next_output_time += output_options.interval.value();
        }
    }

    // Screen output of the final moments and quadrature
    std::cout << '\n';
    std::cout << "Final moments:\n";
    for (std::size_t order = 0; order < state.moments().size(); ++order)
    {
        std::cout << "M" << order << " = " << state.moments()[order] << '\n';
    }
    std::cout << '\n';
    std::cout << "Final quadrature:\n";
    std::cout << "node  weight        abscissa\n";
    const auto& final_quadrature = state.quadrature();
    for (std::size_t node = 0; node < final_quadrature.weights.size(); ++node)
    {
        std::cout << node << "     " << final_quadrature.weights[node] << "   " << final_quadrature.abscissas[node] << '\n';
    }
    if (output_options.write_moments || output_options.write_quadrature)
    {
        std::cout << "\nOutputs written to:\n";
        if (output_options.write_moments)
        {
            std::cout << (output_options.directory / "moment_history.txt").string() << '\n';
        }
        if (output_options.write_quadrature)
        {
            std::cout << (output_options.directory / "quadrature_history.txt").string() << '\n';
        }
    }

    return 0;
}
