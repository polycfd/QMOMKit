#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <utility>
#include <vector>

int main()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {0, 1});
    const qmomkit::RawMomentTensor initial_moments{{4, 4}, {1.0, 0.0, 0.25, 0.0, 2.0, 0.0, 0.5, 0.0, 4.25, 0.0, 1.0625, 0.0, 9.5, 0.0, 2.375, 0.0}};

    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.01;

    // Create the CASE and initialize the STATE.
    //
    // CQMOM reconstructs from a fixed subset of the moment tensor rather than
    // from all of it: for this hierarchy, 10 of the 16 entries. The layout names
    // exactly which moments those are, and extract_cqmom_moments picks them out
    // of the tensor above. A caller whose moments already arrive as that subset
    // can build the case from them directly and skip the tensor entirely.
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::CqmomLeanCase case_definition{qmomkit::extract_cqmom_moments(initial_moments, layout),
                                           {},                                    // conditional-moment options
                                           {},                                    // QMOM options
                                           false,                                 // adaptive branch counts
                                           qmomkit::CqmomDiagnosticLevel::none,   // diagnostics
                                           qmomkit::VandermondeSolver::automatic, // solver policy
                                           {}};                                   // coordinate support
    const qmomkit::CqmomPreparedCase cqmom_prepared_case{std::move(case_definition)};
    qmomkit::CqmomState state = cqmom_prepared_case.initial_state();

    // Define the callback for the phase-space advection velocity. The callback
    // receives every terminal node at once: `coordinates` holds
    // `nodes * dimensions` values packed node-major, so node `i` occupies
    // `[i * dimensions, (i + 1) * dimensions)`. Fill `velocities` with dxi/dt
    // in the same layout, leaving no entry unwritten.
    //
    // This model is two-dimensional, so each node owns exactly two components
    // and both are written explicitly. A model of another dimensionality would
    // loop over its own components instead.
    const qmomkit::PhaseSpaceVelocity velocity =
        [](double, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
    {
        for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
        {
            const std::size_t first = node * dimensions;
            velocities[first] = 1.0 - coordinates[first];
            velocities[first + 1] = -coordinates[first + 1];
        }
    };

    // Define the callback for the moment source
    const qmomkit::CqmomMomentSource source = [](double, const qmomkit::CqmomSourceContext& context)
    {
        constexpr double fragmentation_rate = 0.5;
        std::vector<double> derivative(context.moment_count(), 0.0);
        for (std::size_t index = 0; index < context.moment_count(); ++index)
        {
            const double breakage_coordinate_order = static_cast<double>(context.index(index)[0]);
            const double daughter_factor = std::pow(2.0, 1.0 - breakage_coordinate_order);
            derivative[index] = fragmentation_rate * (daughter_factor - 1.0) * context.value(index);
        }
        return qmomkit::CqmomMomentRate{std::move(derivative)};
    };

    // Define the callback for the moment boundary contribution
    const qmomkit::CqmomMomentBoundaryContribution boundary = [](double, const qmomkit::CqmomSourceContext& context)
    {
        constexpr double injection_rate = 0.1;
        constexpr double injection_size = 0.1;
        std::vector<double> rate(context.moment_count(), 0.0);
        for (std::size_t index = 0; index < context.moment_count(); ++index)
        {
            rate[index] = injection_rate * std::pow(injection_size, static_cast<double>(context.index(index)[0]));
        }
        return qmomkit::CqmomMomentRate{std::move(rate)};
    };

    // Configure the options for the run
    qmomkit::CqmomPhaseSpaceEvolutionOptions evolution_options;
    evolution_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    // Configure the output
    qmomkit::OutputOptions output_options;
    output_options.directory = "output";
    output_options.interval = 0.1;
    const std::size_t precision = 6;

    std::filesystem::create_directories(output_options.directory);

    qmomkit::MomentWriter moment_history_writer =
        qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history.txt", state, output_options, precision);
    qmomkit::QuadratureWriter quadrature_history_writer = qmomkit::QuadratureWriter::for_cqmom(output_options.directory / "quadrature_history.txt",
                                                                                               state.quadrature().number_of_dimensions(), precision);

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
        qmomkit::CqmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_cqmom_phase_space(state, time, step, velocity, source, boundary, evolution_options);

        if (!result.ok())
        {
            std::cerr << "CQMOM phase-space evolution failed: " << result.message << '\n';
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

    // Confirm on the console that the run finished and say where the results are.
    std::cout << "Run complete. Results written to " << output_options.directory.string() << "/:\n"
              << "  moment_history.txt\n"
              << "  quadrature_history.txt\n";

    return 0;
}
