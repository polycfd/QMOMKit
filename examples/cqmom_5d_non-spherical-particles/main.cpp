#include "breakage_model.hpp"
#include "fluid_field.hpp"
#include "particle_model.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <utility>

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <path-to-configuration.qmom>\n";
        return 1;
    }
    const std::filesystem::path configuration_file = argv[1];
    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.005;

    // configuration.qmom defines a five-dimensional phase space (d, chi, u, v, w)
    // with three nodes per dimension, so the quadrature holds at most 3^5 = 243
    // terminal nodes. The callbacks built further below receive all of them in a
    // single call, packed node-major. A faulty file is reported as a message
    // naming the key and line.
    auto configured_case = qmomkit::config::load_cqmom_prepared(configuration_file);
    if (!configured_case)
    {
        std::cerr << configured_case.error() << '\n';
        return 1;
    }
    auto configured_options = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
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
    const qmomkit::CqmomPhaseSpaceEvolutionOptions options = std::move(configured_options).value();
    const qmomkit::OutputOptions output = std::move(configured_output).value();

    // CQMOM consumes only the 606 moments its hierarchy requires, so the sampler is
    // asked for exactly those indices. Both runs start from the same sample.
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(configured_case.value().hierarchy());
    const qmomkit::RawMomentRequest required_moments{layout.required_indices.to_multi_indices()};
    const std::optional<qmomkit::monte_carlo::SelectedSampledInitialization> sampled =
        cqmom_5d_non_spherical_particles::sample_requested_initial_moments(configuration_file, required_moments);
    if (!sampled)
    {
        return 1;
    }
    auto sampled_case = qmomkit::config::with_initial_moments(
        configured_case.value(),
        qmomkit::create_cqmom_moment_data(layout, qmomkit::make_selected_raw_moment_result(required_moments, sampled->moments.values)));
    if (!sampled_case)
    {
        std::cerr << sampled_case.error() << '\n';
        return 1;
    }
    const qmomkit::CqmomPreparedCase cqmom_case = std::move(sampled_case).value();

    // carrier_samples is the host-side carrier flow, owned here by the driver.
    // The velocity callback keeps a reference to it, which is how a host passes
    // local, time-varying field data to a callback without a global. It must
    // outlive that callback. Only the newest level is known at this point, so
    // the first step runs without a reconstruction.
    cqmom_5d_non_spherical_particles::CarrierFlowSamples carrier_samples;
    carrier_samples.current = cqmom_5d_non_spherical_particles::carrier_flow_velocity(initial_time);
    carrier_samples.sample_time = initial_time;
    const qmomkit::PhaseSpaceVelocity velocity = cqmom_5d_non_spherical_particles::particle_phase_space_velocity(
        cqmom_5d_non_spherical_particles::particle_model_parameters(), carrier_samples);
    const qmomkit::CqmomMomentSource source =
        cqmom_5d_non_spherical_particles::binary_breakage_source(cqmom_5d_non_spherical_particles::binary_fragmentation_parameters());
    const qmomkit::RawMomentRequest printed_moments = cqmom_5d_non_spherical_particles::printed_moments();

    try
    {
        // initial_state() throws when the sampled moments admit no quadrature, and
        // the writers throw when a file cannot be written.
        qmomkit::CqmomState state = cqmom_case.initial_state();
        const qmomkit::SelectedRawMomentResult initial_moments = qmomkit::reconstruct_raw_moments(state.quadrature(), printed_moments);
        if (!initial_moments.ok())
        {
            std::cerr << "Initial moment reconstruction failed: " << initial_moments.message << '\n';
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
            moment_writer.emplace(qmomkit::MomentWriter::for_cqmom(output.directory / "moment_history.txt", state, output));
        }
        if (output.write_quadrature)
        {
            quadrature_writer.emplace(
                qmomkit::QuadratureWriter::for_cqmom(output.directory / "quadrature_history.txt", state.hierarchy().number_of_dimensions()));
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
            // The host advances its own field and publishes the new level, keeping
            // the level it replaces. Those two levels are all a coupled host has:
            // it cannot evaluate at QMOMKit's internal stage times, and it does
            // not know the field at time + step. The callback extends the line
            // through them, which is second-order accurate in the step and so not
            // the limiting error of the second-order Strang-split step. Freezing a
            // single value here instead would make the coupling first-order in u
            // and v, the components in which this flow varies.
            if (time > initial_time)
            {
                carrier_samples.previous = carrier_samples.current;
                carrier_samples.current = cqmom_5d_non_spherical_particles::carrier_flow_velocity(time);
                carrier_samples.sample_interval = time - carrier_samples.sample_time;
                carrier_samples.sample_time = time;
            }

            qmomkit::CqmomPhaseSpaceEvolutionResult result =
                qmomkit::advance_cqmom_phase_space(state, time, std::min(time_step, final_time - time), velocity, source, options);
            if (!result.ok())
            {
                std::cerr << "CQMOM phase-space evolution failed: " << result.message << '\n';
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

        const qmomkit::SelectedRawMomentResult final_moments = qmomkit::reconstruct_raw_moments(state.quadrature(), printed_moments);
        if (!final_moments.ok())
        {
            std::cerr << "Final moment reconstruction failed: " << final_moments.message << '\n';
            return 1;
        }
        cqmom_5d_non_spherical_particles::print_moments(std::cout, "QMOMKit run", initial_moments.values, final_moments.values);
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
