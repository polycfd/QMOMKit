#include "host_data.hpp"
#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.01;
    // The locations the host owns, for example the mesh nodes held by one
    // processor. The host advances all of them together with one batch call.
    const std::size_t location_count = 8;

    // The case, the evolution options and the output options come from the same
    // file. A faulty file is reported as a message naming the key and line.
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
    const qmomkit::CqmomPreparedCase cqmom_case = std::move(configured_case).value();
    const qmomkit::CqmomPhaseSpaceEvolutionOptions options = std::move(configured_options).value();
    const qmomkit::OutputOptions output = std::move(configured_output).value();

    // The host packages its data for the locations it owns. Entry i of every array
    // in host_data belongs to location i. The batch calls pass each state's
    // position in the state vector to the callbacks as batch_index, so the
    // callbacks read entry batch_index of host_data. They hold a reference to
    // host_data and read its current contents on every invocation.
    cqmom_2d_batch_host_data::LocalHostData host_data = cqmom_2d_batch_host_data::make_local_host_data(location_count, initial_time);
    const qmomkit::BatchPhaseSpaceVelocity velocity = cqmom_2d_batch_host_data::advection_velocity(host_data);
    const qmomkit::BatchCqmomMomentSource source = cqmom_2d_batch_host_data::binary_fragmentation_source(host_data);

    try
    {
        // All locations start from the configured case and share its prepared
        // branch structure. initial_state() throws when the configured moments
        // admit no quadrature, and the writers throw when a file cannot be written.
        std::vector<qmomkit::CqmomState> states(location_count, cqmom_case.initial_state());
        const double initial_number_density = states.front().moment({0, 0}).value();
        const double initial_mean_velocity = states.front().moment({0, 1}).value() / initial_number_density;

        // One moment and one quadrature history per location.
        std::vector<qmomkit::MomentWriter> moment_writers;
        std::vector<qmomkit::QuadratureWriter> quadrature_writers;
        if (output.write_moments || output.write_quadrature)
        {
            std::filesystem::create_directories(output.directory);
        }
        for (std::size_t location = 0; location < location_count; ++location)
        {
            if (output.write_moments)
            {
                moment_writers.push_back(qmomkit::MomentWriter::for_cqmom(
                    cqmom_2d_batch_host_data::location_output_file(output.directory, location, "moment_history.txt"), states[location], output));
            }
            if (output.write_quadrature)
            {
                quadrature_writers.push_back(qmomkit::QuadratureWriter::for_cqmom(
                    cqmom_2d_batch_host_data::location_output_file(output.directory, location, "quadrature_history.txt"),
                    states[location].hierarchy().number_of_dimensions()));
            }
        }
        const auto write_output = [&](const double time)
        {
            for (std::size_t location = 0; location < location_count; ++location)
            {
                const qmomkit::CqmomState& state = states[location];
                if (output.write_moments)
                {
                    moment_writers[location].write(time, state, state.quadrature());
                }
                if (output.write_quadrature)
                {
                    quadrature_writers[location].write(time, state.quadrature());
                }
            }
        };
        write_output(initial_time);

        const double output_interval = output.interval.value_or(final_time - initial_time);
        double next_output_time = initial_time + output_interval;
        double time = initial_time;
        while (time < final_time - 1.0e-12)
        {
            const double step = std::min(time_step, final_time - time);

            // The host advances its flow field and publishes the new level at all
            // of its locations, keeping the level it replaces.
            if (time > initial_time)
            {
                cqmom_2d_batch_host_data::publish_carrier_flow_level(host_data, time);
            }

            // One call advances every location by one step of the configured
            // splitting scheme. Each entry carries its own status, so every
            // location is checked. What to do with a failed location is the
            // host's decision; this example stops.
            std::vector<qmomkit::CqmomPhaseSpaceEvolutionResult> results =
                qmomkit::advance_cqmom_phase_space_batch_with_index(states, time, step, velocity, source, options);
            for (std::size_t location = 0; location < location_count; ++location)
            {
                if (!results[location].ok())
                {
                    std::cerr << "CQMOM phase-space evolution failed at location " << location << ": " << results[location].message << '\n';
                    return 1;
                }
                states[location] = std::move(results[location].state);
            }
            time += step;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                write_output(time);
                next_output_time += output_interval;
            }
        }

        std::cout << "Advanced " << location_count << " locations to t = " << time << " with one batch call per step.\n\n";
        cqmom_2d_batch_host_data::print_location_summary(host_data, states, initial_time, initial_number_density, initial_mean_velocity, time);
        if (output.write_moments || output.write_quadrature)
        {
            std::cout << "\nPer-location histories written to " << output.directory.string() << "/location_<i>_*.txt\n";
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
