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

namespace
{
constexpr double initial_time = 0.0;
constexpr double final_time = 1.0;
constexpr double time_step = 0.01;

// The locations the host owns, for example the mesh nodes held by one
// processor. The host advances all of them together with one batch call.
constexpr std::size_t location_count = 8;
} // namespace

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";

    try
    {
        using namespace cqmom_2d_batch_host_data;

        // The case, the evolution options and the output options come from the
        // same file. A faulty file is reported as a message naming the key and
        // line.
        auto configured = qmomkit::config::load_cqmom_prepared(configuration_file);
        if (!configured)
        {
            std::cerr << configured.error() << '\n';
            return 1;
        }
        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            std::cerr << output_configured.error() << '\n';
            return 1;
        }
        auto evolution_configured = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        if (!evolution_configured)
        {
            std::cerr << evolution_configured.error() << '\n';
            return 1;
        }
        const qmomkit::CqmomPreparedCase cqmom_prepared_case = std::move(configured).value();
        const qmomkit::OutputOptions output_options = std::move(output_configured).value();
        const qmomkit::CqmomPhaseSpaceEvolutionOptions evolution_options = std::move(evolution_configured).value();

        // The host packages its data for the locations it owns. Entry i of every
        // array in host_data belongs to location i, and states[i] is the CQMOM
        // state of that same location. All states start from the configured case
        // and share its prepared branch structure.
        LocalHostData host_data = make_local_host_data(location_count, initial_time);
        std::vector<qmomkit::CqmomState> states(location_count, cqmom_prepared_case.initial_state());
        const double initial_number_density = states.front().moment({0, 0}).value();
        const double initial_mean_velocity = states.front().moment({0, 1}).value() / initial_number_density;

        // The batch calls pass each state's position in states to the callbacks
        // as batch_index, so the callbacks read entry batch_index of host_data.
        // They hold a reference to host_data and read its current contents on
        // every invocation.
        const qmomkit::BatchPhaseSpaceVelocity velocity = advection_velocity(host_data);
        const qmomkit::BatchCqmomMomentSource source = binary_fragmentation_source(host_data);

        if (output_options.write_moments || output_options.write_quadrature)
        {
            std::filesystem::create_directories(output_options.directory);
        }
        std::vector<qmomkit::MomentWriter> moment_writers;
        std::vector<qmomkit::QuadratureWriter> quadrature_writers;
        for (std::size_t location = 0; location < location_count; ++location)
        {
            const qmomkit::CqmomState& state = states[location];
            if (output_options.write_moments)
            {
                moment_writers.push_back(qmomkit::MomentWriter::for_cqmom(
                    location_output_file(output_options.directory, location, "moment_history.txt"), state, output_options));
                moment_writers.back().write(initial_time, state, state.quadrature());
            }
            if (output_options.write_quadrature)
            {
                quadrature_writers.push_back(qmomkit::QuadratureWriter::for_cqmom(
                    location_output_file(output_options.directory, location, "quadrature_history.txt"), state.quadrature().number_of_dimensions()));
                quadrature_writers.back().write(initial_time, state.quadrature());
            }
        }

        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double next_output_time = initial_time + output_interval;
        double time = initial_time;
        while (time < final_time - 1.0e-12)
        {
            const double step = std::min(time_step, final_time - time);

            // The host advances its flow field and publishes the new level at all
            // of its locations, keeping the level it replaces.
            if (time > initial_time)
            {
                publish_carrier_flow_level(host_data, time);
            }

            // One call advances every location by one Strang-split step.
            std::vector<qmomkit::CqmomPhaseSpaceEvolutionResult> evolved =
                qmomkit::advance_cqmom_phase_space_batch_with_index(states, time, step, velocity, source, evolution_options);

            // Each entry carries its own status, so every location is checked.
            // What to do with a failed location is the host's decision; this
            // example stops.
            for (std::size_t location = 0; location < location_count; ++location)
            {
                if (!evolved[location].ok())
                {
                    std::cerr << "CQMOM phase-space evolution failed at location " << location << ": " << evolved[location].message << '\n';
                    return 1;
                }
                states[location] = std::move(evolved[location].state);
            }
            time += step;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                for (std::size_t location = 0; location < location_count; ++location)
                {
                    const qmomkit::CqmomState& state = states[location];
                    if (output_options.write_moments)
                    {
                        moment_writers[location].write(time, state, state.quadrature());
                    }
                    if (output_options.write_quadrature)
                    {
                        quadrature_writers[location].write(time, state.quadrature());
                    }
                }
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
            }
        }

        std::cout << "Advanced " << location_count << " locations to t = " << time << " with one batch call per step.\n\n";
        print_location_summary(host_data, states, initial_time, initial_number_density, initial_mean_velocity, time);
        if (output_options.write_moments || output_options.write_quadrature)
        {
            std::cout << "\nPer-location histories written to " << output_options.directory.string() << "/location_<i>_*.txt\n";
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
