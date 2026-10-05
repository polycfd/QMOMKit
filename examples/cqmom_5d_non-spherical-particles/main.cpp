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
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
constexpr double initial_time = 0.0;
constexpr double final_time = 1.0;
constexpr double time_step = 0.005;
} // namespace

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <path-to-configuration.qmom>\n";
        return 1;
    }

    const std::filesystem::path configuration_file = argv[1];

    try
    {
        using namespace cqmom_5d_non_spherical_particles;

        // configuration.qmom defines a five-dimensional phase space (d, chi, u, v, w)
        // with three nodes per dimension, so the quadrature holds at most 3^5 = 243
        // terminal nodes. The callbacks built further below receive all of them in a
        // single call, packed node-major.
        auto configured = qmomkit::config::load_cqmom_prepared(configuration_file);
        if (!configured)
        {
            print_configuration_errors(configured.error());
            return 1;
        }
        auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured)
        {
            print_configuration_errors(output_configured.error());
            return 1;
        }
        auto evolution_configured = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        if (!evolution_configured)
        {
            print_configuration_errors(evolution_configured.error());
            return 1;
        }

        const std::optional<std::string> configuration_text = read_configuration_text(configuration_file);
        if (!configuration_text.has_value())
        {
            return 1;
        }
        auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(*configuration_text, configuration_file.string());
        if (!initial_condition)
        {
            print_configuration_errors(initial_condition.error());
            return 1;
        }

        // CQMOM consumes only the 606 moments its hierarchy requires, so the sampler
        // is asked for exactly those indices.
        const qmomkit::CqmomMomentLayout moment_layout = qmomkit::create_cqmom_moment_layout(configured.value().hierarchy());
        const qmomkit::RawMomentRequest required_moments{moment_layout.required_indices.to_multi_indices()};
        const auto sampled = sample_requested_initial_moments(std::move(initial_condition).value(), required_moments);
        if (!sampled.has_value())
        {
            return 1;
        }
        if (!sampled->ok())
        {
            std::cerr << sampled->sampling.message << sampled->moments.message << '\n';
            return 1;
        }

        auto sampled_case_result = qmomkit::config::with_initial_moments(
            configured.value(),
            qmomkit::create_cqmom_moment_data(moment_layout, qmomkit::make_selected_raw_moment_result(required_moments, sampled->moments.values)));
        if (!sampled_case_result)
        {
            print_configuration_errors(sampled_case_result.error());
            return 1;
        }

        qmomkit::CqmomPreparedCase cqmom_prepared_case = std::move(sampled_case_result).value();

        const ParticleModelParameters particle_parameters = particle_model_parameters();
        const BinaryFragmentationParameters breakage_parameters = binary_fragmentation_parameters();

        // carrier_samples is the host-side carrier flow, owned here by the driver.
        // The callback built just below receives a pointer to it, which is how a
        // host passes local, time-varying field data to a callback without a
        // global. It must outlive that callback. Only the newest level is known
        // at this point, so the first step runs without a reconstruction.
        CarrierFlowSamples carrier_samples;
        carrier_samples.current = carrier_flow_velocity(initial_time);
        carrier_samples.sample_time = initial_time;

        // Definition of the callbacks
        const qmomkit::PhaseSpaceVelocity velocity = particle_phase_space_velocity(particle_parameters, &carrier_samples);
        const qmomkit::CqmomMomentSource source = binary_breakage_source(breakage_parameters);

        const qmomkit::OutputOptions output_options = std::move(output_configured).value();
        const qmomkit::CqmomPhaseSpaceEvolutionOptions evolution_options = std::move(evolution_configured).value();

        // terminal_moment_request below lists the moments the final summary prints:
        // of the 606 moments this case needs, only the 12 lowest-order ones, since
        // printing all 606 would be useless. It is built from the moment layout
        // created at the top of main, which is still the correct layout here because
        // inserting the sampled initial moments changed their values, not which
        // moments the case needs.
        const qmomkit::RawMomentRequest terminal_moment_request{select_terminal_moment_indices(moment_layout.required_indices)};

        qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
        const qmomkit::CqmomQuadrature initial_quadrature = state.quadrature();

        if (output_options.write_quadrature || output_options.write_moments)
        {
            ensure_output_directory(output_options.directory);
        }
        std::optional<qmomkit::MomentWriter> moment_history_writer;
        if (output_options.write_moments)
        {
            qmomkit::MomentWriter::for_cqmom(output_options.directory / "initial_moments.txt", state, output_options)
                .write(initial_time, state, state.quadrature());
            moment_history_writer.emplace(qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history.txt", state, output_options));
            moment_history_writer->write(initial_time, state, state.quadrature());
        }
        std::optional<qmomkit::QuadratureWriter> quadrature_history_writer;
        if (output_options.write_quadrature)
        {
            qmomkit::QuadratureWriter::for_cqmom(output_options.directory / "initial_quadrature.txt", initial_quadrature.number_of_dimensions())
                .write(initial_time, initial_quadrature);
            quadrature_history_writer.emplace(
                qmomkit::QuadratureWriter::for_cqmom(output_options.directory / "quadrature_history.txt", initial_quadrature.number_of_dimensions()));
            quadrature_history_writer->write(initial_time, initial_quadrature);
        }

        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double next_output_time = initial_time + output_interval;
        double time = initial_time;

        while (time < final_time - 1.0e-12)
        {
            const double step = std::min(time_step, final_time - time);

            // The host advances its own field and publishes the new level, keeping
            // the level it replaces. Those two levels are all a coupled host has:
            // it cannot evaluate at QMOMKit's internal stage times, and it does
            // not know the field at time + step. The callback extends the line
            // through them, which is second-order accurate in time_step and so no
            // longer the limiting error of the second-order Strang-split step.
            // Freezing a single value here instead would make the coupling
            // first-order in u and v, the components in which this flow varies.
            if (time > initial_time)
            {
                carrier_samples.previous = carrier_samples.current;
                carrier_samples.current = carrier_flow_velocity(time);
                carrier_samples.sample_interval = time - carrier_samples.sample_time;
                carrier_samples.sample_time = time;
            }

            const qmomkit::CqmomPhaseSpaceEvolutionResult evolved =
                qmomkit::advance_cqmom_phase_space(state, time, step, velocity, source, evolution_options);
            if (!evolved.ok())
            {
                std::cerr << "CQMOM phase-space evolution failed: " << evolved.message << '\n';
                return 1;
            }
            state = evolved.state;
            time = evolved.final_time;
            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                if (quadrature_history_writer.has_value())
                {
                    quadrature_history_writer->write(time, state.quadrature());
                }
                if (moment_history_writer.has_value())
                {
                    moment_history_writer->write(time, state, state.quadrature());
                }
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
            }
        }

        // A failed reconstruction still returns moment values, so the status has
        // to be checked before they are reported as a result.
        const qmomkit::SelectedRawMomentResult initial_moments = qmomkit::reconstruct_raw_moments(initial_quadrature, terminal_moment_request);
        if (!initial_moments.ok())
        {
            throw std::runtime_error("Initial CQMOM moment reconstruction failed: " + initial_moments.message);
        }
        const qmomkit::SelectedRawMomentResult final_moments = qmomkit::reconstruct_raw_moments(state.quadrature(), terminal_moment_request);
        if (!final_moments.ok())
        {
            throw std::runtime_error("Final CQMOM moment reconstruction failed: " + final_moments.message);
        }
        print_result_summary(initial_moments, final_moments, output_options.directory, output_options.write_moments, output_options.write_quadrature);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
