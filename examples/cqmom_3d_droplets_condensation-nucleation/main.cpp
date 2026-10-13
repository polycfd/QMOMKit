#include "droplet_model.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
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
    const double final_time = 1.0;
    const double time_step = 0.01;

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
    const qmomkit::CqmomPhaseSpaceEvolutionOptions options = std::move(configured_options).value();
    const qmomkit::OutputOptions output = std::move(configured_output).value();

    // Starting from the moments of the sample rather than from the configured
    // analytical moments is what makes the two runs comparable: both begin from
    // the same finite population. The case keeps its hierarchy and settings.
    const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
        cqmom_3d_droplets_condensation_nucleation::sample_initial_distribution(configuration_file, configured_case.value().hierarchy());
    if (!sampled)
    {
        return 1;
    }
    auto sampled_case = qmomkit::config::with_initial_moments(configured_case.value(), sampled->moments.moments);
    if (!sampled_case)
    {
        std::cerr << sampled_case.error() << '\n';
        return 1;
    }
    const qmomkit::CqmomPreparedCase cqmom_case = std::move(sampled_case).value();

    const qmomkit::PhaseSpaceVelocity velocity = cqmom_3d_droplets_condensation_nucleation::droplet_phase_space_velocity();
    // Nucleation creates droplets at the edge of phase space, so it is a boundary
    // contribution rather than a volumetric source. The composed step still
    // requires a source callback, hence this zero-rate one.
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };
    const qmomkit::RawMomentRequest printed_moments = cqmom_3d_droplets_condensation_nucleation::printed_moments();

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
        // The nucleation flux is scaled by the initial number density M[0,0,0], the
        // first printed moment.
        const qmomkit::CqmomMomentBoundaryContribution nucleation = cqmom_3d_droplets_condensation_nucleation::nucleation_boundary_contribution(
            cqmom_3d_droplets_condensation_nucleation::make_nucleation_model(initial_moments.values.front()));

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
            // Steps end at the start and the end of nucleation. While nucleation is
            // active, the step composes advection with the boundary contribution;
            // otherwise there is nothing to compose, and advection alone is used.
            const double step =
                cqmom_3d_droplets_condensation_nucleation::limit_step_to_nucleation_interval(time, std::min(time_step, final_time - time));
            qmomkit::CqmomPhaseSpaceEvolutionResult result =
                cqmom_3d_droplets_condensation_nucleation::nucleation_is_active(time)
                    ? qmomkit::advance_cqmom_phase_space(state, time, step, velocity, zero_source, nucleation, options)
                    : qmomkit::advance_cqmom_advection(state, time, step, velocity, options);
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
        cqmom_3d_droplets_condensation_nucleation::print_moments(std::cout, "QMOMKit run", initial_moments.values, final_moments.values);
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
