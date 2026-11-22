#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "monte_carlo/output.hpp"
#include "particles_mc.hpp"
#include "physics.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/output.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

int main()
{
    // The particles do not depend on the conditioning order, so this run reads its
    // initial distribution and its output settings from velocity_first.qmom,
    // which size_first.qmom shares.
    const std::filesystem::path configuration_file = "velocity_first.qmom";
    const double initial_time = 0.0;
    const double final_time = 0.01;
    const double time_step = 5.0e-5;

    auto configured_output = qmomkit::config::load_output_options(configuration_file);
    if (!configured_output)
    {
        std::cerr << configured_output.error() << '\n';
        return 1;
    }
    const qmomkit::OutputOptions output = std::move(configured_output).value();

    // The same sample the QMOMKit runs start from, here kept as particles of
    // equal statistical weight.
    const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
        cqmom_2d_hyperbolic_drag_dispersion::sample_initial_distribution(configuration_file);
    if (!sampled)
    {
        return 1;
    }
    qmomkit::monte_carlo::Ensemble ensemble = sampled->sampling.ensemble;
    const double particle_weight = sampled->sampling.number_density / static_cast<double>(ensemble.size());
    const std::vector<std::size_t> moment_extents = sampled->moments.moments.extents;
    std::mt19937_64 dispersion_generator(cqmom_2d_hyperbolic_drag_dispersion::dispersion_seed);

    const qmomkit::PhaseSpaceVelocity velocity = cqmom_2d_hyperbolic_drag_dispersion::drag_velocity();
    const qmomkit::AdaptiveRkOptions drag_options;
    const qmomkit::RawMomentRequest printed_moments = cqmom_2d_hyperbolic_drag_dispersion::printed_moments();
    const auto estimate_printed_moments = [&]()
    {
        const qmomkit::monte_carlo::SelectedEnsembleMomentEstimate estimate = qmomkit::monte_carlo::compute_requested_ensemble_moments(
            ensemble, particle_weight * static_cast<double>(ensemble.size()), printed_moments);
        if (!estimate.ok())
        {
            throw std::runtime_error(estimate.message);
        }
        return estimate.values;
    };

    try
    {
        // The estimates and the writer throw when they fail.
        const std::vector<double> initial_moments = estimate_printed_moments();
        std::optional<qmomkit::MomentWriter> moment_writer;
        if (output.write_moments)
        {
            std::filesystem::create_directories(output.directory);
            moment_writer.emplace(
                qmomkit::monte_carlo::make_moment_history_writer(output.directory / "moment_history_mc.txt", moment_extents, output));
        }
        const auto write_output = [&](const double time)
        {
            if (moment_writer)
            {
                qmomkit::monte_carlo::write_moment_history_row(
                    *moment_writer, time, cqmom_2d_hyperbolic_drag_dispersion::compute_particle_moments(ensemble, particle_weight, moment_extents),
                    ensemble, particle_weight * static_cast<double>(ensemble.size()), output);
            }
        };
        write_output(initial_time);

        const double output_interval = output.interval.value_or(final_time - initial_time);
        double next_output_time = initial_time + output_interval;
        double time = initial_time;
        while (time < final_time - 1.0e-12)
        {
            const double step = std::min(time_step, final_time - time);
            const double half_step = 0.5 * step;

            // Strang splitting, as in the QMOMKit runs: a half step of drag, the
            // dispersion over the whole step, and a second half step of drag.
            qmomkit::monte_carlo::EnsembleEvolutionResult first_drag =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + half_step, drag_options);
            if (!first_drag.ok())
            {
                std::cerr << "First drag half step of the particles failed: " << first_drag.message << '\n';
                return 1;
            }
            ensemble = std::move(first_drag.ensemble);

            cqmom_2d_hyperbolic_drag_dispersion::apply_dispersion(ensemble, step, dispersion_generator);

            qmomkit::monte_carlo::EnsembleEvolutionResult second_drag =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time + half_step, time + step, drag_options);
            if (!second_drag.ok())
            {
                std::cerr << "Second drag half step of the particles failed: " << second_drag.message << '\n';
                return 1;
            }
            ensemble = std::move(second_drag.ensemble);
            time += step;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                write_output(time);
                next_output_time += output_interval;
            }
        }

        cqmom_2d_hyperbolic_drag_dispersion::print_moments(std::cout, "Monte Carlo run", initial_moments, estimate_printed_moments());
        if (moment_writer)
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
