#include "droplet_model.hpp"
#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "monte_carlo/output.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.01;

    // The hierarchy decides which moments the sample estimates and the history
    // records; the node-advection options set the accuracy of the droplet
    // trajectories.
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
    const qmomkit::AdaptiveRkOptions advection_options = configured_options.value().node_options.adaptive_options;
    const qmomkit::OutputOptions output = std::move(configured_output).value();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(configured_case.value().hierarchy());
    const qmomkit::RawMomentRequest layout_moments{layout.required_indices.to_multi_indices()};

    // The same sample the QMOMKit run starts from, here kept as droplets of equal
    // statistical weight.
    const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
        cqmom_3d_droplets_condensation_nucleation::sample_initial_distribution(configuration_file, configured_case.value().hierarchy());
    if (!sampled)
    {
        return 1;
    }
    qmomkit::monte_carlo::Ensemble ensemble = sampled->sampling.ensemble;
    const std::size_t initial_droplet_count = ensemble.size();
    double number_density = sampled->sampling.number_density;
    const cqmom_3d_droplets_condensation_nucleation::NucleationModel nucleation_model =
        cqmom_3d_droplets_condensation_nucleation::make_nucleation_model(number_density);

    const qmomkit::PhaseSpaceVelocity velocity = cqmom_3d_droplets_condensation_nucleation::droplet_phase_space_velocity();
    const qmomkit::RawMomentRequest printed_moments = cqmom_3d_droplets_condensation_nucleation::printed_moments();
    const auto estimate_moments = [&](const qmomkit::RawMomentRequest& request)
    {
        const qmomkit::monte_carlo::SelectedEnsembleMomentEstimate estimate =
            qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, number_density, request);
        if (!estimate.ok())
        {
            throw std::runtime_error(estimate.message);
        }
        return estimate.values;
    };

    try
    {
        // The estimates and the writer throw when they fail.
        const std::vector<double> initial_moments = estimate_moments(printed_moments);
        std::optional<qmomkit::MomentWriter> moment_writer;
        if (output.write_moments)
        {
            std::filesystem::create_directories(output.directory);
            moment_writer.emplace(qmomkit::monte_carlo::make_moment_history_writer(output.directory / "moment_history_mc.txt", layout, output));
        }
        const auto write_output = [&](const double time)
        {
            if (moment_writer)
            {
                const qmomkit::CqmomMomentData moments = qmomkit::create_cqmom_moment_data(
                    layout, qmomkit::make_selected_raw_moment_result(layout_moments, estimate_moments(layout_moments)));
                qmomkit::monte_carlo::write_moment_history_row(*moment_writer, time, moments, ensemble, number_density, output);
            }
        };
        write_output(initial_time);

        const double output_interval = output.interval.value_or(final_time - initial_time);
        double next_output_time = initial_time + output_interval;
        double time = initial_time;
        std::size_t birth_event_index = 0;
        while (time < final_time - 1.0e-12)
        {
            // Steps end at the start and the end of nucleation, as in the QMOMKit
            // run. While nucleation is active, Strang splitting places the births
            // between two advection half-steps.
            const double step =
                cqmom_3d_droplets_condensation_nucleation::limit_step_to_nucleation_interval(time, std::min(time_step, final_time - time));
            if (cqmom_3d_droplets_condensation_nucleation::nucleation_is_active(time))
            {
                const double half_time = time + 0.5 * step;
                qmomkit::monte_carlo::EnsembleEvolutionResult first_advection =
                    qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, half_time, advection_options);
                if (!first_advection.ok())
                {
                    std::cerr << "First advection half-step failed: " << first_advection.message << '\n';
                    return 1;
                }
                const qmomkit::monte_carlo::PopulationBirthResult birth = cqmom_3d_droplets_condensation_nucleation::inject_nucleated_droplets(
                    nucleation_model, first_advection.ensemble, number_density, step,
                    cqmom_3d_droplets_condensation_nucleation::initialization_seed + static_cast<unsigned long long>(birth_event_index++));
                if (!birth.ok())
                {
                    std::cerr << "Droplet nucleation failed: " << birth.message << '\n';
                    return 1;
                }
                number_density = birth.number_density;
                qmomkit::monte_carlo::EnsembleEvolutionResult second_advection =
                    qmomkit::monte_carlo::evolve_ensemble(birth.ensemble, velocity, half_time, time + step, advection_options);
                if (!second_advection.ok())
                {
                    std::cerr << "Second advection half-step failed: " << second_advection.message << '\n';
                    return 1;
                }
                ensemble = std::move(second_advection.ensemble);
            }
            else
            {
                qmomkit::monte_carlo::EnsembleEvolutionResult advection =
                    qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + step, advection_options);
                if (!advection.ok())
                {
                    std::cerr << "Droplet advection failed: " << advection.message << '\n';
                    return 1;
                }
                ensemble = std::move(advection.ensemble);
            }
            time += step;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                write_output(time);
                next_output_time += output_interval;
            }
        }

        cqmom_3d_droplets_condensation_nucleation::print_moments(std::cout, "Monte Carlo run", initial_moments, estimate_moments(printed_moments));
        std::cout << "Droplets: " << initial_droplet_count << " -> " << ensemble.size() << '\n';
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
