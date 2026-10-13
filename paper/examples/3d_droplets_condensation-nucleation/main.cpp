#include "droplet_model.hpp"
#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"

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

    // The CQMOM run starts from the moments of a finite sample of the configured
    // distribution, and the Monte Carlo reference advances that same sample.
    const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
        paper_3d_droplets_condensation_nucleation::sample_initial_distribution(configuration_file, configured_case.value().hierarchy());
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
    qmomkit::monte_carlo::Ensemble ensemble = sampled->sampling.ensemble;
    double number_density = sampled->sampling.number_density;

    const qmomkit::PhaseSpaceVelocity velocity = paper_3d_droplets_condensation_nucleation::droplet_phase_space_velocity();
    const qmomkit::AdaptiveRkOptions& mc_advection_options = options.node_options.adaptive_options;
    // Nucleation creates droplets at the edge of phase space, so it is a boundary
    // contribution rather than a volumetric source. The composed step still
    // requires a source callback, hence this zero-rate one.
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };
    // Both histories hold the moments the hierarchy needs: the QMOMKit history
    // those the state tracks, the Monte Carlo history their estimates from the
    // droplets.
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(cqmom_case.hierarchy());
    const qmomkit::RawMomentRequest layout_moments{layout.required_indices.to_multi_indices()};
    const auto estimate_mc_moments = [&]()
    {
        const qmomkit::monte_carlo::SelectedEnsembleMomentEstimate estimate =
            qmomkit::monte_carlo::compute_requested_ensemble_moments(ensemble, number_density, layout_moments);
        if (!estimate.ok())
        {
            throw std::runtime_error("Monte Carlo moment estimation failed: " + estimate.message);
        }
        return qmomkit::create_cqmom_moment_data(layout, qmomkit::make_selected_raw_moment_result(layout_moments, estimate.values));
    };

    try
    {
        // initial_state() throws when the sampled moments admit no quadrature, and
        // the writers and the Monte Carlo estimates throw when they fail.
        qmomkit::CqmomState state = cqmom_case.initial_state();

        // The nucleation flux is scaled by the initial number density M[0,0,0] of
        // the quadrature.
        const qmomkit::SelectedRawMomentResult initial_number_density =
            qmomkit::reconstruct_raw_moments(state.quadrature(), qmomkit::RawMomentRequest{std::vector<qmomkit::MultiIndex>{{0, 0, 0}}});
        if (!initial_number_density.ok())
        {
            std::cerr << "Initial number density reconstruction failed: " << initial_number_density.message << '\n';
            return 1;
        }
        const paper_3d_droplets_condensation_nucleation::NucleationModel nucleation_model =
            paper_3d_droplets_condensation_nucleation::make_nucleation_model(initial_number_density.values.front());
        const qmomkit::CqmomMomentBoundaryContribution nucleation =
            paper_3d_droplets_condensation_nucleation::nucleation_boundary_contribution(nucleation_model);

        std::optional<qmomkit::MomentWriter> moment_writer;
        std::optional<qmomkit::MomentWriter> mc_moment_writer;
        std::optional<qmomkit::QuadratureWriter> quadrature_writer;
        if (output.write_moments || output.write_quadrature)
        {
            std::filesystem::create_directories(output.directory);
        }
        if (output.write_moments)
        {
            moment_writer.emplace(qmomkit::MomentWriter::for_cqmom(output.directory / "moment_history.txt", state, output));
            mc_moment_writer.emplace(qmomkit::MomentWriter::for_cqmom(output.directory / "moment_history_mc.txt", layout));
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
                mc_moment_writer->write(time, estimate_mc_moments());
            }
            if (quadrature_writer)
            {
                quadrature_writer->write(time, state.quadrature());
            }
        };
        write_output(initial_time);

        // Each birth event samples its droplets with its own seed.
        std::size_t birth_event_index = 0;
        double time = initial_time;
        while (time < final_time)
        {
            // Steps end at the start and the end of nucleation. While nucleation is
            // active, Strang splitting places the births between two advection
            // half-steps, for the quadrature and for the droplets alike.
            const double step =
                paper_3d_droplets_condensation_nucleation::limit_step_to_nucleation_interval(time, std::min(time_step, final_time - time));
            const bool nucleation_active = paper_3d_droplets_condensation_nucleation::nucleation_is_active(time);
            qmomkit::CqmomPhaseSpaceEvolutionResult result =
                nucleation_active ? qmomkit::advance_cqmom_phase_space(state, time, step, velocity, zero_source, nucleation, options)
                                  : qmomkit::advance_cqmom_advection(state, time, step, velocity, options);
            if (!result.ok())
            {
                std::cerr << "CQMOM phase-space evolution failed: " << result.message << '\n';
                return 1;
            }
            state = std::move(result.state);

            if (nucleation_active)
            {
                const double half_time = time + 0.5 * step;
                qmomkit::monte_carlo::EnsembleEvolutionResult first_advection =
                    qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, half_time, mc_advection_options);
                if (!first_advection.ok())
                {
                    std::cerr << "First advection half-step failed: " << first_advection.message << '\n';
                    return 1;
                }
                const qmomkit::monte_carlo::PopulationBirthResult birth = paper_3d_droplets_condensation_nucleation::inject_nucleated_droplets(
                    nucleation_model, first_advection.ensemble, number_density, step,
                    paper_3d_droplets_condensation_nucleation::initialization_seed + static_cast<unsigned long long>(birth_event_index++));
                if (!birth.ok())
                {
                    std::cerr << "Droplet nucleation failed: " << birth.message << '\n';
                    return 1;
                }
                number_density = birth.number_density;
                qmomkit::monte_carlo::EnsembleEvolutionResult second_advection =
                    qmomkit::monte_carlo::evolve_ensemble(birth.ensemble, velocity, half_time, time + step, mc_advection_options);
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
                    qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + step, mc_advection_options);
                if (!advection.ok())
                {
                    std::cerr << "Droplet advection failed: " << advection.message << '\n';
                    return 1;
                }
                ensemble = std::move(advection.ensemble);
            }

            time += step;
            write_output(time);
        }

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
