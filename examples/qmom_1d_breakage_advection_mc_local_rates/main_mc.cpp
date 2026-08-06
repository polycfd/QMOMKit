#include "monte_carlo/evolution.hpp"
#include "particles_mc.hpp"
#include "physics.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/output.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.01;

    const std::optional<std::string> configuration_text = qmom_1d_breakage_advection_mc_local_rates::read_configuration_text(configuration_file);
    if (!configuration_text.has_value())
    {
        return 1;
    }

    auto configured = qmomkit::config::load_qmom(configuration_file);
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

    auto phase_space_evolution_configured = qmomkit::config::load_qmom_phase_space_evolution_options(configuration_file);
    if (!phase_space_evolution_configured)
    {
        std::cerr << phase_space_evolution_configured.error() << '\n';
        return 1;
    }

    const qmomkit::QmomCase qmom_case = std::move(configured).value();
    const qmomkit::OutputOptions output_options = std::move(output_configured).value();
    const qmomkit::QmomPhaseSpaceEvolutionOptions phase_space_evolution_options = std::move(phase_space_evolution_configured).value();
    if (phase_space_evolution_options.scheme != qmomkit::PhaseSpaceEvolutionScheme::strang)
    {
        std::cerr << "The local-rate Monte Carlo reference currently supports only the strang phase-space evolution scheme.\n";
        return 1;
    }

    auto initial_condition = qmomkit::config::parse_qmom_initial_condition(*configuration_text, configuration_file.string());
    if (!initial_condition)
    {
        std::cerr << initial_condition.error() << '\n';
        return 1;
    }

    const std::optional<qmom_1d_breakage_advection_mc_local_rates::SamplingSetup> sampled =
        qmom_1d_breakage_advection_mc_local_rates::sample_initial_distribution(std::move(initial_condition).value());
    if (!sampled.has_value())
    {
        return 1;
    }
    if (!sampled->sampling.ok())
    {
        std::cerr << sampled->sampling.message << '\n';
        return 1;
    }

    const std::size_t maximum_order = sampled->moments.moments.values().size() - 1;

    std::vector<qmom_1d_breakage_advection_mc_local_rates::Particle> particles;
    try
    {
        particles = qmom_1d_breakage_advection_mc_local_rates::make_particles(sampled->sampling.ensemble);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    const qmomkit::monte_carlo::OneDimensionalMomentEstimate& sampled_moments = sampled->moments;
    if (!sampled_moments.ok())
    {
        std::cerr << sampled_moments.message << '\n';
        return 1;
    }

    const double particle_weight = sampled->sampling.number_density / static_cast<double>(sampled->sampling.ensemble.size());
    if (!(particle_weight > 0.0) || !std::isfinite(particle_weight))
    {
        std::cerr << "Monte Carlo particle weight must be positive and finite.\n";
        return 1;
    }
    std::vector<double> final_moments =
        qmom_1d_breakage_advection_mc_local_rates::compute_particle_moments(particles, particle_weight, maximum_order);
    const std::vector<double> initial_moments = final_moments;
    const std::vector<double> sampled_initial_moments = sampled_moments.moments.values();
    double initial_moment_mismatch = 0.0;
    for (std::size_t index = 0; index < std::min(initial_moments.size(), sampled_initial_moments.size()); ++index)
    {
        initial_moment_mismatch = std::max(initial_moment_mismatch, std::abs(initial_moments[index] - sampled_initial_moments[index]));
    }

    std::size_t next_particle_id = particles.size();
    std::mt19937_64 breakage_generator(qmom_1d_breakage_advection_mc_local_rates::breakage_seed);
    std::mt19937_64 nucleation_generator(qmom_1d_breakage_advection_mc_local_rates::nucleation_seed);
    const qmomkit::PhaseSpaceVelocity velocity = qmom_1d_breakage_advection_mc_local_rates::advection_velocity();
    const qmomkit::AdaptiveRkOptions advection_options = phase_space_evolution_options.node_options.adaptive_options;

    std::vector<std::tuple<double, std::size_t, double>> particle_history;
    std::vector<double> mc_times{initial_time};
    std::vector<std::vector<double>> mc_moment_history{final_moments};
    std::size_t total_breakage_events = 0;
    std::size_t total_nucleation_particles = 0;
    double expected_cumulative_nucleation_number = 0.0;
    double realized_cumulative_nucleation_number = 0.0;
    double max_source_first_moment_balance_error = 0.0;
    if (output_options.write_moments)
    {
        for (const qmom_1d_breakage_advection_mc_local_rates::Particle& particle : particles)
        {
            particle_history.emplace_back(initial_time, particle.id, particle.diameter);
        }
    }

    std::optional<qmomkit::MomentWriter> moment_history_writer;
    std::optional<std::ofstream> diagnostics_writer;
    try
    {
        if (output_options.write_moments)
        {
            std::filesystem::create_directories(output_options.directory);
            moment_history_writer.emplace(qmomkit::MomentWriter::for_qmom(output_options.directory / "moment_history_mc.txt", final_moments.size()));
            moment_history_writer->write(initial_time, qmomkit::MomentSet(final_moments));

            diagnostics_writer.emplace(output_options.directory / "mc_diagnostics.txt");
            if (!*diagnostics_writer)
            {
                std::cerr << "Could not open Monte Carlo diagnostics output file.\n";
                return 1;
            }
            *diagnostics_writer << "# time particle_count total_number_density first_moment breakage_event_count nucleation_particle_count "
                                   "expected_cumulative_nucleation_number realized_cumulative_nucleation_number nucleation_m0_increment\n";
            *diagnostics_writer << std::setprecision(17) << initial_time << ' ' << particles.size() << ' ' << final_moments[0] << ' '
                                << final_moments[1] << ' ' << total_breakage_events << ' ' << total_nucleation_particles << ' '
                                << expected_cumulative_nucleation_number << ' ' << realized_cumulative_nucleation_number << ' '
                                << realized_cumulative_nucleation_number << '\n';
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    const double output_interval = output_options.interval.value_or(final_time - initial_time);
    double time = initial_time;
    double next_output_time = initial_time + output_interval;
    while (time < final_time - 1.0e-12)
    {
        const double current_step = std::min(time_step, final_time - time);
        const double half_step = 0.5 * current_step;

        qmomkit::monte_carlo::EnsembleEvolutionResult first_advection = qmomkit::monte_carlo::evolve_ensemble(
            qmom_1d_breakage_advection_mc_local_rates::to_ensemble(particles), velocity, time, time + half_step, advection_options);
        if (!first_advection.ok())
        {
            std::cerr << "First advection particle half-step failed: " << first_advection.message << '\n';
            return 1;
        }

        try
        {
            qmom_1d_breakage_advection_mc_local_rates::update_particle_diameters(particles, first_advection.ensemble);
        }
        catch (const std::exception& error)
        {
            std::cerr << error.what() << '\n';
            return 1;
        }

        const double first_moment_before_source =
            qmom_1d_breakage_advection_mc_local_rates::compute_particle_moments(particles, particle_weight, 1)[1];
        const std::size_t breakage_events =
            qmom_1d_breakage_advection_mc_local_rates::apply_binary_breakage(particles, current_step, breakage_generator, next_particle_id);
        const std::size_t nucleation_particles = qmom_1d_breakage_advection_mc_local_rates::apply_nucleation(particles, current_step, particle_weight,
                                                                                                             nucleation_generator, next_particle_id);
        total_breakage_events += breakage_events;
        total_nucleation_particles += nucleation_particles;
        expected_cumulative_nucleation_number += qmom_1d_breakage_advection_mc_local_rates::nucleation_rate() * current_step;
        realized_cumulative_nucleation_number = static_cast<double>(total_nucleation_particles) * particle_weight;
        const double first_moment_after_source =
            qmom_1d_breakage_advection_mc_local_rates::compute_particle_moments(particles, particle_weight, 1)[1];
        const double realized_nucleation_first_moment =
            static_cast<double>(nucleation_particles) * particle_weight * qmom_1d_breakage_advection_mc_local_rates::nucleation_diameter();
        max_source_first_moment_balance_error =
            std::max(max_source_first_moment_balance_error,
                     std::abs(first_moment_after_source - first_moment_before_source - realized_nucleation_first_moment));

        qmomkit::monte_carlo::EnsembleEvolutionResult second_advection = qmomkit::monte_carlo::evolve_ensemble(
            qmom_1d_breakage_advection_mc_local_rates::to_ensemble(particles), velocity, time + half_step, time + current_step, advection_options);
        if (!second_advection.ok())
        {
            std::cerr << "Second advection particle half-step failed: " << second_advection.message << '\n';
            return 1;
        }

        try
        {
            qmom_1d_breakage_advection_mc_local_rates::update_particle_diameters(particles, second_advection.ensemble);
        }
        catch (const std::exception& error)
        {
            std::cerr << error.what() << '\n';
            return 1;
        }

        time += current_step;

        if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
        {
            final_moments = qmom_1d_breakage_advection_mc_local_rates::compute_particle_moments(particles, particle_weight, maximum_order);
            mc_times.push_back(time);
            mc_moment_history.push_back(final_moments);
            if (moment_history_writer.has_value())
            {
                try
                {
                    moment_history_writer->write(time, qmomkit::MomentSet(final_moments));
                    for (const qmom_1d_breakage_advection_mc_local_rates::Particle& particle : particles)
                    {
                        particle_history.emplace_back(time, particle.id, particle.diameter);
                    }
                    *diagnostics_writer << std::setprecision(17) << time << ' ' << particles.size() << ' ' << final_moments[0] << ' '
                                        << final_moments[1] << ' ' << total_breakage_events << ' ' << total_nucleation_particles << ' '
                                        << expected_cumulative_nucleation_number << ' ' << realized_cumulative_nucleation_number << ' '
                                        << realized_cumulative_nucleation_number << '\n';
                }
                catch (const std::exception& error)
                {
                    std::cerr << error.what() << '\n';
                    return 1;
                }
            }
            while (next_output_time <= time + 1.0e-12)
            {
                next_output_time += output_interval;
            }
        }
    }

    const std::filesystem::path particle_history_file = output_options.directory / "particle_history.txt";
    try
    {
        if (output_options.write_moments)
        {
            std::ofstream particle_output(particle_history_file);
            if (!particle_output)
            {
                std::cerr << "Could not open particle history output file: " << particle_history_file << '\n';
                return 1;
            }
            particle_output << "# time particle_id diameter\n";
            particle_output << std::setprecision(17);
            for (const auto& [sample_time, particle_id, diameter] : particle_history)
            {
                particle_output << sample_time << ' ' << particle_id << ' ' << diameter << '\n';
            }
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    bool qmom_history_available = false;
    bool time_grids_match = false;
    double initial_max_absolute_difference = std::numeric_limits<double>::quiet_NaN();
    double final_max_absolute_difference = std::numeric_limits<double>::quiet_NaN();
    double final_zeroth_moment_difference = std::numeric_limits<double>::quiet_NaN();
    double final_first_moment_difference = std::numeric_limits<double>::quiet_NaN();
    double max_absolute_difference = std::numeric_limits<double>::quiet_NaN();
    double max_relative_difference = std::numeric_limits<double>::quiet_NaN();

    std::ifstream qmom_history_input(output_options.directory / "moment_history.txt");
    if (qmom_history_input)
    {
        std::vector<double> qmom_times;
        std::vector<std::vector<double>> qmom_moment_history;
        std::string line;
        while (std::getline(qmom_history_input, line))
        {
            if (line.empty() || line[0] == '#')
            {
                continue;
            }

            std::istringstream row(line);
            double row_time = 0.0;
            row >> row_time;
            std::vector<double> row_moments;
            double moment = 0.0;
            while (row >> moment)
            {
                row_moments.push_back(moment);
            }
            if (!row_moments.empty())
            {
                qmom_times.push_back(row_time);
                qmom_moment_history.push_back(row_moments);
            }
        }

        qmom_history_available = !qmom_moment_history.empty() && qmom_moment_history.size() == mc_moment_history.size();
        time_grids_match = qmom_history_available;
        max_absolute_difference = 0.0;
        max_relative_difference = 0.0;
        if (qmom_history_available)
        {
            for (std::size_t row = 0; row < qmom_moment_history.size(); ++row)
            {
                if (std::abs(qmom_times[row] - mc_times[row]) > 1.0e-12 || qmom_moment_history[row].size() != mc_moment_history[row].size())
                {
                    time_grids_match = false;
                    qmom_history_available = false;
                    break;
                }
                for (std::size_t index = 0; index < qmom_moment_history[row].size(); ++index)
                {
                    const double difference = std::abs(qmom_moment_history[row][index] - mc_moment_history[row][index]);
                    const double scale = std::max({1.0, std::abs(qmom_moment_history[row][index]), std::abs(mc_moment_history[row][index])});
                    max_absolute_difference = std::max(max_absolute_difference, difference);
                    max_relative_difference = std::max(max_relative_difference, difference / scale);
                    if (row == 0)
                    {
                        initial_max_absolute_difference =
                            std::isnan(initial_max_absolute_difference) ? difference : std::max(initial_max_absolute_difference, difference);
                    }
                    if (row + 1 == qmom_moment_history.size())
                    {
                        final_max_absolute_difference =
                            std::isnan(final_max_absolute_difference) ? difference : std::max(final_max_absolute_difference, difference);
                        if (index == 0)
                        {
                            final_zeroth_moment_difference = difference;
                        }
                        if (index == 1)
                        {
                            final_first_moment_difference = difference;
                        }
                    }
                }
            }
        }
    }

    std::cout << std::setprecision(12);
    qmom_1d_breakage_advection_mc_local_rates::print_common_run_summary(
        std::cout, "Monte Carlo particle reference", configuration_file, output_options.directory,
        qmom_1d_breakage_advection_mc_local_rates::initialization_sample_count, qmom_1d_breakage_advection_mc_local_rates::initialization_seed,
        initial_time, final_time, time_step, output_interval);
    std::cout << "Representation: sampled particles with constant statistical weight\n";
    std::cout << "Configured QMOM node count: " << qmom_case.node_count << '\n';
    std::cout << "Particle weight: " << particle_weight << '\n';
    std::cout << "Breakage seed: " << qmom_1d_breakage_advection_mc_local_rates::breakage_seed << '\n';
    std::cout << "Nucleation seed: " << qmom_1d_breakage_advection_mc_local_rates::nucleation_seed << '\n';
    std::cout << "Fragmentation rate: " << qmom_1d_breakage_advection_mc_local_rates::breakage_event_rate() << '\n';
    std::cout << "Nucleation rate: " << qmom_1d_breakage_advection_mc_local_rates::nucleation_rate() << '\n';
    std::cout << "Nucleation diameter: " << qmom_1d_breakage_advection_mc_local_rates::nucleation_diameter() << '\n';
    std::cout << "Initial sampled-moment mismatch: " << initial_moment_mismatch << '\n';

    qmom_1d_breakage_advection_mc_local_rates::print_moments(std::cout, "Initial moments", initial_moments);
    qmom_1d_breakage_advection_mc_local_rates::print_moments(std::cout, "Final moments", final_moments);

    std::cout << "\nFinal ensemble at t = " << time << ":\n";
    std::cout << "initial_particle_count = " << qmom_1d_breakage_advection_mc_local_rates::initialization_sample_count << '\n';
    std::cout << "final_particle_count = " << particles.size() << '\n';
    std::cout << "breakage_event_count = " << total_breakage_events << '\n';
    std::cout << "nucleation_particle_count = " << total_nucleation_particles << '\n';
    std::cout << "expected_cumulative_nucleation_number = " << expected_cumulative_nucleation_number << '\n';
    std::cout << "realized_cumulative_nucleation_number = " << realized_cumulative_nucleation_number << '\n';
    std::cout << "max_source_first_moment_balance_error = " << max_source_first_moment_balance_error << '\n';

    if (qmom_history_available && time_grids_match)
    {
        std::cout << "\nQMOM comparison from " << (output_options.directory / "moment_history.txt") << ":\n";
        std::cout << "initial_max_absolute_moment_difference = " << initial_max_absolute_difference << '\n';
        std::cout << "final_max_absolute_moment_difference = " << final_max_absolute_difference << '\n';
        std::cout << "final_zeroth_moment_difference = " << final_zeroth_moment_difference << '\n';
        std::cout << "final_first_moment_difference = " << final_first_moment_difference << '\n';
        std::cout << "max_absolute_moment_difference = " << max_absolute_difference << '\n';
        std::cout << "max_relative_moment_difference = " << max_relative_difference << '\n';
    }
    else
    {
        std::cout << "\nQMOM comparison skipped. Run qmom_1d_breakage_advection_mc_local_rates_example first in the same build directory.\n";
    }

    if (output_options.write_moments)
    {
        std::cout << "\nCompleted Monte Carlo particle workflow. Outputs written to " << output_options.directory << '\n';
    }
    else
    {
        std::cout << "\nCompleted Monte Carlo particle workflow. Output disabled by configuration.\n";
    }

    const double expected_nucleation_particles = expected_cumulative_nucleation_number / particle_weight;
    const double nucleation_number_tolerance = 6.0 * std::sqrt(expected_nucleation_particles) * particle_weight;
    if (total_nucleation_particles == 0 ||
        std::abs(realized_cumulative_nucleation_number - expected_cumulative_nucleation_number) > nucleation_number_tolerance)
    {
        std::cerr << "Realized nucleation count is outside the expected fixed-seed Poisson sanity interval.\n";
        return 1;
    }

    return 0;
}
