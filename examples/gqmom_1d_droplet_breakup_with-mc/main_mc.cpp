#include "monte_carlo/moments.hpp"
#include "monte_carlo/output.hpp"
#include "particles_mc.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/output.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <cmath>
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
    const std::filesystem::path configuration_file = "configuration.qmom";
    const int step_count = 800;
    const double time_step = 0.01;

    auto configured_output = qmomkit::config::load_output_options(configuration_file);
    if (!configured_output)
    {
        std::cerr << configured_output.error() << '\n';
        return 1;
    }
    const qmomkit::OutputOptions output = std::move(configured_output).value();

    // The same sample the QMOMKit runs start from, here kept as droplets of equal
    // statistical weight, which breakup preserves.
    std::optional<gqmom_1d_droplet_breakup_with_mc::SamplingSetup> sampled =
        gqmom_1d_droplet_breakup_with_mc::sample_initial_distribution(configuration_file);
    if (!sampled)
    {
        return 1;
    }
    qmomkit::monte_carlo::Ensemble droplets = std::move(sampled->sampling.ensemble);
    const std::size_t initial_droplet_count = droplets.size();
    const double droplet_weight = sampled->sampling.number_density / static_cast<double>(initial_droplet_count);
    const std::size_t maximum_order = sampled->moments.moments.values().size() - 1;
    const std::vector<double> initial_moments = sampled->moments.moments.values();
    std::mt19937_64 breakup_generator(gqmom_1d_droplet_breakup_with_mc::breakup_seed);

    try
    {
        // The writer throws when a file cannot be written, and current_moments
        // when the companion cannot estimate the moments of the droplets.
        std::optional<qmomkit::MomentWriter> moment_writer;
        if (output.write_moments)
        {
            std::filesystem::create_directories(output.directory);
            moment_writer.emplace(
                qmomkit::monte_carlo::make_moment_history_writer(output.directory / "moment_history_mc.txt", initial_moments.size(), output));
        }
        const auto current_moments = [&]()
        {
            qmomkit::monte_carlo::OneDimensionalMomentEstimate estimate =
                qmomkit::monte_carlo::compute_ensemble_moments(droplets, droplet_weight * static_cast<double>(droplets.size()), maximum_order);
            if (!estimate.ok())
            {
                throw std::runtime_error(estimate.message);
            }
            return estimate.moments;
        };
        const auto write_output = [&](const double time)
        {
            if (moment_writer)
            {
                qmomkit::monte_carlo::write_moment_history_row(*moment_writer, time, current_moments(), droplets,
                                                               droplet_weight * static_cast<double>(droplets.size()), output);
            }
        };
        write_output(0.0);

        // Write every output interval of [output], and at the end without one.
        const int steps_per_output = output.interval ? static_cast<int>(std::lround(*output.interval / time_step)) : step_count;
        for (int step = 1; step <= step_count; ++step)
        {
            // Stochastic binary breakup over the step, the counterpart of the
            // source step of the QMOMKit runs.
            gqmom_1d_droplet_breakup_with_mc::apply_binary_breakup(droplets, time_step, breakup_generator);

            if (step % steps_per_output == 0 || step == step_count)
            {
                write_output(step * time_step);
            }
        }

        gqmom_1d_droplet_breakup_with_mc::print_moments(std::cout, "Monte Carlo run", initial_moments, current_moments().values());
        std::cout << "Droplets: " << initial_droplet_count << " -> " << droplets.size() << '\n';
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
