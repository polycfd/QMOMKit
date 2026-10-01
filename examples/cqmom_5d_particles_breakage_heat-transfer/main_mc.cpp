#include "monte_carlo/evolution.hpp"
#include "monte_carlo/output.hpp"
#include "physics.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/output.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <random>
#include <stdexcept>

int main()
{
    using namespace cqmom_5d_particles_breakage_heat_transfer;
    constexpr unsigned long long breakup_seed = 321;
    constexpr double initial_time = 0.0;
    constexpr double final_time = 0.1;
    constexpr double time_step = 0.001;
    try
    {
        const std::filesystem::path configuration_file = "configuration.qmom";
        const auto configured = qmomkit::config::load_cqmom_prepared(configuration_file);
        const auto output = qmomkit::config::load_output_options(configuration_file);
        const auto text = read_configuration_text(configuration_file);
        if (!configured || !output || !text)
        {
            return 1;
        }
        const auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(*text, configuration_file.string());
        if (!initial_condition)
        {
            return 1;
        }
        const auto sampled = sample_initial_population(initial_condition.value(), configured.value().hierarchy());
        if (!sampled || !sampled->ok())
        {
            return 1;
        }
        qmomkit::monte_carlo::Ensemble ensemble = sampled->sampling.ensemble;
        const double particle_weight = sampled->sampling.number_density / static_cast<double>(ensemble.size());
        qmomkit::RawMomentTensor moments = ensemble_moments(ensemble, particle_weight, sampled->moments.moments.extents);
        if (maximum_tensor_difference(sampled->moments.moments, moments) != 0.0)
        {
            throw std::runtime_error("The sampled CQMOM and Monte Carlo initial moments differ.");
        }
        const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(configured.value().hierarchy());
        const qmomkit::RawMomentRequest request{layout.required_indices.to_multi_indices()};
        const qmomkit::SelectedRawMomentResult initial_moments = select_moments(moments, request);

        std::optional<qmomkit::MomentWriter> moment_history_writer;
        if (output.value().write_moments)
        {
            ensure_output_directory(output.value().directory);
            moment_history_writer.emplace(qmomkit::monte_carlo::make_moment_history_writer(output.value().directory / "moment_history_mc.txt",
                                                                                           moments.extents, request.indices, output.value()));
            qmomkit::monte_carlo::write_moment_history_row(*moment_history_writer, initial_time, moments, ensemble,
                                                           particle_weight * static_cast<double>(ensemble.size()), output.value());
        }

        std::mt19937_64 generator(breakup_seed);
        const qmomkit::AdaptiveRkOptions advection;
        const double output_interval = output.value().interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output = initial_time + output_interval;
        while (time < final_time - 1.0e-14)
        {
            const double step = std::min(time_step, final_time - time);
            auto first = qmomkit::monte_carlo::evolve_ensemble(ensemble, particle_velocity(), time, time + 0.5 * step, advection);
            if (!first.ok())
            {
                throw std::runtime_error(first.message);
            }
            ensemble = std::move(first.ensemble);
            apply_breakup(ensemble, step, generator);
            auto second = qmomkit::monte_carlo::evolve_ensemble(ensemble, particle_velocity(), time + 0.5 * step, time + step, advection);
            if (!second.ok())
            {
                throw std::runtime_error(second.message);
            }
            ensemble = std::move(second.ensemble);
            time += step;
            if (time + 1.0e-14 >= next_output || time >= final_time - 1.0e-14)
            {
                moments = ensemble_moments(ensemble, particle_weight, sampled->moments.moments.extents);
                if (moment_history_writer.has_value())
                {
                    qmomkit::monte_carlo::write_moment_history_row(*moment_history_writer, time, moments, ensemble,
                                                                   particle_weight * static_cast<double>(ensemble.size()), output.value());
                }
                while (next_output <= time + 1.0e-14)
                {
                    next_output += output_interval;
                }
            }
        }
        const qmomkit::SelectedRawMomentResult final_moments = select_moments(moments, request);
        std::cout << "5D direct Monte Carlo particle breakage and heat-transfer example\n";
        std::cout << "phase-space ordering: (V,T,u,v,w)\n";
        print_selected_moments("Initial tracked moments:", initial_moments);
        std::cout << '\n';
        print_selected_moments("Final tracked moments:", final_moments);
        if (output.value().write_moments)
        {
            std::cout << "\nOutputs written to:\n";
            std::cout << (output.value().directory / "moment_history_mc.txt").string() << '\n';
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
