#include "monte_carlo/evolution.hpp"
#include "qmomkit/configuration.hpp"
#include "thermal_particles.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>

int main()
{
    using namespace cqmom_thermal_spherical_breakup;
    try
    {
        const std::filesystem::path configuration_file = "configuration.qmom";
        const auto configured = qmomkit::config::load_cqmom(configuration_file);
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
        std::cout << "initial_sample_moment_max_abs_difference=" << maximum_tensor_difference(sampled->moments.moments, moments) << '\n';
        print_diagnostics("MC initial", diagnostics(moments));

        std::ofstream history;
        if (output.value().write_moments)
        {
            ensure_output_directory(output.value().directory);
            history.open(output.value().directory / "diagnostics_mc.txt");
            history << std::setprecision(17)
                    << "# time N volume volume_temperature volume_momentum_u volume_momentum_v volume_momentum_w volume_u2\n";
        }
        const auto write_history = [&history](const double time, const Diagnostics& values)
        {
            if (history)
            {
                history << time << ' ' << values.number_density << ' ' << values.volume << ' ' << values.volume_temperature << ' '
                        << values.momentum_u << ' ' << values.momentum_v << ' ' << values.momentum_w << ' ' << values.volume_u_square << '\n';
            }
        };
        write_history(initial_time, diagnostics(moments));

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
                write_history(time, diagnostics(moments));
                while (next_output <= time + 1.0e-14)
                {
                    next_output += output_interval;
                }
            }
        }
        print_diagnostics("MC final", diagnostics(moments));
        std::cout << "final_particle_count=" << ensemble.size() << '\n';
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
