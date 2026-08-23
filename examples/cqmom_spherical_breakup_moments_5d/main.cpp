#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

int main()
{
    using namespace cqmom_spherical_breakup_moments_5d;
    constexpr double initial_time = 0.0;
    constexpr double final_time = 10.0e-2;
    constexpr double time_step = 10.0e-4;
    try
    {
        const std::filesystem::path configuration_file = "configuration.qmom";
        const auto configured = qmomkit::config::load_cqmom(configuration_file);
        const auto evolution = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        const auto output = qmomkit::config::load_output_options(configuration_file);
        const auto text = read_configuration_text(configuration_file);
        if (!configured || !evolution || !output || !text)
        {
            std::cerr << "Could not load configuration.qmom.\n";
            return 1;
        }
        const auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(*text, configuration_file.string());
        if (!initial_condition)
        {
            std::cerr << initial_condition.error() << '\n';
            return 1;
        }
        const auto sampled = sample_initial_population(initial_condition.value(), configured.value().hierarchy());
        if (!sampled || !sampled->ok())
        {
            std::cerr << "Sampled initialization failed.\n";
            return 1;
        }
        const double particle_weight = sampled->sampling.number_density / static_cast<double>(sampled->sampling.ensemble.size());
        const qmomkit::RawMomentTensor mc_initial = ensemble_moments(sampled->sampling.ensemble, particle_weight, sampled->moments.moments.extents);
        if (maximum_tensor_difference(sampled->moments.moments, mc_initial) != 0.0)
        {
            throw std::runtime_error("The sampled CQMOM and Monte Carlo initial moments differ.");
        }

        const auto sampled_case = qmomkit::config::with_initial_moments(configured.value(), sampled->moments.moments);
        if (!sampled_case)
        {
            std::cerr << sampled_case.error() << '\n';
            return 1;
        }
        const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(sampled_case.value());
        if (!initial.ok())
        {
            std::cerr << "Initial CQMOM inversion failed: " << initial.message << '\n';
            return 1;
        }

        qmomkit::CqmomState state = sampled_case.value().initial_state();
        qmomkit::CqmomQuadrature quadrature = initial.quadrature;
        const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(sampled_case.value().hierarchy());
        const qmomkit::RawMomentRequest request{layout.required_indices};
        const qmomkit::SelectedRawMomentResult initial_moments = qmomkit::reconstruct_raw_moments(quadrature, request);

        const qmomkit::CqmomMomentSource source = breakup_source();
        std::vector<std::pair<double, qmomkit::CqmomQuadrature>> history{{initial_time, quadrature}};
        const double output_interval = output.value().interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output = initial_time + output_interval;
        while (time < final_time - 1.0e-14)
        {
            const qmomkit::CqmomPhaseSpaceEvolutionResult step = qmomkit::advance_cqmom_phase_space(
                state, time, std::min(time_step, final_time - time), particle_velocity(), source, evolution.value());
            if (!step.ok() || !step.quadrature.has_value())
            {
                std::cerr << "CQMOM phase-space evolution failed: " << step.message << '\n';
                return 1;
            }
            state = step.state;
            quadrature = *step.quadrature;
            time = step.final_time;
            if (time + 1.0e-14 >= next_output || time >= final_time - 1.0e-14)
            {
                history.emplace_back(time, quadrature);
                while (next_output <= time + 1.0e-14)
                {
                    next_output += output_interval;
                }
            }
        }
        const qmomkit::SelectedRawMomentResult final_moments = qmomkit::reconstruct_raw_moments(quadrature, request);

        if (output.value().write_moments || output.value().write_quadrature)
        {
            ensure_output_directory(output.value().directory);
        }
        if (output.value().write_moments)
        {
            qmomkit::MomentWriter writer = qmomkit::MomentWriter::for_cqmom(output.value().directory / "moment_history.txt", request);
            for (const auto& [sample_time, sample_quadrature] : history)
            {
                writer.write(sample_time, sample_quadrature);
            }
        }
        if (output.value().write_quadrature)
        {
            qmomkit::QuadratureWriter writer = qmomkit::QuadratureWriter::for_cqmom(output.value().directory / "quadrature_history.txt", 5);
            for (const auto& [sample_time, sample_quadrature] : history)
            {
                writer.write(sample_time, sample_quadrature);
            }
        }

        std::cout << "5D CQMOM thermal spherical-particle breakup example\n";
        std::cout << "phase-space ordering: (V,T,u,v,w)\n";
        print_selected_moments("Initial tracked moments:", initial_moments);
        std::cout << '\n';
        print_selected_moments("Final tracked moments:", final_moments);
        if (output.value().write_moments || output.value().write_quadrature)
        {
            std::cout << "\nOutputs written to:\n";
            if (output.value().write_moments)
            {
                std::cout << (output.value().directory / "moment_history.txt").string() << '\n';
            }
            if (output.value().write_quadrature)
            {
                std::cout << (output.value().directory / "quadrature_history.txt").string() << '\n';
            }
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
