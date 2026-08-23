#include "qmomkit/qmomkit.hpp"
#include "thermal_particles.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <iostream>
#include <utility>
#include <vector>

int main()
{
    using namespace cqmom_thermal_spherical_breakup;
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
        std::cout << "initial_sample_moment_max_abs_difference=" << maximum_tensor_difference(sampled->moments.moments, mc_initial) << '\n';

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
        const qmomkit::RawMomentRequest request = diagnostic_request();
        print_diagnostics("CQMOM initial", diagnostics(qmomkit::reconstruct_raw_moments(quadrature, request)));

        std::size_t generalized_requests = 0;
        const qmomkit::CqmomMomentSource source = breakup_source(generalized_requests);
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
        print_diagnostics("CQMOM final", diagnostics(qmomkit::reconstruct_raw_moments(quadrature, request)));
        std::cout << "fractional_generalized_moment_requests=" << generalized_requests << '\n';

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
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
