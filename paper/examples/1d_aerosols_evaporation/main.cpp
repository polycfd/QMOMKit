#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "monte_carlo/sampling.hpp"
#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

int main()
{
    constexpr std::array<std::size_t, 5> fixed_node_counts{{2, 4, 6, 8, 10}};
    constexpr std::size_t initial_moment_order = 19;
    constexpr std::size_t sample_count = 10000;
    constexpr unsigned long long seed = 12345;
    constexpr double retained_number_density = 2.23e5;
    constexpr double initial_time = 0.0;
    constexpr double final_time = 20.0;
    constexpr double time_step = 0.2;
    const std::filesystem::path configuration_file = "configuration.qmom";

    try
    {
        const qmomkit::config::ConfigurationResult<qmomkit::QmomCase> configured_qmom = qmomkit::config::load_qmom(configuration_file);
        if (!configured_qmom)
        {
            std::cerr << configured_qmom.error() << '\n';
            return 1;
        }

        const qmomkit::config::ConfigurationResult<qmomkit::NodeAdvectionOptions> configured_node_advection =
            qmomkit::config::load_node_advection_options(configuration_file);
        if (!configured_node_advection)
        {
            std::cerr << configured_node_advection.error() << '\n';
            return 1;
        }
        const qmomkit::NodeAdvectionOptions node_advection_options = configured_node_advection.value();
        if (node_advection_options.scheme != qmomkit::NodeAdvectionScheme::rk23)
        {
            std::cerr << "The aerosol evaporation example requires RK23 node advection.\n";
            return 1;
        }

        const qmomkit::config::ConfigurationResult<qmomkit::OutputOptions> configured_output =
            qmomkit::config::load_output_options(configuration_file);
        if (!configured_output)
        {
            std::cerr << configured_output.error() << '\n';
            return 1;
        }

        std::ifstream configuration_input(configuration_file);
        if (!configuration_input)
        {
            std::cerr << "Could not open " << configuration_file << ".\n";
            return 1;
        }
        std::ostringstream configuration_buffer;
        configuration_buffer << configuration_input.rdbuf();
        const qmomkit::config::ConfigurationResult<qmomkit::config::InitialConditionSpecification> initial_configuration =
            qmomkit::config::parse_qmom_initial_condition(configuration_buffer.str(), configuration_file.string());
        if (!initial_configuration)
        {
            std::cerr << initial_configuration.error() << '\n';
            return 1;
        }
        const qmomkit::config::MixtureInitialization* mixture_initialization =
            std::get_if<qmomkit::config::MixtureInitialization>(&initial_configuration.value());
        if (mixture_initialization == nullptr)
        {
            std::cerr << "The aerosol evaporation configuration must define a mixture initialization.\n";
            return 1;
        }

        const qmomkit::monte_carlo::SamplingResult generated_samples =
            qmomkit::monte_carlo::sample_mixture(mixture_initialization->mixture, qmomkit::monte_carlo::SamplingOptions{sample_count, seed});
        if (!generated_samples.ok())
        {
            std::cerr << "Monte Carlo sampling failed: " << generated_samples.message << '\n';
            return 1;
        }
        if (generated_samples.ensemble.size() != sample_count)
        {
            std::cerr << "Monte Carlo sampling did not generate the configured sample count.\n";
            return 1;
        }

        qmomkit::monte_carlo::Ensemble ensemble;
        ensemble.reserve(generated_samples.ensemble.size());
        for (const qmomkit::monte_carlo::Sample& sample : generated_samples.ensemble)
        {
            if (sample.size() != 1)
            {
                std::cerr << "Monte Carlo sampling did not produce one-dimensional samples.\n";
                return 1;
            }
            if (paper_1d_aerosols_evaporation::is_within_modeled_support(sample[0]))
            {
                ensemble.push_back(sample);
            }
        }
        if (ensemble.empty())
        {
            std::cerr << "Support filtering removed every sampled aerosol diameter.\n";
            return 1;
        }

        const qmomkit::monte_carlo::OneDimensionalMomentEstimate initial_moment_estimate =
            qmomkit::monte_carlo::compute_ensemble_moments(ensemble, retained_number_density, initial_moment_order);
        if (!initial_moment_estimate.ok())
        {
            std::cerr << "Monte Carlo initial moment evaluation failed: " << initial_moment_estimate.message << '\n';
            return 1;
        }

        const qmomkit::PhaseSpaceVelocity velocity = paper_1d_aerosols_evaporation::evaporation_velocity();
        const qmomkit::AdaptiveRkOptions mc_advection_options = node_advection_options.adaptive_options;
        const qmomkit::OutputOptions output_options = configured_output.value();
        std::filesystem::create_directories(output_options.directory);

        std::array<std::ofstream, fixed_node_counts.size()> qmom_outputs;
        qmom_outputs[0].open(output_options.directory / "moment_history_qmom_n2.txt");
        qmom_outputs[1].open(output_options.directory / "moment_history_qmom_n4.txt");
        qmom_outputs[2].open(output_options.directory / "moment_history_qmom_n6.txt");
        qmom_outputs[3].open(output_options.directory / "moment_history_qmom_n8.txt");
        qmom_outputs[4].open(output_options.directory / "moment_history_qmom_n10.txt");
        std::ofstream mc_output(output_options.directory / "moment_history_mc.txt");
        std::ofstream initial_sample_output(output_options.directory / "initial_mc_samples.txt");
        std::ofstream initial_quadrature_output(output_options.directory / "initial_quadrature.txt");
        if (!mc_output || !initial_sample_output || !initial_quadrature_output)
        {
            throw std::runtime_error("Could not open aerosol evaporation output files.");
        }
        for (std::ofstream& qmom_output : qmom_outputs)
        {
            if (!qmom_output)
            {
                throw std::runtime_error("Could not open QMOM aerosol evaporation output.");
            }
        }

        for (std::size_t case_index = 0; case_index < fixed_node_counts.size(); ++case_index)
        {
            qmom_outputs[case_index] << "# fixed_qmom_node_count=" << fixed_node_counts[case_index] << '\n';
            qmom_outputs[case_index] << "# M0_unit=m^-3 d32_unit=um diameter_variance_unit=um^2\n";
            qmom_outputs[case_index] << "# time M0 d32 diameter_variance\n";
            qmom_outputs[case_index] << std::setprecision(17);
        }
        mc_output << "# generated_sample_count=" << sample_count << '\n';
        mc_output << "# retained_sample_count=" << ensemble.size() << '\n';
        mc_output << "# retained_number_density=" << retained_number_density << '\n';
        mc_output << "# seed=" << seed << '\n';
        mc_output << "# M0_unit=m^-3 d32_unit=um diameter_variance_unit=um^2\n";
        mc_output << "# time M0 d32 diameter_variance\n";
        initial_sample_output << "# generated_sample_count=" << sample_count << '\n';
        initial_sample_output << "# retained_sample_count=" << ensemble.size() << '\n';
        initial_sample_output << "# retained_number_density=" << retained_number_density << '\n';
        initial_sample_output << "# diameter_unit=um\n";
        initial_sample_output << "# diameter\n";
        initial_quadrature_output << "# fixed_qmom_node_count=10\n";
        initial_quadrature_output << "# diameter_unit=um\n";
        initial_quadrature_output << "# node_index diameter weight\n";
        mc_output << std::setprecision(17);
        initial_sample_output << std::setprecision(17);
        initial_quadrature_output << std::setprecision(17);

        for (const qmomkit::monte_carlo::Sample& sample : ensemble)
        {
            initial_sample_output << sample[0] << '\n';
        }

        std::array<qmomkit::Quadrature1D, fixed_node_counts.size()> qmom_quadratures;
        for (std::size_t case_index = 0; case_index < fixed_node_counts.size(); ++case_index)
        {
            const std::size_t node_count = fixed_node_counts[case_index];
            const qmomkit::QmomCase qmom_case{initial_moment_estimate.moments, node_count, configured_qmom.value().options};
            const qmomkit::QmomResult initial_reconstruction = qmomkit::compute_qmom(qmom_case);
            if (!initial_reconstruction.ok())
            {
                std::cerr << "Initial QMOM reconstruction failed for N=" << node_count << ": " << initial_reconstruction.message << '\n';
                return 1;
            }
            if (initial_reconstruction.weights.size() != node_count)
            {
                std::cerr << "Initial QMOM reconstruction did not retain all fixed nodes for N=" << node_count << ".\n";
                return 1;
            }

            qmom_quadratures[case_index] = {initial_reconstruction.abscissas, initial_reconstruction.weights, qmomkit::QuadratureStatus::success, ""};
            for (const double diameter : qmom_quadratures[case_index].abscissas)
            {
                if (!paper_1d_aerosols_evaporation::is_within_modeled_support(diameter))
                {
                    std::cerr << "Initial QMOM node is outside the modeled support for N=" << node_count << ".\n";
                    return 1;
                }
            }

            const std::vector<double> qmom_moments = qmomkit::reconstruct_raw_moments(qmom_quadratures[case_index], 2 * node_count);
            if (!std::isfinite(qmom_moments[2]) || qmom_moments[2] <= 0.0 || !std::isfinite(qmom_moments[0]) || qmom_moments[0] <= 0.0)
            {
                throw std::runtime_error("Initial QMOM moments are not suitable for reported diameter statistics.");
            }
            const double qmom_d32 = qmom_moments[3] / qmom_moments[2];
            const double qmom_variance =
                qmom_moments[2] / qmom_moments[0] - (qmom_moments[1] / qmom_moments[0]) * (qmom_moments[1] / qmom_moments[0]);
            if (!std::isfinite(qmom_d32) || !std::isfinite(qmom_variance))
            {
                throw std::runtime_error("Initial QMOM diameter statistics are non-finite.");
            }
            qmom_outputs[case_index] << initial_time << ' ' << qmom_moments[0] << ' ' << qmom_d32 << ' ' << qmom_variance << '\n';

            if (node_count == 10)
            {
                for (std::size_t node = 0; node < qmom_quadratures[case_index].weights.size(); ++node)
                {
                    initial_quadrature_output << node << ' ' << qmom_quadratures[case_index].abscissas[node] << ' '
                                              << qmom_quadratures[case_index].weights[node] << '\n';
                }
            }
        }

        if (!std::isfinite(initial_moment_estimate.moments[2]) || initial_moment_estimate.moments[2] <= 0.0)
        {
            throw std::runtime_error("Initial Monte Carlo moments are not suitable for reported diameter statistics.");
        }
        const double initial_mc_d32 = initial_moment_estimate.moments[3] / initial_moment_estimate.moments[2];
        const double initial_mc_variance = initial_moment_estimate.moments[2] / initial_moment_estimate.moments[0] -
                                           (initial_moment_estimate.moments[1] / initial_moment_estimate.moments[0]) *
                                               (initial_moment_estimate.moments[1] / initial_moment_estimate.moments[0]);
        if (!std::isfinite(initial_mc_d32) || !std::isfinite(initial_mc_variance))
        {
            throw std::runtime_error("Initial Monte Carlo diameter statistics are non-finite.");
        }
        mc_output << initial_time << ' ' << initial_moment_estimate.moments[0] << ' ' << initial_mc_d32 << ' ' << initial_mc_variance << '\n';

        double time = initial_time;
        while (time < final_time)
        {
            const double step = std::min(time_step, final_time - time);

            for (std::size_t case_index = 0; case_index < fixed_node_counts.size(); ++case_index)
            {
                const std::size_t node_count = fixed_node_counts[case_index];
                const qmomkit::QmomNodeAdvectionResult qmom_advection =
                    qmomkit::advance_nodes(qmom_quadratures[case_index], time, step, velocity, node_advection_options);
                if (!qmom_advection.ok())
                {
                    std::cerr << "QMOM advection failed for N=" << node_count << ": " << qmom_advection.message << '\n';
                    return 1;
                }
                qmom_quadratures[case_index] = qmom_advection.quadrature;
                paper_1d_aerosols_evaporation::enforce_residual_diameter(qmom_quadratures[case_index]);
                if (qmom_quadratures[case_index].weights.size() != node_count)
                {
                    std::cerr << "QMOM node count changed during fixed-node advection for N=" << node_count << ".\n";
                    return 1;
                }
                for (const double diameter : qmom_quadratures[case_index].abscissas)
                {
                    if (!paper_1d_aerosols_evaporation::is_within_modeled_support(diameter))
                    {
                        std::cerr << "QMOM node left the modeled support for N=" << node_count << ".\n";
                        return 1;
                    }
                }

                // The floor operation changes only abscissas; this
                // reconstruction keeps the reported raw moments consistent
                // with the fixed-node coordinates and unchanged weights.
                const std::vector<double> qmom_moments = qmomkit::reconstruct_raw_moments(qmom_quadratures[case_index], 2 * node_count);
                if (!std::isfinite(qmom_moments[2]) || qmom_moments[2] <= 0.0 || !std::isfinite(qmom_moments[0]) || qmom_moments[0] <= 0.0)
                {
                    throw std::runtime_error("QMOM moments are not suitable for reported diameter statistics.");
                }
                const double qmom_d32 = qmom_moments[3] / qmom_moments[2];
                const double qmom_variance =
                    qmom_moments[2] / qmom_moments[0] - (qmom_moments[1] / qmom_moments[0]) * (qmom_moments[1] / qmom_moments[0]);
                if (!std::isfinite(qmom_d32) || !std::isfinite(qmom_variance))
                {
                    throw std::runtime_error("QMOM diameter statistics are non-finite.");
                }
                qmom_outputs[case_index] << time + step << ' ' << qmom_moments[0] << ' ' << qmom_d32 << ' ' << qmom_variance << '\n';
            }

            const qmomkit::monte_carlo::EnsembleEvolutionResult mc_advection =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + step, mc_advection_options);
            if (!mc_advection.ok())
            {
                std::cerr << "Monte Carlo advection failed: " << mc_advection.message << '\n';
                return 1;
            }
            ensemble = mc_advection.ensemble;
            paper_1d_aerosols_evaporation::enforce_residual_diameter(ensemble);
            for (const qmomkit::monte_carlo::Sample& sample : ensemble)
            {
                if (!paper_1d_aerosols_evaporation::is_within_modeled_support(sample[0]))
                {
                    std::cerr << "Monte Carlo sample left the modeled support.\n";
                    return 1;
                }
            }

            const qmomkit::monte_carlo::OneDimensionalMomentEstimate mc_moment_estimate =
                qmomkit::monte_carlo::compute_ensemble_moments(ensemble, retained_number_density, 3);
            if (!mc_moment_estimate.ok() || !std::isfinite(mc_moment_estimate.moments[2]) || mc_moment_estimate.moments[2] <= 0.0)
            {
                throw std::runtime_error("Monte Carlo moments are not suitable for reported diameter statistics.");
            }
            const double mc_d32 = mc_moment_estimate.moments[3] / mc_moment_estimate.moments[2];
            const double mc_variance =
                mc_moment_estimate.moments[2] / mc_moment_estimate.moments[0] -
                (mc_moment_estimate.moments[1] / mc_moment_estimate.moments[0]) * (mc_moment_estimate.moments[1] / mc_moment_estimate.moments[0]);
            if (!std::isfinite(mc_d32) || !std::isfinite(mc_variance))
            {
                throw std::runtime_error("Monte Carlo diameter statistics are non-finite.");
            }
            time += step;
            mc_output << time << ' ' << mc_moment_estimate.moments[0] << ' ' << mc_d32 << ' ' << mc_variance << '\n';
        }

        std::cout << "generated_sample_count = " << sample_count << '\n';
        std::cout << "retained_sample_count = " << ensemble.size() << '\n';
        std::cout << "retained_number_density = " << retained_number_density << '\n';
        std::cout << qmomkit::format_node_advection_options(node_advection_options);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
