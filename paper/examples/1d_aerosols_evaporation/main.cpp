#include "monte_carlo/evolution.hpp"
#include "monte_carlo/moments.hpp"
#include "monte_carlo/sampling.hpp"
#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace
{

// Writes the number density M0, the Sauter mean diameter d32 = M3/M2 and the
// diameter variance M2/M0 - (M1/M0)^2 of the raw moments M0 to M3.
void write_statistics(std::ostream& output, const double time, const std::vector<double>& moments)
{
    const double mean_diameter = moments[1] / moments[0];
    output << time << ' ' << moments[0] << ' ' << moments[3] / moments[2] << ' ' << moments[2] / moments[0] - mean_diameter * mean_diameter << '\n';
}

} // namespace

int main()
{
    constexpr std::array<std::size_t, 5> fixed_node_counts{{2, 4, 6, 8, 10}};
    constexpr std::size_t sample_count = 10000;
    constexpr unsigned long long seed = 12345;
    constexpr double initial_time = 0.0;
    constexpr double final_time = 20.0;
    constexpr double time_step = 0.2;
    const std::filesystem::path configuration_file = "configuration.qmom";

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
    const qmomkit::config::ConfigurationResult<qmomkit::OutputOptions> configured_output = qmomkit::config::load_output_options(configuration_file);
    if (!configured_output)
    {
        std::cerr << configured_output.error() << '\n';
        return 1;
    }
    const qmomkit::NodeAdvectionOptions& node_advection_options = configured_node_advection.value();
    const qmomkit::OutputOptions& output = configured_output.value();

    // The Monte Carlo companion samples the configured mixture instead of
    // computing its moments exactly.
    std::ifstream configuration_input(configuration_file);
    if (!configuration_input)
    {
        std::cerr << "Could not open configuration file: " << configuration_file << '\n';
        return 1;
    }
    std::ostringstream configuration_text;
    configuration_text << configuration_input.rdbuf();
    const qmomkit::config::ConfigurationResult<qmomkit::config::InitialConditionSpecification> initial_configuration =
        qmomkit::config::parse_qmom_initial_condition(configuration_text.str(), configuration_file.string());
    if (!initial_configuration)
    {
        std::cerr << initial_configuration.error() << '\n';
        return 1;
    }
    const auto* mixture_initialization = std::get_if<qmomkit::config::MixtureInitialization>(&initial_configuration.value());
    if (mixture_initialization == nullptr)
    {
        std::cerr << "The aerosol evaporation study requires a mixture initialization.\n";
        return 1;
    }
    const qmomkit::monte_carlo::SamplingResult generated_samples =
        qmomkit::monte_carlo::sample_mixture(mixture_initialization->mixture, qmomkit::monte_carlo::SamplingOptions{sample_count, seed});
    if (!generated_samples.ok())
    {
        std::cerr << "Monte Carlo sampling failed: " << generated_samples.message << '\n';
        return 1;
    }

    // The study truncates the sampled mixture to the modeled diameter support and
    // assigns the configured total number density to the retained aerosols. The
    // QMOM runs start from the moments of this truncated sample.
    const double number_density = generated_samples.number_density;
    qmomkit::monte_carlo::Ensemble ensemble;
    std::copy_if(generated_samples.ensemble.begin(), generated_samples.ensemble.end(), std::back_inserter(ensemble),
                 [](const qmomkit::monte_carlo::Sample& sample) { return paper_1d_aerosols_evaporation::is_within_modeled_support(sample[0]); });
    const qmomkit::monte_carlo::OneDimensionalMomentEstimate initial_moment_estimate =
        qmomkit::monte_carlo::compute_ensemble_moments(ensemble, number_density, mixture_initialization->maximum_order);
    if (!initial_moment_estimate.ok())
    {
        std::cerr << "Monte Carlo initial moment estimation failed: " << initial_moment_estimate.message << '\n';
        return 1;
    }

    const qmomkit::PhaseSpaceVelocity velocity = paper_1d_aerosols_evaporation::evaporation_velocity();

    try
    {
        std::filesystem::create_directories(output.directory);
        const auto open_output = [&](const std::string& file_name)
        {
            std::ofstream stream(output.directory / file_name);
            if (!stream)
            {
                throw std::runtime_error("Could not open " + (output.directory / file_name).string());
            }
            return stream;
        };

        std::ofstream initial_sample_output = open_output("initial_mc_samples.txt");
        initial_sample_output << "# generated_sample_count=" << sample_count << '\n';
        initial_sample_output << "# retained_sample_count=" << ensemble.size() << '\n';
        initial_sample_output << "# retained_number_density=" << number_density << '\n';
        initial_sample_output << "# diameter_unit=um\n";
        initial_sample_output << "# diameter\n";
        initial_sample_output << std::setprecision(17);
        for (const qmomkit::monte_carlo::Sample& sample : ensemble)
        {
            initial_sample_output << sample[0] << '\n';
        }

        std::ofstream mc_output = open_output("moment_history_mc.txt");
        mc_output << "# generated_sample_count=" << sample_count << '\n';
        mc_output << "# retained_sample_count=" << ensemble.size() << '\n';
        mc_output << "# retained_number_density=" << number_density << '\n';
        mc_output << "# seed=" << seed << '\n';
        mc_output << "# M0_unit=m^-3 d32_unit=um diameter_variance_unit=um^2\n";
        mc_output << "# time M0 d32 diameter_variance\n";
        mc_output << std::setprecision(17);
        write_statistics(mc_output, initial_time, initial_moment_estimate.moments.values());

        // Each fixed node count is a separate QMOM run from the same sampled moments.
        std::array<std::ofstream, fixed_node_counts.size()> qmom_outputs;
        std::array<qmomkit::Quadrature1D, fixed_node_counts.size()> qmom_quadratures;
        for (std::size_t case_index = 0; case_index < fixed_node_counts.size(); ++case_index)
        {
            const std::size_t node_count = fixed_node_counts[case_index];
            const qmomkit::QmomResult initial_reconstruction =
                qmomkit::compute_qmom(qmomkit::QmomCase{initial_moment_estimate.moments, node_count, configured_qmom.value().qmom_options});
            if (!initial_reconstruction.ok())
            {
                std::cerr << "Initial QMOM reconstruction failed for N=" << node_count << ": " << initial_reconstruction.message << '\n';
                return 1;
            }
            qmom_quadratures[case_index].abscissas = initial_reconstruction.abscissas;
            qmom_quadratures[case_index].weights = initial_reconstruction.weights;

            qmom_outputs[case_index] = open_output("moment_history_qmom_n" + std::to_string(node_count) + ".txt");
            qmom_outputs[case_index] << "# fixed_qmom_node_count=" << node_count << '\n';
            qmom_outputs[case_index] << "# M0_unit=m^-3 d32_unit=um diameter_variance_unit=um^2\n";
            qmom_outputs[case_index] << "# time M0 d32 diameter_variance\n";
            qmom_outputs[case_index] << std::setprecision(17);
            write_statistics(qmom_outputs[case_index], initial_time, qmomkit::reconstruct_raw_moments(qmom_quadratures[case_index], 4));
        }

        double time = initial_time;
        while (time < final_time - 1.0e-12)
        {
            const double step = std::min(time_step, final_time - time);

            // Evaporated aerosols keep their residual diameter. The floor changes
            // only abscissas, so the reported moments are reconstructed from the
            // floored nodes and their unchanged weights.
            for (std::size_t case_index = 0; case_index < fixed_node_counts.size(); ++case_index)
            {
                const qmomkit::QmomNodeAdvectionResult qmom_advection =
                    qmomkit::advance_nodes(qmom_quadratures[case_index], time, step, velocity, node_advection_options);
                if (!qmom_advection.ok())
                {
                    std::cerr << "QMOM advection failed for N=" << fixed_node_counts[case_index] << ": " << qmom_advection.message << '\n';
                    return 1;
                }
                qmom_quadratures[case_index] = qmom_advection.quadrature;
                paper_1d_aerosols_evaporation::enforce_residual_diameter(qmom_quadratures[case_index]);
            }

            const qmomkit::monte_carlo::EnsembleEvolutionResult mc_advection =
                qmomkit::monte_carlo::evolve_ensemble(ensemble, velocity, time, time + step, node_advection_options.adaptive_options);
            if (!mc_advection.ok())
            {
                std::cerr << "Monte Carlo advection failed: " << mc_advection.message << '\n';
                return 1;
            }
            ensemble = mc_advection.ensemble;
            paper_1d_aerosols_evaporation::enforce_residual_diameter(ensemble);
            const qmomkit::monte_carlo::OneDimensionalMomentEstimate mc_moment_estimate =
                qmomkit::monte_carlo::compute_ensemble_moments(ensemble, number_density, 3);
            if (!mc_moment_estimate.ok())
            {
                std::cerr << "Monte Carlo moment estimation failed: " << mc_moment_estimate.message << '\n';
                return 1;
            }

            time += step;
            for (std::size_t case_index = 0; case_index < fixed_node_counts.size(); ++case_index)
            {
                write_statistics(qmom_outputs[case_index], time, qmomkit::reconstruct_raw_moments(qmom_quadratures[case_index], 4));
            }
            write_statistics(mc_output, time, mc_moment_estimate.moments.values());
        }

        std::cout << "Output written to " << output.directory.string() << "/\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
