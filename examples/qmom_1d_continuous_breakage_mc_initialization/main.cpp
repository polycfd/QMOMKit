#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <utility>
#include <vector>

namespace
{
constexpr double initial_time = 0.0;
constexpr double final_time = 1.0;
constexpr double time_step = 0.01;
} // namespace

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";

    const std::optional<std::string> configuration_text = qmom_1d_continuous_breakage_mc_initialization::read_configuration_text(configuration_file);
    if (!configuration_text.has_value())
    {
        return 1;
    }

    const auto configured = qmomkit::config::load_qmom(configuration_file);
    if (!configured.ok())
    {
        qmom_1d_continuous_breakage_mc_initialization::print_configuration_errors(configured.errors);
        return 1;
    }

    const auto output_configured = qmomkit::config::load_output_options(configuration_file);
    if (!output_configured.ok())
    {
        qmom_1d_continuous_breakage_mc_initialization::print_configuration_errors(output_configured.errors);
        return 1;
    }

    const qmomkit::QmomCase qmom_case = *configured.value;
    const qmomkit::OutputOptions output_options = *output_configured.value;
    const auto initial_condition = qmomkit::config::parse_qmom_initial_condition(*configuration_text, configuration_file.string());
    if (!initial_condition.ok())
    {
        qmom_1d_continuous_breakage_mc_initialization::print_configuration_errors(initial_condition.errors);
        return 1;
    }

    const std::optional<qmom_1d_continuous_breakage_mc_initialization::SamplingSetup> sampled =
        qmom_1d_continuous_breakage_mc_initialization::sample_initial_distribution(*initial_condition.value);
    if (!sampled.has_value())
    {
        return 1;
    }
    if (!sampled->sampling.ok())
    {
        std::cerr << sampled->sampling.message << '\n';
        return 1;
    }

    const qmomkit::monte_carlo::OneDimensionalMomentEstimate sampled_moments =
        qmom_1d_continuous_breakage_mc_initialization::compute_sampled_moments(*sampled);
    if (!sampled_moments.ok())
    {
        std::cerr << sampled_moments.message << '\n';
        return 1;
    }

    const std::vector<double> initial_moments = sampled_moments.moments.values();
    std::vector<double> moments = initial_moments;
    const qmomkit::QmomCase sampled_case{sampled_moments.moments, qmom_case.node_count, qmom_case.options};
    const qmomkit::QmomResult result = qmomkit::compute_qmom(sampled_case);
    if (!result.ok())
    {
        std::cerr << "QMOM inversion failed: " << result.message << '\n';
        return 1;
    }

    qmomkit::Quadrature1D quadrature{result.abscissas, result.weights};
    const std::vector<double> initial_abscissas = quadrature.abscissas;
    const qmomkit::PhaseSpaceVelocity velocity = qmom_1d_continuous_breakage_mc_initialization::continuous_velocity();
    const qmomkit::MomentSource source = qmom_1d_continuous_breakage_mc_initialization::binary_fragmentation_source();
    qmomkit::NodeAdvectionOptions node_options;
    qmomkit::MomentSourceUpdateOptions source_options;
    std::optional<qmomkit::MomentWriter> moment_history_writer;
    std::optional<qmomkit::QuadratureWriter> quadrature_history_writer;
    try
    {
        if (output_options.write_moments || output_options.write_quadrature)
        {
            std::filesystem::create_directories(output_options.directory);
        }
        if (output_options.write_moments)
        {
            qmomkit::MomentWriter initial_moment_writer =
                qmomkit::MomentWriter::for_qmom(output_options.directory / "initial_moments.txt", initial_moments.size());
            initial_moment_writer.write(initial_time, qmomkit::MomentSet(initial_moments));

            moment_history_writer.emplace(qmomkit::MomentWriter::for_qmom(output_options.directory / "moment_history.txt", initial_moments.size()));
            moment_history_writer->write(initial_time, qmomkit::MomentSet(initial_moments));
        }
        if (output_options.write_quadrature)
        {
            qmomkit::QuadratureWriter initial_quadrature_writer =
                qmomkit::QuadratureWriter::for_qmom(output_options.directory / "initial_quadrature.txt", quadrature.weights.size());
            initial_quadrature_writer.write(initial_time, quadrature);

            quadrature_history_writer.emplace(
                qmomkit::QuadratureWriter::for_qmom(output_options.directory / "quadrature_history.txt", quadrature.weights.size()));
            quadrature_history_writer->write(initial_time, quadrature);
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

        // Strang stage 1: continuous half-step acts on quadrature nodes.
        const qmomkit::QmomNodeAdvectionResult first_continuous = qmomkit::advance_nodes(quadrature, time, half_step, velocity, node_options);
        if (!first_continuous.ok())
        {
            std::cerr << "First continuous half-step failed: " << first_continuous.message << '\n';
            return 1;
        }

        quadrature = first_continuous.quadrature;
        moments = qmomkit::reconstruct_raw_moments(quadrature, moments.size());

        // Strang stage 2: discontinuous fragmentation source acts on moments.
        const qmomkit::MomentSourceUpdateResult source_update =
            qmomkit::advance_moments_with_source(moments, time + half_step, current_step, source, source_options);
        if (!source_update.ok())
        {
            std::cerr << "Discontinuous moment source step failed: " << source_update.message << '\n';
            return 1;
        }

        moments = source_update.moments;

        // Reconstruct QMOM after the source update before advecting nodes again.
        const qmomkit::QmomResult post_source = qmomkit::compute_qmom(qmomkit::MomentSet(moments), qmom_case.node_count, qmom_case.options);
        if (!post_source.ok())
        {
            std::cerr << "Post-source QMOM inversion failed: " << post_source.message << '\n';
            return 1;
        }

        quadrature = qmomkit::Quadrature1D{post_source.abscissas, post_source.weights};

        // Strang stage 3: continuous half-step acts on the reconstructed nodes.
        const qmomkit::QmomNodeAdvectionResult second_continuous =
            qmomkit::advance_nodes(quadrature, time + half_step, half_step, velocity, node_options);
        if (!second_continuous.ok())
        {
            std::cerr << "Second continuous half-step failed: " << second_continuous.message << '\n';
            return 1;
        }

        quadrature = second_continuous.quadrature;
        moments = qmomkit::reconstruct_raw_moments(quadrature, moments.size());
        time += current_step;

        if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
        {
            try
            {
                if (moment_history_writer.has_value())
                {
                    moment_history_writer->write(time, qmomkit::MomentSet(moments));
                }
                if (quadrature_history_writer.has_value())
                {
                    quadrature_history_writer->write(time, quadrature);
                }
            }
            catch (const std::exception& error)
            {
                std::cerr << error.what() << '\n';
                return 1;
            }
            while (next_output_time <= time + 1.0e-12)
            {
                next_output_time += output_interval;
            }
        }
    }

    std::cout << std::setprecision(12);

    qmom_1d_continuous_breakage_mc_initialization::print_common_run_summary(
        std::cout, "QMOM quadrature", configuration_file, output_options.directory,
        qmom_1d_continuous_breakage_mc_initialization::initialization_sample_count,
        qmom_1d_continuous_breakage_mc_initialization::initialization_seed, initial_time, final_time, time_step, output_interval);
    std::cout << "Representation: sampled moments reconstructed as QMOM nodes and weights\n";

    qmom_1d_continuous_breakage_mc_initialization::print_moments(std::cout, "Initial moments", initial_moments);

    std::cout << "\nInitial quadrature:\n";
    std::cout << "node  weight        abscissa\n";
    for (std::size_t node = 0; node < result.weights.size(); ++node)
    {
        std::cout << node << "     " << result.weights[node] << "   " << result.abscissas[node] << '\n';
    }

    qmom_1d_continuous_breakage_mc_initialization::print_moments(std::cout, "Final moments", moments);

    std::cout << "\nFinal quadrature at t = " << time << ":\n";
    std::cout << "node  initial_abscissa  final_abscissa\n";
    for (std::size_t node = 0; node < quadrature.abscissas.size(); ++node)
    {
        std::cout << node << "     " << initial_abscissas[node] << "         " << quadrature.abscissas[node] << '\n';
    }
    if (output_options.write_moments || output_options.write_quadrature)
    {
        std::cout << "\nCompleted QMOM workflow. Outputs written to " << output_options.directory << '\n';
    }
    else
    {
        std::cout << "\nCompleted QMOM workflow. Output disabled by configuration.\n";
    }

    return 0;
}
