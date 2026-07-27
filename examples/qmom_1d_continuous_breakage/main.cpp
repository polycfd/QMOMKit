#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <vector>

namespace
{
constexpr double initial_time = 0.0;
constexpr double final_time = 1.0;
constexpr double time_step = 0.01;

void print_moments(const char* label, const std::vector<double>& moments)
{
    std::cout << label << '\n';
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        std::cout << "M" << order << " = " << moments[order] << '\n';
    }
}

void print_quadrature(const char* label, const qmomkit::Quadrature1D& quadrature)
{
    std::cout << label << '\n';
    std::cout << "node  weight        abscissa\n";
    for (std::size_t node = 0; node < quadrature.weights.size(); ++node)
    {
        std::cout << node << "     " << quadrature.weights[node] << "   " << quadrature.abscissas[node] << '\n';
    }
}
} // namespace

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";

    const auto configured = qmomkit::config::load_qmom(configuration_file);
    if (!configured.ok())
    {
        qmom_1d_continuous_breakage::print_configuration_errors(configured.errors);
        return 1;
    }

    const auto output_configured = qmomkit::config::load_output_options(configuration_file);
    if (!output_configured.ok())
    {
        qmom_1d_continuous_breakage::print_configuration_errors(output_configured.errors);
        return 1;
    }

    const qmomkit::QmomCase qmom_case = *configured.value;
    const qmomkit::OutputOptions output_options = *output_configured.value;
    const std::vector<double> initial_moments = qmom_case.moments.values();
    std::vector<double> moments = initial_moments;
    const qmomkit::QmomResult result = qmomkit::compute_qmom(qmom_case);
    if (!result.ok())
    {
        std::cerr << "QMOM inversion failed: " << result.message << '\n';
        return 1;
    }

    qmomkit::Quadrature1D quadrature{result.abscissas, result.weights};
    const qmomkit::PhaseSpaceVelocity velocity = qmom_1d_continuous_breakage::continuous_velocity();
    const qmomkit::MomentSource source = qmom_1d_continuous_breakage::binary_fragmentation_source();
    qmomkit::NodeAdvectionOptions node_options;
    qmomkit::MomentSourceUpdateOptions source_options;

    std::optional<qmomkit::MomentWriter> moment_writer;
    std::optional<qmomkit::QuadratureWriter> quadrature_writer;
    try
    {
        if (output_options.write_moments || output_options.write_quadrature)
        {
            std::filesystem::create_directories(output_options.directory);
        }
        if (output_options.write_moments)
        {
            moment_writer.emplace(qmomkit::MomentWriter::for_qmom(output_options.directory / "moments.txt", moments.size()));
            moment_writer->write(initial_time, qmomkit::MomentSet(moments));
        }
        if (output_options.write_quadrature)
        {
            quadrature_writer.emplace(qmomkit::QuadratureWriter::for_qmom(output_options.directory / "quadrature.txt", quadrature.weights.size()));
            quadrature_writer->write(initial_time, quadrature);
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    double time = initial_time;
    const double output_interval = output_options.interval.value_or(final_time - initial_time);
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
                if (moment_writer.has_value())
                {
                    moment_writer->write(time, qmomkit::MomentSet(moments));
                }
                if (quadrature_writer.has_value())
                {
                    quadrature_writer->write(time, quadrature);
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
    std::cout << "1D QMOM continuous + discontinuous evolution example\n";

    std::cout << '\n';
    print_moments("Initial moments:", initial_moments);

    std::cout << '\n';
    print_moments("Final moments:", moments);

    std::cout << '\n';
    print_quadrature("Initial quadrature:", qmomkit::Quadrature1D{result.abscissas, result.weights});

    std::cout << '\n';
    print_quadrature("Final quadrature:", quadrature);
    std::cout << "final_time = " << time << '\n';

    if (output_options.write_moments || output_options.write_quadrature)
    {
        std::cout << "\nOutputs written to:\n";
        if (output_options.write_moments)
        {
            std::cout << (output_options.directory / "moments.txt").string() << '\n';
        }
        if (output_options.write_quadrature)
        {
            std::cout << (output_options.directory / "quadrature.txt").string() << '\n';
        }
    }
    else
    {
        std::cout << "\nOutput disabled by configuration.\n";
    }

    return 0;
}
