#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const int step_count = 800;
    const double time_step = 0.01;

    // The case and the output options come from the same file. A faulty file is
    // reported as a message naming the key and line.
    auto configured_case = qmomkit::config::load_qmom(configuration_file);
    if (!configured_case)
    {
        std::cerr << configured_case.error() << '\n';
        return 1;
    }
    auto configured_output = qmomkit::config::load_output_options(configuration_file);
    if (!configured_output)
    {
        std::cerr << configured_output.error() << '\n';
        return 1;
    }
    qmomkit::QmomCase gqmom_case = std::move(configured_case).value();
    const qmomkit::OutputOptions output = std::move(configured_output).value();

    // Starting from the moments of the sample rather than from the configured
    // analytical moments is what makes the runs comparable with the Monte Carlo
    // run: all begin from the same finite population. Only the moments of the
    // case are replaced; its node counts and options are kept.
    const std::optional<gqmom_1d_droplet_breakup_with_mc::SamplingSetup> sampled =
        gqmom_1d_droplet_breakup_with_mc::sample_initial_distribution(configuration_file);
    if (!sampled)
    {
        return 1;
    }
    gqmom_case.moments = sampled->moments.moments;

    // The same case with classical QMOM, for comparison: N nodes from the same
    // moments instead of the N_G nodes of GQMOM.
    qmomkit::QmomCase qmom_case = gqmom_case;
    qmom_case.qmom_options.family = qmomkit::QmomFamily::classical;

    // The model is one moment source, defined in physics.cpp, which both states
    // evaluate on the nodes of their own quadratures.
    const qmomkit::QmomMomentSource source = gqmom_1d_droplet_breakup_with_mc::binary_breakup_source();

    try
    {
        // initial_state() throws when the sampled moments admit no quadrature,
        // and the writers throw when a file cannot be written.
        qmomkit::QmomState gqmom = gqmom_case.initial_state();
        qmomkit::QmomState qmom = qmom_case.initial_state();
        const std::vector<double> gqmom_initial = gqmom.moment_values();
        const std::vector<double> qmom_initial = qmom.moment_values();

        std::optional<qmomkit::MomentWriter> gqmom_moments;
        std::optional<qmomkit::MomentWriter> qmom_moments;
        std::optional<qmomkit::QuadratureWriter> gqmom_quadrature;
        std::optional<qmomkit::QuadratureWriter> qmom_quadrature;
        if (output.write_moments || output.write_quadrature)
        {
            std::filesystem::create_directories(output.directory);
        }
        if (output.write_moments)
        {
            gqmom_moments.emplace(qmomkit::MomentWriter::for_qmom(output.directory / "gqmom_moment_history.txt", gqmom, output));
            qmom_moments.emplace(qmomkit::MomentWriter::for_qmom(output.directory / "qmom_moment_history.txt", qmom, output));
        }
        if (output.write_quadrature)
        {
            gqmom_quadrature.emplace(
                qmomkit::QuadratureWriter::for_qmom(output.directory / "gqmom_quadrature_history.txt", gqmom.maximum_node_count()));
            qmom_quadrature.emplace(qmomkit::QuadratureWriter::for_qmom(output.directory / "qmom_quadrature_history.txt", qmom.maximum_node_count()));
        }
        const auto write_output = [&](const double time)
        {
            if (gqmom_moments)
            {
                gqmom_moments->write(time, gqmom, gqmom.quadrature());
                qmom_moments->write(time, qmom, qmom.quadrature());
            }
            if (gqmom_quadrature)
            {
                gqmom_quadrature->write(time, gqmom.quadrature());
                qmom_quadrature->write(time, qmom.quadrature());
            }
        };
        write_output(0.0);

        // Write every output interval of [output], and at the end without one.
        const int steps_per_output = output.interval ? static_cast<int>(std::lround(*output.interval / time_step)) : step_count;
        for (int step = 1; step <= step_count; ++step)
        {
            // One source step of each state. There is no node advection: breakup
            // changes the number of droplets, not the volume of a droplet between
            // breakups.
            const double time = (step - 1) * time_step;
            qmomkit::QmomPhaseSpaceEvolutionResult gqmom_result = qmomkit::advance_qmom_moment_source(gqmom, time, time_step, source);
            if (!gqmom_result.ok())
            {
                std::cerr << "GQMOM source step failed: " << gqmom_result.message << '\n';
                return 1;
            }
            qmomkit::QmomPhaseSpaceEvolutionResult qmom_result = qmomkit::advance_qmom_moment_source(qmom, time, time_step, source);
            if (!qmom_result.ok())
            {
                std::cerr << "QMOM source step failed: " << qmom_result.message << '\n';
                return 1;
            }
            gqmom = std::move(gqmom_result.state);
            qmom = std::move(qmom_result.state);

            if (step % steps_per_output == 0 || step == step_count)
            {
                write_output(step * time_step);
            }
        }

        // Classical QMOM evolves only m0..m(2N-1). The state also holds m(2N), but
        // as the moment of its N nodes, so it is not printed.
        const auto evolved = [&qmom_case](const std::vector<double>& moments)
        { return std::vector<double>(moments.begin(), moments.begin() + static_cast<std::ptrdiff_t>(2 * qmom_case.node_count)); };
        gqmom_1d_droplet_breakup_with_mc::print_moments(std::cout, "GQMOM run", gqmom_initial, gqmom.moment_values());
        gqmom_1d_droplet_breakup_with_mc::print_moments(std::cout, "QMOM run", evolved(qmom_initial), evolved(qmom.moment_values()));
        if (gqmom_moments || gqmom_quadrature)
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
