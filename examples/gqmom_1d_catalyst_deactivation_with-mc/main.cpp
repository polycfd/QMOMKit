#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const int step_count = 100;
    const double time_step = 0.1;

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
    const std::optional<gqmom_1d_catalyst_deactivation_with_mc::SamplingSetup> sampled =
        gqmom_1d_catalyst_deactivation_with_mc::sample_initial_distribution(configuration_file);
    if (!sampled)
    {
        return 1;
    }
    gqmom_case.moments = sampled->moments.moments;

    // The same case with classical QMOM, for comparison: N nodes from the same
    // moments instead of the N_G nodes of GQMOM.
    qmomkit::QmomCase qmom_case = gqmom_case;
    qmom_case.qmom_options.family = qmomkit::QmomFamily::classical;

    // The model, defined in physics.cpp, is the coking velocity, which moves
    // the nodes, and the removal source, which both states evaluate on the
    // nodes of their own quadratures.
    const qmomkit::PhaseSpaceVelocity velocity = gqmom_1d_catalyst_deactivation_with_mc::coking_velocity();
    const qmomkit::QmomMomentSource source = gqmom_1d_catalyst_deactivation_with_mc::removal_source();

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
            // One Strang-split step of each state: coking over half the step,
            // removal over the step, and coking over the other half.
            const double time = (step - 1) * time_step;
            qmomkit::QmomPhaseSpaceEvolutionResult gqmom_result = qmomkit::advance_qmom_phase_space(gqmom, time, time_step, velocity, source);
            if (!gqmom_result.ok())
            {
                std::cerr << "GQMOM step failed: " << gqmom_result.message << '\n';
                return 1;
            }
            qmomkit::QmomPhaseSpaceEvolutionResult qmom_result = qmomkit::advance_qmom_phase_space(qmom, time, time_step, velocity, source);
            if (!qmom_result.ok())
            {
                std::cerr << "QMOM step failed: " << qmom_result.message << '\n';
                return 1;
            }
            gqmom = std::move(gqmom_result.state);
            qmom = std::move(qmom_result.state);

            if (step % steps_per_output == 0 || step == step_count)
            {
                write_output(step * time_step);
            }
        }

        // The initial and final raw moments side by side, as main_mc.cpp prints
        // them. Classical QMOM evolves only m0..m(2N-1). The state also holds
        // m(2N), but as the moment of its N nodes, so it is not printed.
        const auto print_moments =
            [](const char* run, const std::vector<double>& initial_moments, const std::vector<double>& final_moments, const std::size_t moment_count)
        {
            std::cout << run << ", raw moments of the coked-fraction distribution:\n";
            std::cout << "      " << std::setw(14) << "initial" << std::setw(14) << "final" << '\n';
            for (std::size_t order = 0; order < moment_count; ++order)
            {
                std::cout << "  M" << order << "  " << std::setw(14) << initial_moments[order] << std::setw(14) << final_moments[order] << '\n';
            }
        };
        print_moments("GQMOM run", gqmom_initial, gqmom.moment_values(), gqmom_initial.size());
        print_moments("QMOM run", qmom_initial, qmom.moment_values(), 2 * qmom_case.node_count);
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
