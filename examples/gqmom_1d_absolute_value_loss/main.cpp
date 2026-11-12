#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"

#include <cmath>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const int step_count = 200;
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
    const qmomkit::QmomCase gqmom_case = std::move(configured_case).value();
    const qmomkit::OutputOptions output = std::move(configured_output).value();

    // The same case with classical QMOM, for comparison: N nodes from the same
    // moments instead of the N_G nodes of GQMOM.
    qmomkit::QmomCase qmom_case = gqmom_case;
    qmom_case.qmom_options.family = qmomkit::QmomFamily::classical;

    // The model is one moment source, defined in physics.cpp, which both
    // states evaluate on the nodes of their own quadratures.
    const qmomkit::QmomMomentSource source = gqmom_1d_absolute_value_loss::absolute_value_loss();

    try
    {
        // initial_state() throws when the configured moments admit no quadrature,
        // and the writers throw when a file cannot be written.
        qmomkit::QmomState gqmom = gqmom_case.initial_state();
        qmomkit::QmomState qmom = qmom_case.initial_state();

        // The exact solution starts from the normal density of the configuration,
        // whose number density, mean and variance are those of the initial moments.
        const std::vector<double> initial = gqmom.moment_values();
        const double number_density = initial[0];
        const double mean = initial[1] / number_density;
        const double standard_deviation = std::sqrt((initial[2] / number_density) - (mean * mean));

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
            gqmom_moments.emplace(qmomkit::MomentWriter::for_qmom(output.directory / "gqmom_moments.txt", gqmom, output));
            qmom_moments.emplace(qmomkit::MomentWriter::for_qmom(output.directory / "qmom_moments.txt", qmom, output));
        }
        if (output.write_quadrature)
        {
            gqmom_quadrature.emplace(qmomkit::QuadratureWriter::for_qmom(output.directory / "gqmom_quadrature.txt", gqmom.maximum_node_count()));
            qmom_quadrature.emplace(qmomkit::QuadratureWriter::for_qmom(output.directory / "qmom_quadrature.txt", qmom.maximum_node_count()));
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
            // One source step of each state. There is no node advection: the
            // loss changes the number of elements, not their coordinates.
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
        const double final_time = step_count * time_step;

        // The number density m0 and the mean of xi^2, m2/m0, of the exact
        // solution and of both states, at the start and at the end.
        const gqmom_1d_absolute_value_loss::ExactMoments exact =
            gqmom_1d_absolute_value_loss::exact_moments(number_density, mean, standard_deviation, final_time);
        const auto print_row = [](const std::string& label, const double m0, const double m2)
        { std::cout << "  " << std::left << std::setw(8) << label << std::setw(12) << m0 << m2 / m0 << '\n'; };
        std::cout << "t = 0: m0 " << number_density << ", m2/m0 " << initial[2] / number_density << '\n';
        std::cout << "t = " << final_time << ":\n";
        std::cout << "          m0          m2/m0\n";
        print_row("exact", exact.m0, exact.m2);
        print_row("GQMOM", gqmom.moment_values()[0], gqmom.moment_values()[2]);
        print_row("QMOM", qmom.moment_values()[0], qmom.moment_values()[2]);
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
