#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.01;

    try
    {
        const auto configured = qmomkit::config::load_cqmom(configuration_file);
        if (!configured.ok())
        {
            cqmom_2d_continuous_breakage::print_configuration_errors(configured.errors);
            return 1;
        }
        const auto output_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_configured.ok())
        {
            cqmom_2d_continuous_breakage::print_configuration_errors(output_configured.errors);
            return 1;
        }
        const qmomkit::OutputOptions output_options = *output_configured.value;

        const std::optional<std::string> configuration_text = cqmom_2d_continuous_breakage::read_configuration_text(configuration_file);
        if (!configuration_text.has_value())
            return 1;

        const auto initial_condition = qmomkit::config::parse_cqmom_initial_condition(*configuration_text, configuration_file.string());
        if (!initial_condition.ok())
        {
            cqmom_2d_continuous_breakage::print_configuration_errors(initial_condition.errors);
            return 1;
        }

        const std::optional<qmomkit::monte_carlo::TensorSampledInitialization> sampled =
            cqmom_2d_continuous_breakage::sample_initial_distribution(*initial_condition.value);
        if (!sampled.has_value())
            return 1;
        if (!sampled->ok())
        {
            std::cerr << sampled->sampling.message << sampled->moments.message << '\n';
            return 1;
        }

        const std::size_t initial_section = configuration_text->find("[initial]");
        if (initial_section == std::string::npos)
        {
            std::cerr << "Could not find the initial condition section in " << configuration_file << ".\n";
            return 1;
        }

        std::ostringstream sampled_case_text;
        sampled_case_text << std::setprecision(17) << configuration_text->substr(0, initial_section) << "\n[moments]\nextents = ";
        for (std::size_t dimension = 0; dimension < sampled->moments.moments.extents.size(); ++dimension)
        {
            if (dimension > 0)
            {
                sampled_case_text << ", ";
            }
            sampled_case_text << sampled->moments.moments.extents[dimension];
        }
        sampled_case_text << '\n';
        for (std::size_t flat_index = 0; flat_index < sampled->moments.moments.values.size(); ++flat_index)
        {
            const qmomkit::MultiIndex order = qmomkit::unflatten_multi_index(flat_index, sampled->moments.moments.extents);
            sampled_case_text << "M[";
            for (std::size_t dimension = 0; dimension < order.size(); ++dimension)
            {
                if (dimension > 0)
                {
                    sampled_case_text << ',';
                }
                sampled_case_text << order[dimension];
            }
            sampled_case_text << "] = " << sampled->moments.moments.values[flat_index] << '\n';
        }

        const auto sampled_case_result = qmomkit::config::parse_cqmom(sampled_case_text.str(), "sampled-cqmom-initialization");
        if (!sampled_case_result.ok())
        {
            cqmom_2d_continuous_breakage::print_configuration_errors(sampled_case_result.errors);
            return 1;
        }
        const qmomkit::CqmomCase cqmom_case = *sampled_case_result.value;

        const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(cqmom_case);
        if (!cqmom.ok())
        {
            std::cerr << "CQMOM inversion failed: " << cqmom.message << '\n';
            return 1;
        }

        qmomkit::CqmomState state = cqmom_case.initial_state();
        qmomkit::CqmomQuadrature quadrature = cqmom.quadrature;
        const qmomkit::CqmomMomentLayout moment_output_layout = qmomkit::create_cqmom_moment_layout(cqmom_case.hierarchy());
        const qmomkit::RawMomentRequest moment_output_request{moment_output_layout.required_indices};
        const qmomkit::PhaseSpaceVelocity velocity = cqmom_2d_continuous_breakage::continuous_velocity();
        const qmomkit::CqmomMomentSource source = cqmom_2d_continuous_breakage::binary_fragmentation_source();

        qmomkit::CqmomStateTransportOptions transport_options;

        std::vector<std::pair<double, qmomkit::CqmomQuadrature>> quadrature_history{{initial_time, quadrature}};

        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output_time = initial_time + output_interval;
        while (time < final_time - 1.0e-12)
        {
            const double current_step = std::min(time_step, final_time - time);
            const qmomkit::CqmomStateTransportResult transported =
                qmomkit::advance_transport_cqmom(state, time, current_step, velocity, source, transport_options);
            if (!transported.ok())
            {
                std::cerr << "CQMOM transport failed: " << transported.message << '\n';
                return 1;
            }
            if (!transported.quadrature.has_value())
            {
                std::cerr << "CQMOM transport did not return a quadrature.\n";
                return 1;
            }

            state = transported.state;
            quadrature = *transported.quadrature;
            time = transported.final_time;

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                quadrature_history.emplace_back(time, quadrature);
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
            }
        }

        if (output_options.write_quadrature)
        {
            cqmom_2d_continuous_breakage::ensure_output_directory(output_options.directory);
        }
        if (output_options.write_moments)
        {
            cqmom_2d_continuous_breakage::ensure_output_directory(output_options.directory);
            qmomkit::MomentWriter::for_cqmom(output_options.directory / "initial_moments.txt", moment_output_request)
                .write(initial_time, cqmom.quadrature);

            qmomkit::MomentWriter moment_history_writer =
                qmomkit::MomentWriter::for_cqmom(output_options.directory / "moment_history.txt", moment_output_request);
            for (const auto& [history_time, history_quadrature] : quadrature_history)
            {
                moment_history_writer.write(history_time, history_quadrature);
            }
        }
        if (output_options.write_quadrature)
        {
            qmomkit::QuadratureWriter::for_cqmom(output_options.directory / "initial_quadrature.txt", cqmom.quadrature.number_of_dimensions(),
                                                 cqmom.quadrature.weights().size())
                .write(initial_time, cqmom.quadrature);

            qmomkit::QuadratureWriter quadrature_history_writer = qmomkit::QuadratureWriter::for_cqmom(
                output_options.directory / "quadrature_history.txt", quadrature_history.front().second.number_of_dimensions(),
                quadrature_history.front().second.weights().size());
            for (const auto& [history_time, history_quadrature] : quadrature_history)
            {
                quadrature_history_writer.write(history_time, history_quadrature);
            }
        }
        const qmomkit::SelectedRawMomentResult initial_moments = qmomkit::reconstruct_raw_moments(cqmom.quadrature, moment_output_request);
        if (!initial_moments.ok())
        {
            std::cerr << "Initial CQMOM moment reconstruction failed: " << initial_moments.message << '\n';
            return 1;
        }
        const qmomkit::SelectedRawMomentResult final_moments = qmomkit::reconstruct_raw_moments(quadrature, moment_output_request);
        if (!final_moments.ok())
        {
            std::cerr << "Final CQMOM moment reconstruction failed: " << final_moments.message << '\n';
            return 1;
        }
        cqmom_2d_continuous_breakage::print_result_summary(
            initial_moments, final_moments, quadrature, output_options.directory, cqmom_2d_continuous_breakage::initialization_sample_count,
            cqmom_2d_continuous_breakage::initialization_seed, time_step, final_time, output_options.write_moments, output_options.write_quadrature);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
