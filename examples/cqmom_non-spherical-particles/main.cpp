#include "example_workflow.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <utility>
#include <vector>

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";

    try
    {
        std::optional<cqmom_non_spherical_particles::PreparedCqmomExample> prepared =
            cqmom_non_spherical_particles::prepare_cqmom_example(configuration_file);
        if (!prepared.has_value())
        {
            return 1;
        }
        cqmom_non_spherical_particles::PreparedCqmomExample example = std::move(*prepared);
        qmomkit::CqmomState state = std::move(example.state);
        qmomkit::CqmomQuadrature quadrature = example.initial_quadrature;
        cqmom_non_spherical_particles::print_transport_settings(example);

        std::vector<std::pair<double, qmomkit::CqmomQuadrature>> quadrature_history{{cqmom_non_spherical_particles::cqmom_initial_time, quadrature}};

        const double output_interval = example.output_options.interval.value_or(cqmom_non_spherical_particles::cqmom_final_time -
                                                                                cqmom_non_spherical_particles::cqmom_initial_time);
        double time = cqmom_non_spherical_particles::cqmom_initial_time;
        double next_output_time = cqmom_non_spherical_particles::cqmom_initial_time + output_interval;
        // Advance one macro step with the integration and operator-splitting schemes
        // selected in configuration.qmom. The continuous operator advances particle
        // trajectories, while the discontinuous operator applies binary fragmentation.
        while (time < cqmom_non_spherical_particles::cqmom_final_time - 1.0e-12)
        {
            const double current_step =
                std::min(cqmom_non_spherical_particles::cqmom_time_step, cqmom_non_spherical_particles::cqmom_final_time - time);
            const qmomkit::CqmomTransportResult transported =
                qmomkit::advance_transport_cqmom(state, time, current_step, example.velocity, example.source, example.transport_options);
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
            cqmom_non_spherical_particles::print_bounded_support_event_summaries(transported.diagnostics, time);

            if (time + 1.0e-12 >= next_output_time || time >= cqmom_non_spherical_particles::cqmom_final_time - 1.0e-12)
            {
                quadrature_history.emplace_back(time, quadrature);
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
            }
        }

        cqmom_non_spherical_particles::write_cqmom_outputs_and_summary(example, quadrature_history, quadrature,
                                                                       cqmom_non_spherical_particles::cqmom_final_time);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
