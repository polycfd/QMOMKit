#include "local_rate_workflow.hpp"

#include "physics.hpp"
#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

int qmom_1d_breakage_advection_mc_local_rates::run_local_rate_workflow()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.01;

    const std::optional<std::string> configuration_text = qmom_1d_breakage_advection_mc_local_rates::read_configuration_text(configuration_file);
    if (!configuration_text.has_value())
    {
        return 1;
    }

    auto configured = qmomkit::config::load_qmom(configuration_file);
    if (!configured)
    {
        std::cerr << configured.error() << '\n';
        return 1;
    }

    auto output_configured = qmomkit::config::load_output_options(configuration_file);
    if (!output_configured)
    {
        std::cerr << output_configured.error() << '\n';
        return 1;
    }

    auto phase_space_evolution_configured = qmomkit::config::load_qmom_phase_space_evolution_options(configuration_file);
    if (!phase_space_evolution_configured)
    {
        std::cerr << phase_space_evolution_configured.error() << '\n';
        return 1;
    }

    const qmomkit::QmomCase qmom_case = std::move(configured).value();
    const qmomkit::OutputOptions output_options = std::move(output_configured).value();
    qmomkit::QmomPhaseSpaceEvolutionOptions phase_space_evolution_options = std::move(phase_space_evolution_configured).value();
    if (phase_space_evolution_options.scheme != qmomkit::PhaseSpaceEvolutionScheme::strang)
    {
        std::cerr << "The local-rate example currently supports only the strang phase-space evolution scheme.\n";
        return 1;
    }
    auto initial_condition = qmomkit::config::parse_qmom_initial_condition(*configuration_text, configuration_file.string());
    if (!initial_condition)
    {
        std::cerr << initial_condition.error() << '\n';
        return 1;
    }

    const std::optional<qmom_1d_breakage_advection_mc_local_rates::SamplingSetup> sampled =
        qmom_1d_breakage_advection_mc_local_rates::sample_initial_distribution(std::move(initial_condition).value());
    if (!sampled.has_value())
    {
        return 1;
    }
    if (!sampled->sampling.ok())
    {
        std::cerr << sampled->sampling.message << '\n';
        return 1;
    }

    const qmomkit::monte_carlo::OneDimensionalMomentEstimate& sampled_moments = sampled->moments;
    if (!sampled_moments.ok())
    {
        std::cerr << sampled_moments.message << '\n';
        return 1;
    }

    const std::vector<double> initial_moments = sampled_moments.moments.values();
    std::vector<double> moments = initial_moments;

    const qmomkit::QmomCase sampled_case{sampled_moments.moments, qmom_case.node_count, qmom_case.options};
    phase_space_evolution_options.qmom_options = sampled_case.options;
    const qmomkit::QmomResult initial_result = qmomkit::compute_qmom(sampled_case);
    if (!initial_result.ok())
    {
        std::cerr << "Initial QMOM inversion failed: " << initial_result.message << '\n';
        std::cerr << "  QmomStatus: " << qmom_1d_breakage_advection_mc_local_rates::to_string(initial_result.status) << '\n';
        return 1;
    }

    qmomkit::Quadrature1D quadrature{initial_result.abscissas, initial_result.weights};

    // In a CG/DG finite-element code, the application would own spatial
    // coefficients and assemble a local raw-moment vector at each quadrature
    // point. This homogeneous example keeps only that local vector and lets
    // the surrounding application code advance it.
    const qmomkit::PhaseSpaceVelocity velocity = qmom_1d_breakage_advection_mc_local_rates::advection_velocity();
    const qmomkit::MomentSource source = qmom_1d_breakage_advection_mc_local_rates::binary_fragmentation_source();
    const qmomkit::MomentBoundaryContribution boundary = qmom_1d_breakage_advection_mc_local_rates::boundary_contribution();

    const auto moments_are_invertible = [&](const std::vector<double>& candidate)
    { return qmomkit::compute_qmom(qmomkit::MomentSet(candidate), qmom_case.node_count, qmom_case.options).ok(); };

    qmomkit::SspRk3Options rate_step_options = phase_space_evolution_options.source_options.ssp_rk3_options;
    // Reject a tentative SSP-RK3 state if QMOM cannot reconstruct a
    // quadrature from its moments. This moment-source path uses SSP-RK3 retry
    // controls, so a rejected state fails the substep. A CG/DG solver
    // may apply additional limiters or admissibility checks.
    rate_step_options.stage_accept = [&](qmomkit::SspRk3Stage, const std::vector<double>& candidate) { return moments_are_invertible(candidate); };

    // QMOMKit returns the instantaneous raw-moment rate dM/dt.
    // The application remains responsible for advancing the moments in time.

    // Advection-induced moment rate
    const qmomkit::OdeRhs advection_rhs = [&](const double stage_time, const std::vector<double>& stage_moments)
    {
        const qmomkit::QmomAdvectionMomentRateResult rates = qmomkit::evaluate_qmom_advection_moment_rates(
            stage_moments, sampled_case.node_count, stage_time, velocity, phase_space_evolution_options.qmom_options);
        if (!rates.ok())
        {
            std::ostringstream message;
            message << "Advection moment-rate evaluation failed at t = " << stage_time
                    << ": advection-rate status = " << qmom_1d_breakage_advection_mc_local_rates::to_string(rates.status)
                    << ", QmomStatus = " << qmom_1d_breakage_advection_mc_local_rates::to_string(rates.inversion_status)
                    << ", message = " << rates.message;
            throw std::runtime_error(message.str());
        }
        return rates.rates;
    };

    // Moment-source rate
    // The source callback already returns raw-moment rates, so this evaluation does not require a QMOM inversion.
    const qmomkit::OdeRhs source_rhs = [&](const double stage_time, const std::vector<double>& stage_moments)
    {
        const qmomkit::QmomMomentSourceRateResult rates = qmomkit::evaluate_qmom_moment_source_rates(stage_time, stage_moments, source, boundary);
        if (!rates.ok())
        {
            std::ostringstream message;
            message << "Moment-source rate evaluation failed at t = " << stage_time
                    << ": moment-source rate status = " << qmom_1d_breakage_advection_mc_local_rates::to_string(rates.status)
                    << ", message = " << rates.message;
            throw std::runtime_error(message.str());
        }
        return rates.rates;
    };

    qmom_1d_breakage_advection_mc_local_rates::OutputFiles output_files;
    try
    {
        output_files = qmom_1d_breakage_advection_mc_local_rates::open_output_files(output_options, initial_time, initial_moments, quadrature);
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

        // Advance one time step:
        //   1. advection evolution over dt/2,
        //   2. moment-source evolution over dt,
        //   3. advection evolution over dt/2.
        qmomkit::SspRk3Result first_advection;
        try
        {
            first_advection = qmomkit::ssp_rk3_step(moments, time, half_step, advection_rhs, rate_step_options);
        }
        catch (const std::exception& error)
        {
            std::cerr << "First advection local-rate half-step failed: " << error.what() << '\n';
            return 1;
        }
        if (!first_advection.ok())
        {
            std::cerr << "First advection local-rate half-step failed: SSP-RK3 status = "
                      << qmom_1d_breakage_advection_mc_local_rates::to_string(first_advection.status) << ", retries = " << first_advection.retries
                      << ", attempted time step = " << first_advection.time_step_used << '\n';
            return 1;
        }
        moments = first_advection.state;

        qmomkit::SspRk3Result source_update;
        try
        {
            source_update = qmomkit::ssp_rk3_step(moments, time + half_step, current_step, source_rhs, rate_step_options);
        }
        catch (const std::exception& error)
        {
            std::cerr << "Moment-source local stage failed: " << error.what() << '\n';
            return 1;
        }
        if (!source_update.ok())
        {
            std::cerr << "Moment-source local stage failed: SSP-RK3 status = "
                      << qmom_1d_breakage_advection_mc_local_rates::to_string(source_update.status) << ", retries = " << source_update.retries
                      << ", attempted time step = " << source_update.time_step_used << '\n';
            return 1;
        }
        moments = source_update.state;

        qmomkit::SspRk3Result second_advection;
        try
        {
            second_advection = qmomkit::ssp_rk3_step(moments, time + half_step, half_step, advection_rhs, rate_step_options);
        }
        catch (const std::exception& error)
        {
            std::cerr << "Second advection local-rate half-step failed: " << error.what() << '\n';
            return 1;
        }
        if (!second_advection.ok())
        {
            std::cerr << "Second advection local-rate half-step failed: SSP-RK3 status = "
                      << qmom_1d_breakage_advection_mc_local_rates::to_string(second_advection.status) << ", retries = " << second_advection.retries
                      << ", attempted time step = " << second_advection.time_step_used << '\n';
            return 1;
        }
        moments = second_advection.state;
        time += current_step;

        // Moments remain the evolved state. Reconstruct quadrature only when
        // weights and abscissas are needed for output.
        const qmomkit::QmomResult output_quadrature = qmomkit::compute_qmom(qmomkit::MomentSet(moments), qmom_case.node_count, qmom_case.options);
        if (!output_quadrature.ok())
        {
            std::cerr << "QMOM inversion for quadrature output failed at t = " << time << ": " << output_quadrature.message << '\n';
            std::cerr << "  QmomStatus: " << qmom_1d_breakage_advection_mc_local_rates::to_string(output_quadrature.status) << '\n';
            return 1;
        }
        quadrature = qmomkit::Quadrature1D{output_quadrature.abscissas, output_quadrature.weights};

        if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
        {
            try
            {
                qmom_1d_breakage_advection_mc_local_rates::write_output_state(output_files, time, moments, quadrature);
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

    qmom_1d_breakage_advection_mc_local_rates::print_run_report(
        std::cout, configuration_file, output_options, qmom_1d_breakage_advection_mc_local_rates::initialization_sample_count,
        qmom_1d_breakage_advection_mc_local_rates::initialization_seed, initial_time, final_time, time_step, output_interval, initial_moments,
        qmomkit::Quadrature1D{initial_result.abscissas, initial_result.weights}, moments, quadrature, time);

    return 0;
}
