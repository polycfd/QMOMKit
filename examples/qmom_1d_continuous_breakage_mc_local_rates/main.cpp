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

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";
    const double initial_time = 0.0;
    const double final_time = 1.0;
    const double time_step = 0.01;

    const std::optional<std::string> configuration_text = qmom_1d_continuous_breakage_mc_local_rates::read_configuration_text(configuration_file);
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

    auto transport_configured = qmomkit::config::load_qmom_transport_options(configuration_file);
    if (!transport_configured)
    {
        std::cerr << transport_configured.error() << '\n';
        return 1;
    }

    const qmomkit::QmomCase qmom_case = std::move(configured).value();
    const qmomkit::OutputOptions output_options = std::move(output_configured).value();
    qmomkit::QmomTransportOptions transport_options = std::move(transport_configured).value();
    auto initial_condition = qmomkit::config::parse_qmom_initial_condition(*configuration_text, configuration_file.string());
    if (!initial_condition)
    {
        std::cerr << initial_condition.error() << '\n';
        return 1;
    }

    const std::optional<qmom_1d_continuous_breakage_mc_local_rates::SamplingSetup> sampled =
        qmom_1d_continuous_breakage_mc_local_rates::sample_initial_distribution(std::move(initial_condition).value());
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
        qmom_1d_continuous_breakage_mc_local_rates::compute_sampled_moments(*sampled);
    if (!sampled_moments.ok())
    {
        std::cerr << sampled_moments.message << '\n';
        return 1;
    }

    const std::vector<double> initial_moments = sampled_moments.moments.values();
    std::vector<double> moments = initial_moments;

    const qmomkit::QmomCase sampled_case{sampled_moments.moments, qmom_case.node_count, qmom_case.options};
    transport_options.qmom_options = sampled_case.options;
    const qmomkit::QmomResult initial_result = qmomkit::compute_qmom(sampled_case);
    if (!initial_result.ok())
    {
        std::cerr << "Initial QMOM inversion failed: " << initial_result.message << '\n';
        std::cerr << "  QmomStatus: " << qmom_1d_continuous_breakage_mc_local_rates::to_string(initial_result.status) << '\n';
        return 1;
    }

    qmomkit::Quadrature1D quadrature{initial_result.abscissas, initial_result.weights};

    // In a CG/DG finite-element code, the application would own spatial
    // coefficients and assemble a local raw-moment vector at each quadrature
    // point. This homogeneous example keeps only that local vector and lets
    // the surrounding application code advance it.
    const qmomkit::PhaseSpaceVelocity velocity = qmom_1d_continuous_breakage_mc_local_rates::continuous_velocity();
    const qmomkit::MomentSource source = qmom_1d_continuous_breakage_mc_local_rates::binary_fragmentation_source();

    const auto moments_are_invertible = [&](const std::vector<double>& candidate)
    { return qmomkit::compute_qmom(qmomkit::MomentSet(candidate), qmom_case.node_count, qmom_case.options).ok(); };

    qmomkit::SspRk3Options rate_step_options = transport_options.source_options.ssp_rk3_options;
    // Reject a tentative SSP-RK3 state if QMOM cannot reconstruct a
    // quadrature from its moments. This example keeps the default SSP-RK3
    // retry controls, so a rejected state fails the substep. A CG/DG solver
    // may apply additional limiters or admissibility checks.
    rate_step_options.stage_accept = [&](qmomkit::SspRk3Stage, const std::vector<double>& candidate) { return moments_are_invertible(candidate); };

    // QMOMKit returns the instantaneous raw-moment rate dM/dt.
    // The application remains responsible for advancing the moments in time.

    // Continuous phase-space advection rate
    const qmomkit::OdeRhs continuous_rhs = [&](const double stage_time, const std::vector<double>& stage_moments)
    {
        const qmomkit::QmomContinuousMomentRateResult rates = qmomkit::evaluate_qmom_continuous_moment_rates(
            stage_moments, sampled_case.node_count, stage_time, velocity, transport_options.qmom_options);
        if (!rates.ok())
        {
            std::ostringstream message;
            message << "Continuous moment-rate evaluation failed at t = " << stage_time
                    << ": continuous-rate status = " << qmom_1d_continuous_breakage_mc_local_rates::to_string(rates.status)
                    << ", QmomStatus = " << qmom_1d_continuous_breakage_mc_local_rates::to_string(rates.inversion_status)
                    << ", message = " << rates.message;
            throw std::runtime_error(message.str());
        }
        return rates.rates;
    };

    // Discontinuous phase-space source rate
    // The source callback already returns raw-moment rates, so this evaluation does not require a QMOM inversion.
    const qmomkit::OdeRhs source_rhs = [&](const double stage_time, const std::vector<double>& stage_moments)
    {
        const qmomkit::QmomDiscontinuousMomentRateResult rates = qmomkit::evaluate_qmom_discontinuous_moment_rates(stage_time, stage_moments, source);
        if (!rates.ok())
        {
            std::ostringstream message;
            message << "Discontinuous moment-rate evaluation failed at t = " << stage_time
                    << ": discontinuous-rate status = " << qmom_1d_continuous_breakage_mc_local_rates::to_string(rates.status)
                    << ", message = " << rates.message;
            throw std::runtime_error(message.str());
        }
        return rates.rates;
    };

    qmom_1d_continuous_breakage_mc_local_rates::OutputFiles output_files;
    try
    {
        output_files = qmom_1d_continuous_breakage_mc_local_rates::open_output_files(output_options, initial_time, initial_moments, quadrature);
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
        //   1. continuous evolution over dt/2,
        //   2. discontinuous source evolution over dt,
        //   3. continuous evolution over dt/2.
        qmomkit::SspRk3Result first_continuous;
        try
        {
            first_continuous = qmomkit::ssp_rk3_step(moments, time, half_step, continuous_rhs, rate_step_options);
        }
        catch (const std::exception& error)
        {
            std::cerr << "First continuous local-rate half-step failed: " << error.what() << '\n';
            return 1;
        }
        if (!first_continuous.ok())
        {
            std::cerr << "First continuous local-rate half-step failed: SSP-RK3 status = "
                      << qmom_1d_continuous_breakage_mc_local_rates::to_string(first_continuous.status) << ", retries = " << first_continuous.retries
                      << ", attempted time step = " << first_continuous.time_step_used << '\n';
            return 1;
        }
        moments = first_continuous.state;

        qmomkit::SspRk3Result source_update;
        try
        {
            source_update = qmomkit::ssp_rk3_step(moments, time + half_step, current_step, source_rhs, rate_step_options);
        }
        catch (const std::exception& error)
        {
            std::cerr << "Discontinuous local source step failed: " << error.what() << '\n';
            return 1;
        }
        if (!source_update.ok())
        {
            std::cerr << "Discontinuous local source step failed: SSP-RK3 status = "
                      << qmom_1d_continuous_breakage_mc_local_rates::to_string(source_update.status) << ", retries = " << source_update.retries
                      << ", attempted time step = " << source_update.time_step_used << '\n';
            return 1;
        }
        moments = source_update.state;

        qmomkit::SspRk3Result second_continuous;
        try
        {
            second_continuous = qmomkit::ssp_rk3_step(moments, time + half_step, half_step, continuous_rhs, rate_step_options);
        }
        catch (const std::exception& error)
        {
            std::cerr << "Second continuous local-rate half-step failed: " << error.what() << '\n';
            return 1;
        }
        if (!second_continuous.ok())
        {
            std::cerr << "Second continuous local-rate half-step failed: SSP-RK3 status = "
                      << qmom_1d_continuous_breakage_mc_local_rates::to_string(second_continuous.status)
                      << ", retries = " << second_continuous.retries << ", attempted time step = " << second_continuous.time_step_used << '\n';
            return 1;
        }
        moments = second_continuous.state;
        time += current_step;

        // Moments remain the evolved state. Reconstruct quadrature only when
        // weights and abscissas are needed for output.
        const qmomkit::QmomResult output_quadrature = qmomkit::compute_qmom(qmomkit::MomentSet(moments), qmom_case.node_count, qmom_case.options);
        if (!output_quadrature.ok())
        {
            std::cerr << "QMOM inversion for quadrature output failed at t = " << time << ": " << output_quadrature.message << '\n';
            std::cerr << "  QmomStatus: " << qmom_1d_continuous_breakage_mc_local_rates::to_string(output_quadrature.status) << '\n';
            return 1;
        }
        quadrature = qmomkit::Quadrature1D{output_quadrature.abscissas, output_quadrature.weights};

        if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
        {
            try
            {
                qmom_1d_continuous_breakage_mc_local_rates::write_output_state(output_files, time, moments, quadrature);
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

    qmom_1d_continuous_breakage_mc_local_rates::print_run_report(
        std::cout, configuration_file, output_options, qmom_1d_continuous_breakage_mc_local_rates::initialization_sample_count,
        qmom_1d_continuous_breakage_mc_local_rates::initialization_seed, initial_time, final_time, time_step, output_interval, initial_moments,
        qmomkit::Quadrature1D{initial_result.abscissas, initial_result.weights}, moments, quadrature, time);

    return 0;
}
