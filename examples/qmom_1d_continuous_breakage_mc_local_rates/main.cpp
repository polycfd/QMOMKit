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
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
constexpr double initial_time = 0.0;
constexpr double final_time = 1.0;
constexpr double time_step = 0.01;

class ContinuousRateFailure : public std::runtime_error
{
public:
    ContinuousRateFailure(const double time, const qmomkit::ContinuousMomentRateStatus status, const qmomkit::QmomStatus inversion_status,
                          std::string message)
        : std::runtime_error(message.empty() ? "continuous moment-rate evaluation failed" : std::move(message)), time_(time), status_(status),
          inversion_status_(inversion_status)
    {
    }

    [[nodiscard]] double time() const noexcept
    {
        return time_;
    }

    [[nodiscard]] qmomkit::ContinuousMomentRateStatus status() const noexcept
    {
        return status_;
    }

    [[nodiscard]] qmomkit::QmomStatus inversion_status() const noexcept
    {
        return inversion_status_;
    }

private:
    double time_ = 0.0;
    qmomkit::ContinuousMomentRateStatus status_ = qmomkit::ContinuousMomentRateStatus::success;
    qmomkit::QmomStatus inversion_status_ = qmomkit::QmomStatus::success;
};

class DiscontinuousRateFailure : public std::runtime_error
{
public:
    DiscontinuousRateFailure(const double time, const qmomkit::DiscontinuousMomentRateStatus status, std::string message)
        : std::runtime_error(message.empty() ? "discontinuous moment-rate evaluation failed" : std::move(message)), time_(time), status_(status)
    {
    }

    [[nodiscard]] double time() const noexcept
    {
        return time_;
    }

    [[nodiscard]] qmomkit::DiscontinuousMomentRateStatus status() const noexcept
    {
        return status_;
    }

private:
    double time_ = 0.0;
    qmomkit::DiscontinuousMomentRateStatus status_ = qmomkit::DiscontinuousMomentRateStatus::success;
};

[[nodiscard]] bool qmom_inversion_succeeds(const std::vector<double>& candidate, const std::size_t node_count,
                                           const qmomkit::QmomOptions& qmom_options)
{
    const qmomkit::QmomResult inversion = qmomkit::compute_qmom(qmomkit::MomentSet(candidate), node_count, qmom_options);
    return inversion.ok();
}

[[nodiscard]] qmomkit::SspRk3Options make_application_stage_policy(const std::size_t node_count, const qmomkit::QmomOptions& qmom_options)
{
    qmomkit::SspRk3Options options;

    // This example uses successful QMOM inversion as one application-selected
    // stage and final acceptance policy. That keeps later continuous-rate calls
    // on invertible moment states, but it is not built into the rate evaluators
    // and is not a universal definition of acceptability. A spatial solver
    // might add bounds, limiting, or other coefficient-level checks here.
    options.stage_accept = [node_count, qmom_options](qmomkit::SspRk3Stage, const std::vector<double>& candidate)
    { return qmom_inversion_succeeds(candidate, node_count, qmom_options); };

    return options;
}

[[nodiscard]] std::optional<std::vector<double>> run_local_ssp_step(const std::vector<double>& moments, const double time,
                                                                    const double requested_time_step, const qmomkit::OdeRhs& rhs,
                                                                    const qmomkit::SspRk3Options& options, const std::string& operation)
{
    try
    {
        const qmomkit::SspRk3Result step = qmomkit::ssp_rk3_step(moments, time, requested_time_step, rhs, options);
        if (!step.ok())
        {
            qmom_1d_continuous_breakage_mc_local_rates::print_ssp_failure(step, operation, requested_time_step);
            return std::nullopt;
        }
        return step.state;
    }
    catch (const ContinuousRateFailure& error)
    {
        qmom_1d_continuous_breakage_mc_local_rates::print_continuous_rate_failure(operation, error.time(), error.status(), error.inversion_status(),
                                                                                  error.what());
    }
    catch (const DiscontinuousRateFailure& error)
    {
        qmom_1d_continuous_breakage_mc_local_rates::print_discontinuous_rate_failure(operation, error.time(), error.status(), error.what());
    }
    catch (const std::exception& error)
    {
        std::cerr << "SSP-RK3 step threw while integrating " << operation << ":\n";
        std::cerr << "  message: " << error.what() << '\n';
        std::cerr << "  requested time step: " << requested_time_step << '\n';
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<qmomkit::Quadrature1D> reconstruct_for_output(const std::vector<double>& moments, const std::size_t node_count,
                                                                          const qmomkit::QmomOptions& qmom_options, const double time)
{
    const qmomkit::QmomResult result = qmomkit::compute_qmom(qmomkit::MomentSet(moments), node_count, qmom_options);
    if (!result.ok())
    {
        std::cerr << "QMOM inversion for quadrature output failed at t = " << time << ": " << result.message << '\n';
        std::cerr << "  QmomStatus: " << qmom_1d_continuous_breakage_mc_local_rates::to_string(result.status) << '\n';
        return std::nullopt;
    }
    return qmomkit::Quadrature1D{result.abscissas, result.weights};
}

} // namespace

int main()
{
    const std::filesystem::path configuration_file = "configuration.qmom";

    const std::optional<std::string> configuration_text = qmom_1d_continuous_breakage_mc_local_rates::read_configuration_text(configuration_file);
    if (!configuration_text.has_value())
    {
        return 1;
    }

    const auto configured = qmomkit::config::load_qmom(configuration_file);
    if (!configured.ok())
    {
        qmom_1d_continuous_breakage_mc_local_rates::print_configuration_errors(configured.errors);
        return 1;
    }

    const auto output_configured = qmomkit::config::load_output_options(configuration_file);
    if (!output_configured.ok())
    {
        qmom_1d_continuous_breakage_mc_local_rates::print_configuration_errors(output_configured.errors);
        return 1;
    }

    const qmomkit::QmomCase qmom_case = *configured.value;
    const qmomkit::OutputOptions output_options = *output_configured.value;
    const auto initial_condition = qmomkit::config::parse_qmom_initial_condition(*configuration_text, configuration_file.string());
    if (!initial_condition.ok())
    {
        qmom_1d_continuous_breakage_mc_local_rates::print_configuration_errors(initial_condition.errors);
        return 1;
    }

    const std::optional<qmom_1d_continuous_breakage_mc_local_rates::SamplingSetup> sampled =
        qmom_1d_continuous_breakage_mc_local_rates::sample_initial_distribution(*initial_condition.value);
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

    const qmomkit::QmomResult initial_result = qmomkit::compute_qmom(qmomkit::MomentSet(initial_moments), qmom_case.node_count, qmom_case.options);
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
    const qmomkit::SspRk3Options rate_step_options = make_application_stage_policy(qmom_case.node_count, qmom_case.options);

    // The continuous local-rate API returns an instantaneous raw-moment rate.
    // It inverts the supplied moments, evaluates the phase-space velocity at
    // quadrature nodes, and reconstructs moment rates. It does not advance
    // time; this example converts a failed public result into the RHS failure
    // mechanism used by the application-owned SSP-RK3 step.
    const qmomkit::OdeRhs continuous_rhs = [&](const double stage_time, const std::vector<double>& stage_moments)
    {
        const qmomkit::QmomContinuousMomentRateResult rates =
            qmomkit::evaluate_qmom_continuous_moment_rates(stage_moments, qmom_case.node_count, stage_time, velocity, qmom_case.options);
        if (!rates.ok())
        {
            throw ContinuousRateFailure(stage_time, rates.status, rates.inversion_status, rates.message);
        }
        return rates.rates;
    };

    // The discontinuous local-rate API calls a source model that already
    // returns raw-moment rates. No quadrature inversion is part of this rate
    // evaluation; the application still owns the RK stages and failure
    // propagation around it.
    const qmomkit::OdeRhs source_rhs = [&](const double stage_time, const std::vector<double>& stage_moments)
    {
        const qmomkit::QmomDiscontinuousMomentRateResult rates = qmomkit::evaluate_qmom_discontinuous_moment_rates(stage_time, stage_moments, source);
        if (!rates.ok())
        {
            throw DiscontinuousRateFailure(stage_time, rates.status, rates.message);
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

        // QMOMKit supplies local rates; this driver owns the SSP-RK3 substeps.
        // A finite-element code could use the same RHS adapters inside a
        // coefficient-level time integrator at each spatial quadrature point.
        //
        // The established continuous/source sequence is:
        //   continuous half-step -> discontinuous full step -> continuous half-step.
        // Each accepted substep updates the application-owned moment vector.
        std::optional<std::vector<double>> updated =
            run_local_ssp_step(moments, time, half_step, continuous_rhs, rate_step_options, "the first continuous local-rate half-step");
        if (!updated.has_value())
        {
            return 1;
        }
        moments = std::move(*updated);

        updated = run_local_ssp_step(moments, time + half_step, current_step, source_rhs, rate_step_options, "the discontinuous local source step");
        if (!updated.has_value())
        {
            return 1;
        }
        moments = std::move(*updated);

        updated =
            run_local_ssp_step(moments, time + half_step, half_step, continuous_rhs, rate_step_options, "the second continuous local-rate half-step");
        if (!updated.has_value())
        {
            return 1;
        }
        moments = std::move(*updated);
        time += current_step;

        // Moments remain the evolved state. Reconstruct quadrature only when
        // weights and abscissas are needed for output.
        const std::optional<qmomkit::Quadrature1D> output_quadrature = reconstruct_for_output(moments, qmom_case.node_count, qmom_case.options, time);
        if (!output_quadrature.has_value())
        {
            return 1;
        }
        quadrature = *output_quadrature;

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
