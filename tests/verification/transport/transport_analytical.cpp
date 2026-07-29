#include "../../../src/internal/cqmom_case.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace qmomkit
{

struct TestDenseCqmomTransportOptions : CqmomTransportOptions
{
    ConditionalMomentOptions conditional_options;
    QmomOptions qmom_options;
    bool adaptive_branch_counts = false;
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::detailed;
    VandermondeSolver vandermonde_solver = VandermondeSolver::automatic;
};

struct TestLeanCqmomTransportOptions : CqmomTransportOptions
{
    ConditionalMomentOptions conditional_options;
    QmomOptions qmom_options;
    bool adaptive_branch_counts = false;
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::detailed;
    VandermondeSolver vandermonde_solver = VandermondeSolver::automatic;
};

struct TestLeanCqmomTransportState
{
    CqmomMomentData moments;
};

struct TestDenseCqmomTransportResult
{
    RawMomentTensor moments;
    std::optional<CqmomQuadrature> quadrature;
    TransportStatus status = TransportStatus::success;
    TransportDiagnostics diagnostics;
    double final_time = 0.0;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == TransportStatus::success;
    }
};

struct TestLeanCqmomTransportResult
{
    TestLeanCqmomTransportState state;
    std::optional<CqmomQuadrature> quadrature;
    TransportStatus status = TransportStatus::success;
    TransportDiagnostics diagnostics;
    double final_time = 0.0;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == TransportStatus::success;
    }
};

TestDenseCqmomTransportResult advance_transport_cqmom(const RawMomentTensor& moments, const PhaseSpaceHierarchy& hierarchy, const double t,
                                                      const double time_step, const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                      const TestDenseCqmomTransportOptions& options)
{
    CqmomDenseCase dense_case{hierarchy,
                              moments,
                              options.conditional_options,
                              options.qmom_options,
                              options.adaptive_branch_counts,
                              options.diagnostic_level,
                              options.vandermonde_solver};
    CqmomState state = internal::create_cqmom_state(dense_case);
    const CqmomTransportResult result =
        advance_transport_cqmom(state, t, time_step, velocity, source, static_cast<const CqmomTransportOptions&>(options));
    const CqmomState& result_state = result.ok() ? result.state : state;
    return TestDenseCqmomTransportResult{RawMomentTensor{moments.extents, result_state.moment_values()},
                                         result.quadrature,
                                         result.status,
                                         result.diagnostics,
                                         result.final_time,
                                         result.message};
}

TestLeanCqmomTransportResult advance_transport_cqmom(const TestLeanCqmomTransportState& lean_state, const double t, const double time_step,
                                                     const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                     const TestLeanCqmomTransportOptions& options)
{
    CqmomLeanCase lean_case{lean_state.moments,       options.conditional_options, options.qmom_options, options.adaptive_branch_counts,
                            options.diagnostic_level, options.vandermonde_solver};
    CqmomState state = internal::create_cqmom_state(lean_case);
    const CqmomTransportResult result =
        advance_transport_cqmom(state, t, time_step, velocity, source, static_cast<const CqmomTransportOptions&>(options));
    const CqmomState& result_state = result.ok() ? result.state : state;
    return TestLeanCqmomTransportResult{TestLeanCqmomTransportState{CqmomMomentData{lean_state.moments.layout, result_state.moment_values()}},
                                        result.quadrature,
                                        result.status,
                                        result.diagnostics,
                                        result.final_time,
                                        result.message};
}

} // namespace qmomkit

namespace
{

constexpr double commuting_a = 0.2;
constexpr double commuting_lambda_scale = 0.05;
constexpr double noncommuting_b = 0.45;
constexpr double noncommuting_gamma = 0.17;
constexpr double cqmom_a = 0.16;
constexpr double cqmom_c = -0.11;
constexpr double cqmom_lambda_scale = 0.025;
constexpr double cqmom_native_source_gamma = -0.12;
constexpr double cqmom_combined_source_gamma = -0.09;

struct ErrorMetrics
{
    double moment_absolute = 0.0;
    double moment_scaled = 0.0;
    double reconstruction_absolute = 0.0;
    double reconstruction_scaled = 0.0;
    double conservation = 0.0;
};

struct RunSummary
{
    qmomkit::TransportScheme scheme = qmomkit::TransportScheme::strang;
    double dt = 0.0;
    std::size_t steps = 0;
    double final_time = 0.0;
    double error_norm = 0.0;
    double observed_order = std::numeric_limits<double>::quiet_NaN();
    ErrorMetrics errors;
    std::size_t retry_count = 0;
    bool success = false;
    bool output_row = false;
};

struct VerificationResult
{
    bool passed = true;
    std::string message;
    std::vector<RunSummary> runs;
};

struct QmomStepRun
{
    qmomkit::QmomTransportResult result;
    bool diagnostics_successful = true;
    std::size_t retry_count = 0;
};

struct CqmomNativeSourceAudit
{
    bool passed = true;
    std::string message;
    std::size_t invocation_count = 0;
    bool checked_initial_context_values = false;
    bool checked_split_context_values = false;
};

void append_failure(VerificationResult& result, const std::string& message)
{
    result.passed = false;
    if (!result.message.empty())
    {
        result.message += '\n';
    }
    result.message += message;
}

double integer_power(const double base, const std::size_t exponent)
{
    double value = 1.0;
    for (std::size_t order = 0; order < exponent; ++order)
    {
        value *= base;
    }
    return value;
}

void record_audit_failure(CqmomNativeSourceAudit& audit, const std::string& message)
{
    audit.passed = false;
    if (!audit.message.empty())
    {
        audit.message += '\n';
    }
    audit.message += message;
}

bool contains_text(const std::string& text, const std::string& fragment)
{
    return text.find(fragment) != std::string::npos;
}

bool vectors_equal(const std::vector<double>& lhs, const std::vector<double>& rhs, const double tolerance)
{
    if (lhs.size() != rhs.size())
    {
        return false;
    }
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        if (std::abs(lhs[i] - rhs[i]) > tolerance)
        {
            return false;
        }
    }
    return true;
}

void update_errors(ErrorMetrics& errors, const std::vector<double>& computed, const std::vector<double>& exact)
{
    for (std::size_t i = 0; i < computed.size(); ++i)
    {
        const double absolute = std::abs(computed[i] - exact[i]);
        errors.moment_absolute = std::max(errors.moment_absolute, absolute);
        errors.moment_scaled = std::max(errors.moment_scaled, absolute / std::max(1.0, std::abs(exact[i])));
    }
    if (!computed.empty() && !exact.empty())
    {
        errors.conservation = std::abs(computed.front() - exact.front());
    }
}

void update_reconstruction_errors(ErrorMetrics& errors, const std::vector<double>& reconstructed, const std::vector<double>& exact)
{
    for (std::size_t i = 0; i < reconstructed.size(); ++i)
    {
        const double absolute = std::abs(reconstructed[i] - exact[i]);
        errors.reconstruction_absolute = std::max(errors.reconstruction_absolute, absolute);
        errors.reconstruction_scaled = std::max(errors.reconstruction_scaled, absolute / std::max(1.0, std::abs(exact[i])));
    }
}

double scaled_infinity_error(const std::vector<double>& computed, const std::vector<double>& exact)
{
    double error = 0.0;
    for (std::size_t i = 0; i < computed.size(); ++i)
    {
        error = std::max(error, std::abs(computed[i] - exact[i]) / std::max(1.0, std::abs(exact[i])));
    }
    return error;
}

void print_vector_comparison(const std::string& title, const std::vector<double>& computed, const std::vector<double>& exact)
{
    std::cout << title << '\n';
    for (std::size_t i = 0; i < computed.size(); ++i)
    {
        std::cout << "  m[" << i << "] computed=" << computed[i] << " exact=" << exact[i] << " error=" << computed[i] - exact[i] << '\n';
    }
}

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::QmomCase>& configured)
{
    std::cerr << configured.error() << '\n';
}

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::CqmomDenseCase>& configured)
{
    std::cerr << configured.error() << '\n';
}

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::CqmomLeanCase>& configured)
{
    std::cerr << configured.error() << '\n';
}

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::CqmomCase>& configured)
{
    std::cerr << configured.error() << '\n';
}

qmomkit::QmomTransportOptions qmom_options(const qmomkit::QmomOptions& qmom_options, const double rk_relative_tolerance,
                                           const double rk_absolute_tolerance, const double max_step,
                                           const qmomkit::TransportScheme scheme = qmomkit::TransportScheme::strang)
{
    qmomkit::QmomTransportOptions options;
    options.qmom_options = qmom_options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::adaptive_rk23;
    options.node_options.rk23_options.relative_tolerance = rk_relative_tolerance;
    options.node_options.rk23_options.absolute_tolerance = rk_absolute_tolerance;
    options.node_options.rk23_options.max_step = max_step;
    options.source_options.ssp_rk3_options.estimate_error = true;
    options.scheme = scheme;
    return options;
}

const char* transport_scheme_name(const qmomkit::TransportScheme scheme)
{
    switch (scheme)
    {
    case qmomkit::TransportScheme::strang:
        return "strang";
    case qmomkit::TransportScheme::cd:
        return "cd";
    case qmomkit::TransportScheme::dc:
        return "dc";
    }
    return "unknown";
}

qmomkit::TestDenseCqmomTransportOptions cqmom_options(const qmomkit::CqmomDenseCase& dense_case)
{
    qmomkit::TestDenseCqmomTransportOptions options;
    options.conditional_options = dense_case.conditional_options;
    options.qmom_options = dense_case.qmom_options;
    options.adaptive_branch_counts = dense_case.adaptive_branch_counts;
    options.diagnostic_level = dense_case.diagnostic_level;
    options.vandermonde_solver = dense_case.vandermonde_solver;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::adaptive_rk23;
    options.node_options.rk23_options.relative_tolerance = 1.0e-12;
    options.node_options.rk23_options.absolute_tolerance = 1.0e-14;
    options.node_options.rk23_options.max_step = 0.01;
    options.source_options.ssp_rk3_options.estimate_error = true;
    options.scheme = qmomkit::TransportScheme::strang;
    return options;
}

qmomkit::CqmomTransportOptions cqmom_state_options()
{
    qmomkit::CqmomTransportOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::adaptive_rk23;
    options.node_options.rk23_options.relative_tolerance = 1.0e-12;
    options.node_options.rk23_options.absolute_tolerance = 1.0e-14;
    options.node_options.rk23_options.max_step = 0.01;
    options.source_options.ssp_rk3_options.estimate_error = true;
    options.scheme = qmomkit::TransportScheme::strang;
    return options;
}

qmomkit::TestLeanCqmomTransportOptions cqmom_lean_options(const qmomkit::CqmomLeanCase& lean_case)
{
    qmomkit::TestLeanCqmomTransportOptions options;
    options.conditional_options = lean_case.conditional_options;
    options.qmom_options = lean_case.qmom_options;
    options.adaptive_branch_counts = lean_case.adaptive_branch_counts;
    options.diagnostic_level = lean_case.diagnostic_level;
    options.vandermonde_solver = lean_case.vandermonde_solver;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::adaptive_rk23;
    options.node_options.rk23_options.relative_tolerance = 1.0e-12;
    options.node_options.rk23_options.absolute_tolerance = 1.0e-14;
    options.node_options.rk23_options.max_step = 0.01;
    options.source_options.ssp_rk3_options.estimate_error = true;
    options.scheme = qmomkit::TransportScheme::strang;
    return options;
}

bool diagnostics_are_successful(const qmomkit::TransportDiagnostics& diagnostics)
{
    bool continuous_success = false;
    bool reconstruction_success = false;
    switch (diagnostics.scheme)
    {
    case qmomkit::TransportScheme::strang:
        continuous_success = diagnostics.first_continuous_status == qmomkit::NodeAdvectionStatus::success &&
                             diagnostics.second_continuous_status == qmomkit::NodeAdvectionStatus::success;
        reconstruction_success = diagnostics.first_reconstruction_diagnostics.success && diagnostics.final_reconstruction_diagnostics.success;
        break;
    case qmomkit::TransportScheme::cd:
        continuous_success = diagnostics.first_continuous_status == qmomkit::NodeAdvectionStatus::success;
        reconstruction_success = diagnostics.first_reconstruction_diagnostics.success;
        break;
    case qmomkit::TransportScheme::dc:
        continuous_success = diagnostics.second_continuous_status == qmomkit::NodeAdvectionStatus::success;
        reconstruction_success = diagnostics.final_reconstruction_diagnostics.success;
        break;
    }
    return diagnostics.failed_stage == qmomkit::TransportStage::none && diagnostics.initial_qmom_status == qmomkit::QmomStatus::success &&
           diagnostics.post_source_qmom_status == qmomkit::QmomStatus::success && continuous_success &&
           diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::success && diagnostics.source_diagnostics.accepted &&
           !diagnostics.source_diagnostics.rejected && diagnostics.source_retries == 0 && reconstruction_success;
}

bool cqmom_diagnostics_are_successful(const qmomkit::TransportDiagnostics& diagnostics)
{
    return diagnostics.failed_stage == qmomkit::TransportStage::none && diagnostics.initial_cqmom_status == qmomkit::CqmomStatus::success &&
           diagnostics.post_source_cqmom_status == qmomkit::CqmomStatus::success &&
           diagnostics.first_continuous_status == qmomkit::NodeAdvectionStatus::success &&
           diagnostics.second_continuous_status == qmomkit::NodeAdvectionStatus::success &&
           diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::success && diagnostics.source_diagnostics.accepted &&
           !diagnostics.source_diagnostics.rejected && diagnostics.source_retries == 0 && diagnostics.first_reconstruction_diagnostics.success &&
           diagnostics.final_reconstruction_diagnostics.success;
}

bool cqmom_diagnostics_request_solver(const std::vector<qmomkit::CqmomDiagnostics>& diagnostics, const std::size_t dimensions,
                                      const qmomkit::VandermondeSolver solver)
{
    bool saw_conditional_branch = false;
    for (const qmomkit::CqmomDiagnostics& diagnostic : diagnostics)
    {
        if (diagnostic.conditioning_level + 1 >= dimensions)
        {
            continue;
        }
        saw_conditional_branch = true;
        if (diagnostic.conditional_vandermonde_requested_solver != solver || diagnostic.conditional_vandermonde_selected_solver != solver)
        {
            return false;
        }
    }
    return saw_conditional_branch;
}

bool cqmom_diagnostics_enable_adaptive_branch_counts(const std::vector<qmomkit::CqmomDiagnostics>& diagnostics)
{
    if (diagnostics.empty())
    {
        return false;
    }
    for (const qmomkit::CqmomDiagnostics& diagnostic : diagnostics)
    {
        if (!diagnostic.adaptive_branch_counts_enabled)
        {
            return false;
        }
    }
    return true;
}

void verify_helper_based_common_cqmom_policies(const std::string& case_name, const qmomkit::TransportDiagnostics& diagnostics,
                                               const std::size_t dimensions, VerificationResult& verification)
{
    if (diagnostics.cqmom_diagnostic_level != qmomkit::CqmomDiagnosticLevel::basic)
    {
        append_failure(verification, case_name + " helper did not propagate CqmomDiagnosticLevel::basic.");
    }
    if (diagnostics.cqmom_vandermonde_solver != qmomkit::VandermondeSolver::dense)
    {
        append_failure(verification, case_name + " helper did not propagate VandermondeSolver::dense.");
    }
    if (!cqmom_diagnostics_enable_adaptive_branch_counts(diagnostics.initial_cqmom_diagnostics))
    {
        append_failure(verification, case_name + " helper did not propagate adaptive branch-count diagnostics.");
    }
    if (!cqmom_diagnostics_request_solver(diagnostics.initial_cqmom_diagnostics, dimensions, qmomkit::VandermondeSolver::dense))
    {
        append_failure(verification, case_name + " helper did not use the dense Vandermonde solver in conditional reconstruction.");
    }
    if (diagnostics.cqmom_node_count_transitions.empty())
    {
        append_failure(verification, case_name + " helper did not expose adaptive active-count transition diagnostics.");
    }
}

std::vector<double> commuting_exact(const std::vector<double>& initial, const double elapsed)
{
    // Commuting QMOM reference: affine node velocity dx/dt = a x and diagonal
    // moment source dM_k/dt = -lambda_k M_k act on each raw moment by scalar
    // exponentials, so the split and unsplit solutions coincide.
    std::vector<double> exact(initial.size(), 0.0);
    for (std::size_t order = 0; order < initial.size(); ++order)
    {
        const double lambda = commuting_lambda_scale * static_cast<double>(order);
        exact[order] = initial[order] * std::exp((static_cast<double>(order) * commuting_a - lambda) * elapsed);
    }
    return exact;
}

std::vector<double> commuting_source_exact(const std::vector<double>& initial, const double elapsed)
{
    std::vector<double> exact(initial.size(), 0.0);
    for (std::size_t order = 0; order < initial.size(); ++order)
    {
        const double lambda = commuting_lambda_scale * static_cast<double>(order);
        exact[order] = initial[order] * std::exp(-lambda * elapsed);
    }
    return exact;
}

std::vector<double> commuting_velocity(double, const std::vector<double>& xi)
{
    return {commuting_a * xi[0]};
}

std::vector<double> commuting_source(double, const std::vector<double>& moments)
{
    std::vector<double> rhs(moments.size(), 0.0);
    for (std::size_t order = 1; order < moments.size(); ++order)
    {
        rhs[order] = -commuting_lambda_scale * static_cast<double>(order) * moments[order];
    }
    return rhs;
}

std::vector<double> zero_qmom_source(double, const std::vector<double>& moments)
{
    return std::vector<double>(moments.size(), 0.0);
}

qmomkit::MomentBoundaryRate qmom_boundary_injection_removal(double, const std::vector<double>& moments)
{
    constexpr double injection_rate = 0.08;
    constexpr double injection_coordinate = 0.25;
    constexpr double removal_rate = -0.03;
    constexpr double removal_coordinate = 1.1;
    std::vector<double> rate(moments.size(), 0.0);
    for (std::size_t order = 0; order < rate.size(); ++order)
    {
        rate[order] = injection_rate * integer_power(injection_coordinate, order) + removal_rate * integer_power(removal_coordinate, order);
    }
    return qmomkit::MomentBoundaryRate{rate};
}

std::vector<double> qmom_boundary_rate_values(const std::size_t count)
{
    return qmom_boundary_injection_removal(0.0, std::vector<double>(count, 0.0)).values;
}

std::vector<double> qmom_constant_source_values(const std::size_t count)
{
    std::vector<double> rate(count, 0.0);
    for (std::size_t order = 0; order < count; ++order)
    {
        rate[order] = 0.01 * static_cast<double>(order + 1);
    }
    return rate;
}

std::vector<double> qmom_constant_source(double, const std::vector<double>& moments)
{
    return qmom_constant_source_values(moments.size());
}

std::vector<double> noncommuting_exact(const std::vector<double>& initial, const double elapsed)
{
    // Noncommuting QMOM reference: constant node velocity couples M_k to
    // M_{k-1}, while the source damps each order diagonally. The closed-form
    // triangular exponential below is independent of the transport splitting
    // implementation and exposes first- versus second-order behavior.
    const std::size_t count = initial.size();
    std::vector<double> diagonal(count, 0.0);
    for (std::size_t order = 0; order < count; ++order)
    {
        diagonal[order] = -noncommuting_gamma * static_cast<double>(order);
    }

    std::vector<std::vector<double>> coefficients(count, std::vector<double>(count, 0.0));
    coefficients[0][0] = initial[0];
    for (std::size_t order = 1; order < count; ++order)
    {
        double lower_sum = 0.0;
        for (std::size_t basis = 0; basis < order; ++basis)
        {
            coefficients[order][basis] =
                static_cast<double>(order) * noncommuting_b * coefficients[order - 1][basis] / (diagonal[basis] - diagonal[order]);
            lower_sum += coefficients[order][basis];
        }
        coefficients[order][order] = initial[order] - lower_sum;
    }

    std::vector<double> exact(count, 0.0);
    for (std::size_t order = 0; order < count; ++order)
    {
        for (std::size_t basis = 0; basis <= order; ++basis)
        {
            exact[order] += coefficients[order][basis] * std::exp(diagonal[basis] * elapsed);
        }
    }
    return exact;
}

std::vector<double> noncommuting_velocity(double, const std::vector<double>&)
{
    return {noncommuting_b};
}

std::vector<double> noncommuting_source(double, const std::vector<double>& moments)
{
    std::vector<double> rhs(moments.size(), 0.0);
    for (std::size_t order = 1; order < moments.size(); ++order)
    {
        rhs[order] = -noncommuting_gamma * static_cast<double>(order) * moments[order];
    }
    return rhs;
}

QmomStepRun run_qmom_steps(const std::vector<double>& initial, const std::size_t node_count, const double t0, const double dt,
                           const std::size_t steps, const qmomkit::PhaseSpaceVelocity& velocity, const qmomkit::MomentSource& source,
                           const qmomkit::QmomTransportOptions& options)
{
    std::vector<double> moments = initial;
    double time = t0;
    qmomkit::QmomTransportResult result;
    bool diagnostics_successful = true;
    std::size_t retry_count = 0;
    for (std::size_t step = 0; step < steps; ++step)
    {
        result = qmomkit::advance_transport_qmom(moments, node_count, time, dt, velocity, source, options);
        if (!result.ok())
        {
            retry_count += result.diagnostics.source_retries;
            return {result, diagnostics_successful, retry_count};
        }
        diagnostics_successful = diagnostics_successful && diagnostics_are_successful(result.diagnostics);
        retry_count += result.diagnostics.source_retries;
        moments = result.moments;
        time = result.final_time;
    }
    return {result, diagnostics_successful, retry_count};
}

RunSummary validate_qmom_run(VerificationResult& verification, const std::vector<double>& initial, const std::size_t node_count, const double t0,
                             const double dt, const std::size_t steps, const qmomkit::PhaseSpaceVelocity& velocity,
                             const qmomkit::MomentSource& source, const qmomkit::QmomTransportOptions& options, const std::vector<double>& exact,
                             const double moment_tolerance, const double reconstruction_tolerance)
{
    const QmomStepRun step_run = run_qmom_steps(initial, node_count, t0, dt, steps, velocity, source, options);
    const qmomkit::QmomTransportResult& result = step_run.result;
    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = dt;
    summary.steps = steps;
    summary.retry_count = step_run.retry_count;
    if (!result.ok())
    {
        summary.final_time = result.final_time;
        append_failure(verification, "transport integrator failed: " + result.message);
        return summary;
    }
    summary.final_time = result.final_time;
    summary.success = true;
    if (!step_run.diagnostics_successful)
    {
        summary.success = false;
        append_failure(verification, "transport diagnostics did not report successful completion of every substage.");
    }
    if (std::abs(result.final_time - (t0 + dt * static_cast<double>(steps))) > 1.0e-12)
    {
        summary.success = false;
        append_failure(verification, "transport final time does not match the requested final time.");
    }

    update_errors(summary.errors, result.moments, exact);
    summary.error_norm = scaled_infinity_error(result.moments, exact);
    const qmomkit::Quadrature1D final_quadrature{result.quadrature.abscissas, result.quadrature.weights, qmomkit::QuadratureStatus::success, ""};
    const std::vector<double> reconstructed = qmomkit::reconstruct_raw_moments(final_quadrature, exact.size());
    update_reconstruction_errors(summary.errors, reconstructed, exact);
    print_vector_comparison("Final QMOM moments", result.moments, exact);
    if (summary.errors.moment_scaled > moment_tolerance)
    {
        summary.success = false;
        append_failure(verification, "final QMOM moments exceed the analytical tolerance.");
    }
    if (summary.errors.reconstruction_scaled > reconstruction_tolerance)
    {
        summary.success = false;
        append_failure(verification, "final QMOM quadrature reconstruction exceeds the analytical tolerance.");
    }
    return summary;
}

void validate_commuting_component_checks(VerificationResult& verification, const qmomkit::QmomCase& case_)
{
    const qmomkit::QmomResult initial = qmomkit::compute_qmom(case_);
    if (!initial.ok())
    {
        append_failure(verification, "initial QMOM inversion failed during component checks.");
        return;
    }

    qmomkit::NodeAdvectionOptions node_options;
    node_options.scheme = qmomkit::NodeAdvectionScheme::adaptive_rk23;
    node_options.rk23_options.relative_tolerance = 1.0e-12;
    node_options.rk23_options.absolute_tolerance = 1.0e-14;
    node_options.rk23_options.max_step = 0.005;
    constexpr double elapsed = 0.03;
    const qmomkit::QmomNodeAdvectionResult advanced = qmomkit::advance_nodes(initial, 0.0, elapsed, commuting_velocity, node_options);
    if (!advanced.ok())
    {
        append_failure(verification, "continuous operator component check failed.");
        return;
    }
    const double scale = std::exp(commuting_a * elapsed);
    for (std::size_t node = 0; node < initial.weights.size(); ++node)
    {
        if (std::abs(advanced.quadrature.weights[node] - initial.weights[node]) > 1.0e-13 ||
            std::abs(advanced.quadrature.abscissas[node] - scale * initial.abscissas[node]) > 1.0e-10)
        {
            append_failure(verification, "continuous operator did not preserve weights while scaling abscissas.");
            break;
        }
    }

    constexpr double source_probe_step = 0.002;
    const qmomkit::MomentSourceUpdateResult sourced =
        qmomkit::advance_moments_with_source(case_.moments.values(), 0.0, source_probe_step, commuting_source, {});
    if (!sourced.ok())
    {
        append_failure(verification, "source operator component check failed.");
        return;
    }
    const std::vector<double> exact = commuting_source_exact(case_.moments.values(), source_probe_step);
    if (scaled_infinity_error(sourced.moments, exact) > 2.0e-9)
    {
        append_failure(verification, "source operator diagonal-law check exceeds tolerance.");
    }
}

VerificationResult run_commuting_qmom(const std::filesystem::path& configuration_path)
{
    // Verify QMOM transport when continuous advection and source evolution
    // commute. The reference is exact scalar exponential moment evolution, so
    // this case protects substage diagnostics, mass preservation, and the
    // high-accuracy path used when splitting error should be negligible.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_qmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load commuting QMOM configuration.");
        return verification;
    }

    const qmomkit::QmomCase& qmom_case = configured.value();
    validate_commuting_component_checks(verification, qmom_case);
    const std::vector<double> initial = qmom_case.moments.values();
    const double t0 = 0.0;
    const auto options = qmom_options(qmom_case.options, 1.0e-12, 1.0e-14, 0.01);

    verification.runs.push_back(validate_qmom_run(verification, initial, qmom_case.node_count, t0, 0.08, 1, commuting_velocity, commuting_source,
                                                  options, commuting_exact(initial, 0.08), 2.0e-6, 2.0e-6));
    verification.runs.push_back(validate_qmom_run(verification, initial, qmom_case.node_count, t0, 0.04, 6, commuting_velocity, commuting_source,
                                                  options, commuting_exact(initial, 0.24), 2.0e-6, 2.0e-6));

    const std::vector<double> refinement_steps{0.12, 0.06, 0.03};
    double previous_error = 0.0;
    for (const double dt : refinement_steps)
    {
        const std::size_t steps = static_cast<std::size_t>(std::lround(0.24 / dt));
        RunSummary summary = validate_qmom_run(verification, initial, qmom_case.node_count, t0, dt, steps, commuting_velocity, commuting_source,
                                               options, commuting_exact(initial, 0.24), 8.0e-6, 8.0e-6);
        summary.output_row = true;
        if (previous_error > 0.0 && summary.error_norm > 0.0)
        {
            summary.observed_order = std::log(previous_error / summary.error_norm) / std::log(2.0);
            if (summary.observed_order < 2.4)
            {
                append_failure(verification, "commuting refinement did not decrease at the expected high-order source-integration rate.");
            }
        }
        previous_error = summary.error_norm;
        verification.runs.push_back(summary);
    }
    return verification;
}

VerificationResult run_qmom_boundary_contribution(const std::filesystem::path& configuration_path)
{
    VerificationResult verification;
    const auto configured = qmomkit::config::load_qmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load QMOM boundary contribution configuration.");
        return verification;
    }

    const qmomkit::QmomCase& qmom_case = configured.value();
    const std::vector<double> initial = qmom_case.moments.values();
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    qmomkit::QmomTransportOptions options = qmom_options(qmom_case.options, 1.0e-12, 1.0e-14, 0.01, qmomkit::TransportScheme::dc);

    constexpr double boundary_dt = 0.04;
    const qmomkit::QmomTransportResult boundary_result = qmomkit::advance_transport_qmom(
        initial, qmom_case.node_count, 0.0, boundary_dt, zero_velocity, zero_qmom_source, qmom_boundary_injection_removal, options);
    RunSummary boundary_summary;
    boundary_summary.scheme = options.scheme;
    boundary_summary.dt = boundary_dt;
    boundary_summary.steps = 1;
    boundary_summary.final_time = boundary_result.final_time;
    boundary_summary.retry_count = boundary_result.diagnostics.source_retries;
    boundary_summary.success = boundary_result.ok() && diagnostics_are_successful(boundary_result.diagnostics);
    std::vector<double> boundary_exact = initial;
    const std::vector<double> boundary_rate = qmom_boundary_rate_values(initial.size());
    for (std::size_t i = 0; i < boundary_exact.size(); ++i)
    {
        boundary_exact[i] += boundary_dt * boundary_rate[i];
    }
    if (!boundary_result.ok())
    {
        append_failure(verification, "QMOM boundary injection/removal transport failed: " + boundary_result.message);
    }
    else
    {
        update_errors(boundary_summary.errors, boundary_result.moments, boundary_exact);
        const qmomkit::Quadrature1D final_quadrature{boundary_result.quadrature.abscissas, boundary_result.quadrature.weights,
                                                     qmomkit::QuadratureStatus::success, ""};
        update_reconstruction_errors(boundary_summary.errors, qmomkit::reconstruct_raw_moments(final_quadrature, boundary_exact.size()),
                                     boundary_exact);
        boundary_summary.error_norm = scaled_infinity_error(boundary_result.moments, boundary_exact);
        print_vector_comparison("Final QMOM boundary contribution moments", boundary_result.moments, boundary_exact);
        if (boundary_summary.errors.moment_scaled > 5.0e-11 || boundary_summary.errors.reconstruction_scaled > 5.0e-11)
        {
            boundary_summary.success = false;
            append_failure(verification, "QMOM boundary injection/removal moments exceeded the analytical tolerance.");
        }
    }
    verification.runs.push_back(boundary_summary);

    constexpr double combined_dt = 0.03;
    const qmomkit::QmomTransportResult combined_result = qmomkit::advance_transport_qmom(
        initial, qmom_case.node_count, 0.0, combined_dt, zero_velocity, qmom_constant_source, qmom_boundary_injection_removal, options);
    RunSummary combined_summary;
    combined_summary.scheme = options.scheme;
    combined_summary.dt = combined_dt;
    combined_summary.steps = 1;
    combined_summary.final_time = combined_result.final_time;
    combined_summary.retry_count = combined_result.diagnostics.source_retries;
    combined_summary.success = combined_result.ok() && diagnostics_are_successful(combined_result.diagnostics);
    std::vector<double> combined_exact = initial;
    const std::vector<double> source_rate = qmom_constant_source_values(initial.size());
    for (std::size_t i = 0; i < combined_exact.size(); ++i)
    {
        combined_exact[i] += combined_dt * (source_rate[i] + boundary_rate[i]);
    }
    if (!combined_result.ok())
    {
        append_failure(verification, "QMOM combined source plus boundary contribution transport failed: " + combined_result.message);
    }
    else
    {
        update_errors(combined_summary.errors, combined_result.moments, combined_exact);
        const qmomkit::Quadrature1D final_quadrature{combined_result.quadrature.abscissas, combined_result.quadrature.weights,
                                                     qmomkit::QuadratureStatus::success, ""};
        update_reconstruction_errors(combined_summary.errors, qmomkit::reconstruct_raw_moments(final_quadrature, combined_exact.size()),
                                     combined_exact);
        combined_summary.error_norm = scaled_infinity_error(combined_result.moments, combined_exact);
        print_vector_comparison("Final QMOM combined boundary/source moments", combined_result.moments, combined_exact);
        if (combined_summary.errors.moment_scaled > 5.0e-11 || combined_summary.errors.reconstruction_scaled > 5.0e-11)
        {
            combined_summary.success = false;
            append_failure(verification, "QMOM combined source plus boundary contribution exceeded the analytical tolerance.");
        }
    }
    verification.runs.push_back(combined_summary);

    struct BoundaryFailureCase
    {
        std::string name;
        qmomkit::MomentBoundaryContribution boundary;
        std::string expected_message_fragment;
    };
    const std::vector<BoundaryFailureCase> failure_cases{{"wrong size", [](double, const std::vector<double>& moments)
                                                          { return qmomkit::MomentBoundaryRate{std::vector<double>(moments.size() + 1, 0.0)}; },
                                                          "boundary contribution"},
                                                         {"NaN",
                                                          [](double, const std::vector<double>& moments)
                                                          {
                                                              std::vector<double> rate(moments.size(), 0.0);
                                                              rate.front() = std::numeric_limits<double>::quiet_NaN();
                                                              return qmomkit::MomentBoundaryRate{rate};
                                                          },
                                                          "boundary contribution"},
                                                         {"infinity",
                                                          [](double, const std::vector<double>& moments)
                                                          {
                                                              std::vector<double> rate(moments.size(), 0.0);
                                                              rate.front() = std::numeric_limits<double>::infinity();
                                                              return qmomkit::MomentBoundaryRate{rate};
                                                          },
                                                          "boundary contribution"}};

    for (const BoundaryFailureCase& failure_case : failure_cases)
    {
        const qmomkit::QmomTransportResult failure_result = qmomkit::advance_transport_qmom(initial, qmom_case.node_count, 0.0, 0.01, zero_velocity,
                                                                                            zero_qmom_source, failure_case.boundary, options);
        RunSummary failure_summary;
        failure_summary.scheme = options.scheme;
        failure_summary.dt = 0.01;
        failure_summary.steps = 1;
        failure_summary.final_time = failure_result.final_time;
        failure_summary.retry_count = failure_result.diagnostics.source_retries;
        failure_summary.success = !failure_result.ok() && failure_result.status == qmomkit::TransportStatus::source_failed &&
                                  failure_result.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step &&
                                  failure_result.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::source_failed &&
                                  contains_text(failure_result.message, failure_case.expected_message_fragment);
        if (!failure_summary.success)
        {
            append_failure(verification, "QMOM boundary failure case '" + failure_case.name + "' did not report source-stage boundary failure.");
        }
        verification.runs.push_back(failure_summary);
    }

    const qmomkit::MomentBoundaryContribution nonrealizable_boundary = [](double, const std::vector<double>& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        rate.front() = -100.0;
        return qmomkit::MomentBoundaryRate{rate};
    };
    const qmomkit::QmomTransportResult nonrealizable_result =
        qmomkit::advance_transport_qmom(initial, qmom_case.node_count, 0.0, 0.02, zero_velocity, zero_qmom_source, nonrealizable_boundary, options);
    RunSummary nonrealizable_summary;
    nonrealizable_summary.scheme = options.scheme;
    nonrealizable_summary.dt = 0.02;
    nonrealizable_summary.steps = 1;
    nonrealizable_summary.final_time = nonrealizable_result.final_time;
    nonrealizable_summary.retry_count = nonrealizable_result.diagnostics.source_retries;
    nonrealizable_summary.success = !nonrealizable_result.ok() && nonrealizable_result.status == qmomkit::TransportStatus::inversion_failed &&
                                    nonrealizable_result.diagnostics.failed_stage == qmomkit::TransportStage::post_source_inversion &&
                                    nonrealizable_result.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::success &&
                                    nonrealizable_result.diagnostics.post_source_qmom_status != qmomkit::QmomStatus::success;
    if (!nonrealizable_summary.success)
    {
        append_failure(verification, "QMOM boundary non-realizable update did not fail at post-source inversion.");
    }
    verification.runs.push_back(nonrealizable_summary);
    return verification;
}

VerificationResult run_noncommuting_qmom(const std::filesystem::path& configuration_path)
{
    // Verify the operator-order contract on a manufactured noncommuting QMOM
    // system. Strang should approach second-order convergence, while the fixed
    // CD and DC compositions should approach first order and differ at coarse
    // timestep.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_qmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load noncommuting QMOM configuration.");
        return verification;
    }

    const qmomkit::QmomCase& qmom_case = configured.value();
    const std::vector<double> initial = qmom_case.moments.values();
    const double t0 = 0.0;
    const double final_time = 1.2;
    const std::vector<double> exact = noncommuting_exact(initial, final_time);
    const std::vector<double> refinement_steps{0.3, 0.15, 0.075, 0.0375, 0.01875};

    const std::vector<qmomkit::TransportScheme> schemes{qmomkit::TransportScheme::strang, qmomkit::TransportScheme::cd, qmomkit::TransportScheme::dc};
    for (const qmomkit::TransportScheme scheme : schemes)
    {
        const auto options = qmom_options(qmom_case.options, 1.0e-12, 1.0e-14, 0.01, scheme);
        const std::size_t first_run_index = verification.runs.size();
        double previous_error = 0.0;
        for (const double dt : refinement_steps)
        {
            const std::size_t steps = static_cast<std::size_t>(std::lround(final_time / dt));
            const double tolerance = (scheme == qmomkit::TransportScheme::strang) ? 8.0e-4 : 3.0e-2;
            RunSummary summary = validate_qmom_run(verification, initial, qmom_case.node_count, t0, dt, steps, noncommuting_velocity,
                                                   noncommuting_source, options, exact, tolerance, tolerance);
            summary.output_row = true;
            if (previous_error > 0.0 && summary.error_norm > 0.0)
            {
                summary.observed_order = std::log(previous_error / summary.error_norm) / std::log(2.0);
            }
            previous_error = summary.error_norm;
            verification.runs.push_back(summary);
        }

        if (verification.runs.size() >= first_run_index + 5)
        {
            const double penultimate_order = verification.runs[first_run_index + 3].observed_order;
            const double finest_order = verification.runs[first_run_index + 4].observed_order;
            if (scheme == qmomkit::TransportScheme::strang)
            {
                if (penultimate_order < 1.75 || finest_order < 1.85)
                {
                    append_failure(verification, "noncommuting Strang refinement did not show asymptotic second-order convergence.");
                }
            }
            else if (penultimate_order < 0.75 || finest_order < 0.85)
            {
                append_failure(verification, "noncommuting first-order refinement did not show asymptotic first-order convergence.");
            }
            if (verification.runs[first_run_index + 4].errors.conservation > 1.0e-12)
            {
                append_failure(verification, "noncommuting run did not conserve m0.");
            }
        }
    }

    const auto cd_options = qmom_options(qmom_case.options, 1.0e-12, 1.0e-14, 0.01, qmomkit::TransportScheme::cd);
    const auto dc_options = qmom_options(qmom_case.options, 1.0e-12, 1.0e-14, 0.01, qmomkit::TransportScheme::dc);
    const QmomStepRun cd_coarse = run_qmom_steps(initial, qmom_case.node_count, t0, refinement_steps.front(),
                                                 static_cast<std::size_t>(std::lround(final_time / refinement_steps.front())), noncommuting_velocity,
                                                 noncommuting_source, cd_options);
    const QmomStepRun dc_coarse = run_qmom_steps(initial, qmom_case.node_count, t0, refinement_steps.front(),
                                                 static_cast<std::size_t>(std::lround(final_time / refinement_steps.front())), noncommuting_velocity,
                                                 noncommuting_source, dc_options);
    if (cd_coarse.result.ok() && dc_coarse.result.ok() && scaled_infinity_error(cd_coarse.result.moments, dc_coarse.result.moments) < 1.0e-4)
    {
        append_failure(verification, "noncommuting CD and DC schemes did not differ at finite timestep.");
    }
    return verification;
}

void verify_qmom_continuous_failure_result(const std::string& case_name, const qmomkit::QmomTransportResult& result,
                                           const std::vector<double>& initial_moments, const std::string& expected_message_fragment,
                                           VerificationResult& verification)
{
    // A velocity callback contract violation should stop in the continuous
    // node-advection stage. The returned moments/quadrature are diagnostic
    // state from before the failed substep, not a completed final state.
    if (result.ok())
    {
        append_failure(verification, case_name + " unexpectedly reported QMOM transport success.");
    }
    if (result.status != qmomkit::TransportStatus::continuous_failed)
    {
        append_failure(verification, case_name + " did not report TransportStatus::continuous_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::first_continuous_half_step)
    {
        append_failure(verification, case_name + " did not identify the first continuous half-step as the failed stage.");
    }
    if (result.diagnostics.first_continuous_status != qmomkit::NodeAdvectionStatus::velocity_failed)
    {
        append_failure(verification, case_name + " did not report NodeAdvectionStatus::velocity_failed.");
    }
    if (!contains_text(result.message, expected_message_fragment) || !contains_text(result.diagnostics.message, expected_message_fragment) ||
        !contains_text(result.diagnostics.first_continuous_diagnostics.message, expected_message_fragment))
    {
        append_failure(verification, case_name + " did not preserve a useful velocity-callback failure message.");
    }
    if (std::abs(result.final_time) > 1.0e-14)
    {
        append_failure(verification, case_name + " advanced final_time despite continuous transport failure.");
    }
    if (!vectors_equal(result.moments, initial_moments, 0.0))
    {
        append_failure(verification, case_name + " returned moments changed by a failed continuous substep.");
    }
    if (result.diagnostics.source_diagnostics.accepted || result.diagnostics.source_diagnostics.rejected ||
        result.diagnostics.first_reconstruction_diagnostics.success || result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, case_name + " reported later transport-stage diagnostics after the failed continuous substep.");
    }
}

void verify_qmom_rk23_minimum_step_failure_result(const qmomkit::QmomTransportResult& result, const std::vector<double>& initial_moments,
                                                  VerificationResult& verification)
{
    // The velocity callback is finite and dimensionally valid. Failure is
    // caused by adaptive RK23 rejecting the first trial step and being unable
    // to reduce below the configured minimum step.
    if (result.ok())
    {
        append_failure(verification, "QMOM RK23 minimum-step case unexpectedly reported transport success.");
    }
    if (result.status != qmomkit::TransportStatus::continuous_failed)
    {
        append_failure(verification, "QMOM RK23 minimum-step case did not report TransportStatus::continuous_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::first_continuous_half_step)
    {
        append_failure(verification, "QMOM RK23 minimum-step case did not identify the first continuous half-step as the failed stage.");
    }
    if (result.diagnostics.first_continuous_status != qmomkit::NodeAdvectionStatus::integration_failed ||
        result.diagnostics.first_continuous_diagnostics.adaptive_rk23_status != qmomkit::AdaptiveRk23Status::step_size_underflow)
    {
        append_failure(verification, "QMOM RK23 minimum-step case did not preserve the RK23 step-size-underflow status.");
    }
    if (result.diagnostics.first_continuous_diagnostics.rejected_step_count == 0 || result.diagnostics.total_node_advection_rejected_steps == 0 ||
        result.diagnostics.first_continuous_diagnostics.step_count != 0)
    {
        append_failure(verification, "QMOM RK23 minimum-step case did not preserve rejected-step diagnostics.");
    }
    if (!contains_text(result.message, "RK23 required a step smaller than min_step.") ||
        !contains_text(result.diagnostics.message, "RK23 required a step smaller than min_step."))
    {
        append_failure(verification, "QMOM RK23 minimum-step case did not preserve the minimum-step message.");
    }
    if (std::abs(result.final_time) > 1.0e-14 || std::abs(result.diagnostics.first_continuous_diagnostics.final_time) > 1.0e-14)
    {
        append_failure(verification, "QMOM RK23 minimum-step case advanced time despite rejecting all candidate steps.");
    }
    if (!vectors_equal(result.moments, initial_moments, 0.0) || result.diagnostics.source_diagnostics.accepted ||
        result.diagnostics.first_reconstruction_diagnostics.success || result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, "QMOM RK23 minimum-step case reported misleading successful output after failure.");
    }
}

void verify_qmom_post_source_inversion_failure_result(const qmomkit::QmomTransportResult& result, const double expected_failure_time,
                                                      VerificationResult& verification)
{
    // The source callback is finite, correctly sized, and successfully
    // integrated, but it drives m0 negative. The following QMOM inversion must
    // reject the evolved moments instead of returning a successful final
    // quadrature.
    if (result.ok())
    {
        append_failure(verification, "QMOM post-source inversion failure unexpectedly reported transport success.");
    }
    if (result.status != qmomkit::TransportStatus::inversion_failed)
    {
        append_failure(verification, "QMOM post-source inversion failure did not report TransportStatus::inversion_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::post_source_inversion)
    {
        append_failure(verification, "QMOM post-source inversion failure did not identify post_source_inversion as the failed stage.");
    }
    if (result.diagnostics.source_status != qmomkit::MomentSourceUpdateStatus::success ||
        result.diagnostics.source_diagnostics.ssp_rk3_status != qmomkit::SspRk3Status::success || !result.diagnostics.source_diagnostics.accepted)
    {
        append_failure(verification, "QMOM post-source inversion failure did not preserve the successful source-step diagnostics.");
    }
    if (result.diagnostics.post_source_qmom_status == qmomkit::QmomStatus::success || result.quadrature.ok())
    {
        append_failure(verification, "QMOM post-source inversion failure reported a successful post-source quadrature.");
    }
    if (result.message.empty() || result.diagnostics.message.empty() ||
        (!contains_text(result.message, "m0") && !contains_text(result.message, "realizable")))
    {
        append_failure(verification, "QMOM post-source inversion failure did not provide a meaningful reconstruction message.");
    }
    if (std::abs(result.final_time - expected_failure_time) > 1.0e-14)
    {
        append_failure(verification, "QMOM post-source inversion failure did not stop at the source-stage time.");
    }
    if (result.diagnostics.second_continuous_diagnostics.step_count > 0 || result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, "QMOM post-source inversion failure reported later-stage completion after the failed inversion.");
    }
}

void verify_qmom_initial_inversion_failure_result(const qmomkit::QmomTransportResult& result, const std::vector<double>& input_moments,
                                                  VerificationResult& verification)
{
    // This represents a user starting transport from a finite but
    // non-realizable QMOM moment sequence. Transport should fail before any
    // velocity or source callback can be used.
    if (result.ok())
    {
        append_failure(verification, "QMOM initial inversion failure unexpectedly reported transport success.");
    }
    if (result.status != qmomkit::TransportStatus::inversion_failed)
    {
        append_failure(verification, "QMOM initial inversion failure did not report TransportStatus::inversion_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::initial_inversion)
    {
        append_failure(verification, "QMOM initial inversion failure did not identify initial_inversion as the failed stage.");
    }
    if (result.diagnostics.initial_qmom_status == qmomkit::QmomStatus::success)
    {
        append_failure(verification, "QMOM initial inversion failure preserved a successful initial_qmom_status.");
    }
    if (result.message.empty() || result.diagnostics.message.empty() ||
        (!contains_text(result.message, "m0") && !contains_text(result.message, "realizable")))
    {
        append_failure(verification, "QMOM initial inversion failure did not provide a meaningful realizability message.");
    }
    if (std::abs(result.final_time) > 1.0e-14)
    {
        append_failure(verification, "QMOM initial inversion failure advanced final_time.");
    }
    if (!vectors_equal(result.moments, input_moments, 0.0))
    {
        append_failure(verification, "QMOM initial inversion failure changed the input moments.");
    }
    if (result.quadrature.ok() || result.diagnostics.source_diagnostics.accepted || result.diagnostics.first_reconstruction_diagnostics.success ||
        result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, "QMOM initial inversion failure reported later-stage or successful quadrature output.");
    }
}

VerificationResult run_qmom_transport_failures(const std::filesystem::path& configuration_path)
{
    // Exercise QMOM transport failure propagation after a valid initial
    // inversion: velocity callback errors, RK23 step-size underflow, and a
    // source step that creates a non-realizable post-source moment sequence.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_qmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load QMOM transport failure configuration.");
        return verification;
    }

    const qmomkit::QmomCase& qmom_case = configured.value();
    const std::vector<double> initial = qmom_case.moments.values();
    const auto options = qmom_options(qmom_case.options, 1.0e-12, 1.0e-14, 0.01);
    constexpr double dt = 0.08;

    struct FailureCase
    {
        std::string name;
        qmomkit::PhaseSpaceVelocity velocity;
        std::string expected_message_fragment;
    };

    const std::vector<FailureCase> failure_cases{{"wrong velocity dimension",
                                                  [](double, const std::vector<double>&)
                                                  {
                                                      // Simulates a user callback returning a vector that does not
                                                      // match the configured one-dimensional phase space.
                                                      return std::vector<double>{0.0, 0.0};
                                                  },
                                                  "Phase-space velocity callback returned the wrong dimension."},
                                                 {"non-finite velocity",
                                                  [](double, const std::vector<double>&)
                                                  {
                                                      // Simulates a physical model producing an invalid node velocity.
                                                      return std::vector<double>{std::numeric_limits<double>::quiet_NaN()};
                                                  },
                                                  "Phase-space velocity callback returned a non-finite value."}};

    for (const FailureCase& failure_case : failure_cases)
    {
        const qmomkit::QmomTransportResult result =
            qmomkit::advance_transport_qmom(initial, qmom_case.node_count, 0.0, dt, failure_case.velocity, zero_qmom_source, options);
        verify_qmom_continuous_failure_result(failure_case.name, result, initial, failure_case.expected_message_fragment, verification);

        RunSummary summary;
        summary.scheme = options.scheme;
        summary.dt = dt;
        summary.steps = 1;
        summary.final_time = result.final_time;
        summary.retry_count = result.diagnostics.source_retries;
        summary.success = !result.ok() && result.status == qmomkit::TransportStatus::continuous_failed &&
                          result.diagnostics.failed_stage == qmomkit::TransportStage::first_continuous_half_step &&
                          result.diagnostics.first_continuous_status == qmomkit::NodeAdvectionStatus::velocity_failed &&
                          vectors_equal(result.moments, initial, 0.0);
        verification.runs.push_back(summary);
    }

    qmomkit::QmomTransportOptions rk_failure_options = options;
    rk_failure_options.node_options.scheme = qmomkit::NodeAdvectionScheme::adaptive_rk23;
    rk_failure_options.node_options.rk23_options.initial_step = 0.04;
    rk_failure_options.node_options.rk23_options.max_step = 0.04;
    rk_failure_options.node_options.rk23_options.min_step = 0.02;
    const auto stiff_velocity = [](double, const std::vector<double>& xi)
    {
        // This is a valid velocity callback, but the deliberately stiff
        // exponential growth creates an oversized RK23 local-error estimate for
        // the first half-step. With min_step = dt/4, the rejected candidate
        // cannot be retried at a small enough step.
        return std::vector<double>{1000.0 * xi[0]};
    };
    const qmomkit::QmomTransportResult rk_failure =
        qmomkit::advance_transport_qmom(initial, qmom_case.node_count, 0.0, dt, stiff_velocity, zero_qmom_source, rk_failure_options);
    verify_qmom_rk23_minimum_step_failure_result(rk_failure, initial, verification);

    RunSummary rk_summary;
    rk_summary.scheme = rk_failure_options.scheme;
    rk_summary.dt = dt;
    rk_summary.steps = 1;
    rk_summary.final_time = rk_failure.final_time;
    rk_summary.retry_count = rk_failure.diagnostics.source_retries;
    rk_summary.success =
        !rk_failure.ok() && rk_failure.status == qmomkit::TransportStatus::continuous_failed &&
        rk_failure.diagnostics.failed_stage == qmomkit::TransportStage::first_continuous_half_step &&
        rk_failure.diagnostics.first_continuous_status == qmomkit::NodeAdvectionStatus::integration_failed &&
        rk_failure.diagnostics.first_continuous_diagnostics.adaptive_rk23_status == qmomkit::AdaptiveRk23Status::step_size_underflow &&
        rk_failure.diagnostics.first_continuous_diagnostics.rejected_step_count > 0 &&
        rk_failure.diagnostics.first_continuous_diagnostics.step_count == 0;
    verification.runs.push_back(rk_summary);

    const auto nonrealizable_source = [](double, const std::vector<double>& moments)
    {
        std::vector<double> rates(moments.size(), 0.0);
        if (!rates.empty())
        {
            rates.front() = -20.0;
        }
        return rates;
    };
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::QmomTransportResult post_source_inversion =
        qmomkit::advance_transport_qmom(initial, qmom_case.node_count, 0.0, dt, zero_velocity, nonrealizable_source, options);
    verify_qmom_post_source_inversion_failure_result(post_source_inversion, 0.5 * dt, verification);

    RunSummary reconstruction_summary;
    reconstruction_summary.scheme = options.scheme;
    reconstruction_summary.dt = dt;
    reconstruction_summary.steps = 1;
    reconstruction_summary.final_time = post_source_inversion.final_time;
    reconstruction_summary.retry_count = post_source_inversion.diagnostics.source_retries;
    reconstruction_summary.success = !post_source_inversion.ok() && post_source_inversion.status == qmomkit::TransportStatus::inversion_failed &&
                                     post_source_inversion.diagnostics.failed_stage == qmomkit::TransportStage::post_source_inversion &&
                                     post_source_inversion.diagnostics.post_source_qmom_status != qmomkit::QmomStatus::success;
    verification.runs.push_back(reconstruction_summary);
    return verification;
}

VerificationResult run_qmom_transport_inversion_failures(const std::filesystem::path& configuration_path)
{
    // Initial-inversion failure simulates a user starting transport from an
    // impossible moment sequence. No callback should run and the result should
    // preserve the input moments as diagnostic, non-success data.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_qmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load QMOM transport inversion-failure configuration.");
        return verification;
    }

    const qmomkit::QmomCase& qmom_case = configured.value();
    const std::vector<double> input = qmom_case.moments.values();
    const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(
        input, qmom_case.node_count, 0.0, 0.08, commuting_velocity, zero_qmom_source, qmom_options(qmom_case.options, 1.0e-12, 1.0e-14, 0.01));
    verify_qmom_initial_inversion_failure_result(result, input, verification);

    RunSummary summary;
    summary.scheme = qmomkit::TransportScheme::strang;
    summary.dt = 0.08;
    summary.steps = 1;
    summary.final_time = result.final_time;
    summary.retry_count = result.diagnostics.source_retries;
    summary.success = !result.ok() && result.status == qmomkit::TransportStatus::inversion_failed &&
                      result.diagnostics.failed_stage == qmomkit::TransportStage::initial_inversion &&
                      result.diagnostics.initial_qmom_status != qmomkit::QmomStatus::success;
    verification.runs.push_back(summary);
    return verification;
}

std::vector<double> cqmom_exact_values(const qmomkit::RawMomentTensor& initial, const double elapsed)
{
    // CQMOM commuting reference: diagonal phase-space scaling and diagonal
    // source damping give each tensor moment M_ij an exact scalar exponential
    // factor based on its multi-index.
    std::vector<double> exact(initial.values.size(), 0.0);
    for (std::size_t flat = 0; flat < initial.values.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, initial.extents);
        const double lambda = cqmom_lambda_scale * static_cast<double>(index[0] + index[1]);
        const double exponent = static_cast<double>(index[0]) * cqmom_a + static_cast<double>(index[1]) * cqmom_c - lambda;
        exact[flat] = initial.values[flat] * std::exp(exponent * elapsed);
    }
    return exact;
}

std::vector<double> cqmom_continuous_exact_values(const qmomkit::RawMomentTensor& initial, const double elapsed)
{
    std::vector<double> exact(initial.values.size(), 0.0);
    for (std::size_t flat = 0; flat < initial.values.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, initial.extents);
        const double exponent = static_cast<double>(index[0]) * cqmom_a + static_cast<double>(index[1]) * cqmom_c;
        exact[flat] = initial.values[flat] * std::exp(exponent * elapsed);
    }
    return exact;
}

std::vector<double> cqmom_velocity(double, const std::vector<double>& xi)
{
    return {cqmom_a * xi[0], cqmom_c * xi[1]};
}

qmomkit::CqmomMomentRate cqmom_source_for_context(double, const qmomkit::CqmomSourceContext& moments)
{
    std::vector<double> rhs(moments.size(), 0.0);
    for (std::size_t entry = 0; entry < moments.size(); ++entry)
    {
        const qmomkit::MultiIndex index = moments.index(entry);
        const double lambda = cqmom_lambda_scale * static_cast<double>(index[0] + index[1]);
        rhs[entry] = -lambda * moments.value(entry);
    }
    return qmomkit::CqmomMomentRate{std::move(rhs)};
}

qmomkit::CqmomMomentRate cqmom_state_zero_source(double, const qmomkit::CqmomSourceContext& moments)
{
    return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)};
}

std::vector<qmomkit::MultiIndex> cqmom_native_source_context_indices()
{
    // For a fixed 2x2 hierarchy with conditioning_order = {0, 1}, the
    // user-facing CQMOM source context exposes the required raw moments in this
    // phase-space multi-index order. The callback verifies this order directly
    // because returned CqmomMomentRate values are mapped by position.
    return {{0, 0}, {1, 0}, {2, 0}, {3, 0}, {0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}, {1, 3}};
}

std::vector<double> cqmom_native_source_terminal_x()
{
    return {1.0, 1.0, 3.0, 3.0};
}

std::vector<double> cqmom_native_source_terminal_y()
{
    return {-1.0, 2.0, -1.0, 2.0};
}

std::vector<double> cqmom_native_source_initial_weights()
{
    return {0.4 * 0.25, 0.4 * 0.75, 0.6 * 0.25, 0.6 * 0.75};
}

std::vector<double> cqmom_native_source_weight_rates()
{
    return {0.02, -0.01, 0.03, -0.015};
}

double cqmom_native_source_weight_rate_moment(const qmomkit::MultiIndex& index)
{
    double moment = 0.0;
    const std::vector<double> x = cqmom_native_source_terminal_x();
    const std::vector<double> y = cqmom_native_source_terminal_y();
    const std::vector<double> rates = cqmom_native_source_weight_rates();
    for (std::size_t node = 0; node < rates.size(); ++node)
    {
        moment += rates[node] * integer_power(x[node], index[0]) * integer_power(y[node], index[1]);
    }
    return moment;
}

std::vector<double> cqmom_native_source_required_values(const qmomkit::RawMomentTensor& moments)
{
    std::vector<double> values;
    for (const qmomkit::MultiIndex& index : cqmom_native_source_context_indices())
    {
        values.push_back(moments.values[qmomkit::flatten_multi_index(index, moments.extents)]);
    }
    return values;
}

double ssp_rk3_linear_affine_factor(const double z)
{
    return 1.0 + z + 0.5 * z * z + (z * z * z) / 6.0;
}

double ssp_rk3_constant_source_factor(const double z)
{
    return 1.0 + 0.5 * z + (z * z) / 6.0;
}

std::vector<double> cqmom_native_source_exact_values(const qmomkit::RawMomentTensor& initial, const double time_step)
{
    // Native CQMOM source reference for the manufactured affine source
    // dM/dt = gamma M + constant. The source integrator uses a single SSP-RK3
    // step, so the expected value uses the method's stability polynomial rather
    // than the exact exponential solution.
    std::vector<double> exact(initial.values.size(), 0.0);
    const double z = cqmom_native_source_gamma * time_step;
    const double linear_factor = ssp_rk3_linear_affine_factor(z);
    const double source_factor = time_step * ssp_rk3_constant_source_factor(z);
    for (std::size_t flat = 0; flat < initial.values.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, initial.extents);
        exact[flat] = linear_factor * initial.values[flat] + source_factor * cqmom_native_source_weight_rate_moment(index);
    }
    return exact;
}

std::vector<double> cqmom_combined_expected_source_context_values(const qmomkit::RawMomentTensor& initial, const double half_step)
{
    // In Strang mode the native CQMOM source sees the state after the first
    // continuous half-step. For dx/dt = a x and dy/dt = c y, raw moment M_ij
    // is scaled by exp((i a + j c) half_step).
    std::vector<double> expected;
    for (const qmomkit::MultiIndex& index : cqmom_native_source_context_indices())
    {
        const double exponent = static_cast<double>(index[0]) * cqmom_a + static_cast<double>(index[1]) * cqmom_c;
        const double initial_value = initial.values[qmomkit::flatten_multi_index(index, initial.extents)];
        expected.push_back(initial_value * std::exp(exponent * half_step));
    }
    return expected;
}

std::vector<double> cqmom_combined_native_source_exact_values(const qmomkit::RawMomentTensor& initial, const double time_step)
{
    // The manufactured source is dM_alpha/dt = gamma M_alpha. The transport
    // framework advances that source with one SSP-RK3 step, whose stability
    // polynomial is P(z) = 1 + z + z^2/2 + z^3/6. The split reference is
    // therefore continuous(dt/2), source(dt), continuous(dt/2).
    std::vector<double> exact(initial.values.size(), 0.0);
    const double source_factor = ssp_rk3_linear_affine_factor(cqmom_combined_source_gamma * time_step);
    for (std::size_t flat = 0; flat < initial.values.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, initial.extents);
        const double continuous_exponent = static_cast<double>(index[0]) * cqmom_a + static_cast<double>(index[1]) * cqmom_c;
        exact[flat] = initial.values[flat] * std::exp(continuous_exponent * time_step) * source_factor;
    }
    return exact;
}

qmomkit::CqmomMomentRate cqmom_native_source_callback(const std::vector<double>& expected_initial_values, CqmomNativeSourceAudit& audit, double,
                                                      const qmomkit::CqmomSourceContext& context)
{
    // The native source callback doubles as a public API audit. It checks the
    // context hierarchy, required-index order, and initial values before
    // returning a manufactured affine source in that same public order.
    ++audit.invocation_count;
    const std::vector<qmomkit::MultiIndex> expected_indices = cqmom_native_source_context_indices();
    if (context.size() != expected_indices.size())
    {
        record_audit_failure(audit, "CQMOM source context size did not match the expected transported moment count.");
    }
    if (context.hierarchy().node_counts() != std::vector<std::size_t>{2, 2})
    {
        record_audit_failure(audit, "CQMOM source context hierarchy node_counts were not {2, 2}.");
    }
    if (context.hierarchy().conditioning_order() != std::vector<std::size_t>{0, 1})
    {
        record_audit_failure(audit, "CQMOM source context conditioning_order was not {0, 1}.");
    }

    const std::size_t entries = std::min(context.size(), expected_indices.size());
    std::vector<double> rates(context.size(), 0.0);
    for (std::size_t entry = 0; entry < entries; ++entry)
    {
        const qmomkit::MultiIndex& index = context.index(entry);
        if (index != expected_indices[entry])
        {
            record_audit_failure(audit, "CQMOM source context index order did not match the expected public moment order.");
            continue;
        }
        if (audit.invocation_count == 1 && entry < expected_initial_values.size() &&
            std::abs(context.value(entry) - expected_initial_values[entry]) > 1.0e-12)
        {
            record_audit_failure(audit, "CQMOM source context initial value did not match the corresponding transported moment.");
        }

        rates[entry] = cqmom_native_source_gamma * context.value(entry) + cqmom_native_source_weight_rate_moment(index);
    }
    if (audit.invocation_count == 1)
    {
        audit.checked_initial_context_values = true;
    }
    return qmomkit::CqmomMomentRate{rates};
}

qmomkit::CqmomMomentRate cqmom_combined_native_source_callback(const std::vector<double>& expected_source_values, CqmomNativeSourceAudit& audit,
                                                               double, const qmomkit::CqmomSourceContext& context)
{
    // In the combined Strang case, the callback should see moments after the
    // first continuous half-step. Verifying those context values protects the
    // public timing of native CQMOM physics callbacks.
    ++audit.invocation_count;
    const std::vector<qmomkit::MultiIndex> expected_indices = cqmom_native_source_context_indices();
    if (context.size() != expected_indices.size())
    {
        record_audit_failure(audit, "combined CQMOM source context size did not match the expected transported moment count.");
    }
    if (context.hierarchy().node_counts() != std::vector<std::size_t>{2, 2})
    {
        record_audit_failure(audit, "combined CQMOM source context hierarchy node_counts were not {2, 2}.");
    }
    if (context.hierarchy().conditioning_order() != std::vector<std::size_t>{0, 1})
    {
        record_audit_failure(audit, "combined CQMOM source context conditioning_order was not {0, 1}.");
    }

    const std::size_t entries = std::min(context.size(), expected_indices.size());
    std::vector<double> rates(context.size(), 0.0);
    for (std::size_t entry = 0; entry < entries; ++entry)
    {
        const qmomkit::MultiIndex& index = context.index(entry);
        if (index != expected_indices[entry])
        {
            record_audit_failure(audit, "combined CQMOM source context index order did not match the expected public moment order.");
            continue;
        }
        if (audit.invocation_count == 1 && entry < expected_source_values.size() &&
            std::abs(context.value(entry) - expected_source_values[entry]) > 5.0e-11)
        {
            record_audit_failure(audit, "combined CQMOM source context value did not match the post-continuous-half-step moment.");
        }

        rates[entry] = cqmom_combined_source_gamma * context.value(entry);
    }
    if (audit.invocation_count == 1)
    {
        audit.checked_split_context_values = true;
    }
    return qmomkit::CqmomMomentRate{rates};
}

std::vector<qmomkit::MultiIndex> cqmom_boundary_context_indices()
{
    return {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {1, 0}, {2, 0}, {3, 0}, {1, 1}, {2, 1}, {3, 1}};
}

std::vector<double> cqmom_boundary_weight_rates(const qmomkit::CqmomQuadrature& quadrature)
{
    std::vector<double> rates(quadrature.weights().size(), 0.0);
    for (std::size_t node = 0; node < rates.size(); ++node)
    {
        rates[node] = 0.01 * static_cast<double>(node + 1);
    }
    return rates;
}

double cqmom_boundary_weight_rate_moment(const qmomkit::CqmomQuadrature& quadrature, const std::vector<double>& rates,
                                         const qmomkit::MultiIndex& index)
{
    double moment_rate = 0.0;
    for (std::size_t node = 0; node < rates.size(); ++node)
    {
        double contribution = rates[node];
        for (std::size_t dimension = 0; dimension < quadrature.abscissas().size(); ++dimension)
        {
            contribution *= integer_power(quadrature.abscissas()[dimension][node], index[dimension]);
        }
        moment_rate += contribution;
    }
    return moment_rate;
}

qmomkit::CqmomMomentRate cqmom_boundary_contribution_callback(const qmomkit::CqmomQuadrature& initial_quadrature, const std::vector<double>& rates,
                                                              CqmomNativeSourceAudit& audit, double, const qmomkit::CqmomSourceContext& context)
{
    ++audit.invocation_count;
    const std::vector<qmomkit::MultiIndex> expected_indices = cqmom_boundary_context_indices();
    if (context.size() != expected_indices.size())
    {
        record_audit_failure(audit, "CQMOM boundary context size did not match the expected transported moment count.");
    }
    if (context.hierarchy().node_counts() != std::vector<std::size_t>{2, 2})
    {
        record_audit_failure(audit, "CQMOM boundary context hierarchy node_counts were not {2, 2}.");
    }
    if (context.hierarchy().conditioning_order() != std::vector<std::size_t>{1, 0})
    {
        record_audit_failure(audit, "CQMOM boundary context conditioning_order was not {1, 0}.");
    }

    const std::size_t entries = std::min(context.size(), expected_indices.size());
    std::vector<double> boundary_rates(context.size(), 0.0);
    for (std::size_t entry = 0; entry < entries; ++entry)
    {
        const qmomkit::MultiIndex& index = context.index(entry);
        if (index != expected_indices[entry])
        {
            record_audit_failure(audit, "CQMOM boundary context index order did not match the expected public moment order.");
            continue;
        }
        boundary_rates[entry] = cqmom_boundary_weight_rate_moment(initial_quadrature, rates, index);
    }
    if (audit.invocation_count == 1)
    {
        audit.checked_initial_context_values = true;
    }
    return qmomkit::CqmomMomentRate{boundary_rates};
}

VerificationResult run_commuting_cqmom_2d(const std::filesystem::path& configuration_path)
{
    // Verify dense CQMOM transport on a two-dimensional commuting analytical
    // problem. CQMOM source callbacks act on the transport moment layout, so the
    // final transported moments and quadrature reconstruction are compared with
    // the exact tensor exponential reference on that layout.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom_dense(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load commuting CQMOM configuration.");
        return verification;
    }

    const qmomkit::CqmomDenseCase& dense_case = configured.value();
    const qmomkit::CqmomMomentSource source = cqmom_source_for_context;
    const auto options = cqmom_options(dense_case);
    const double dt = 0.04;
    const std::size_t steps = 5;
    const double elapsed = dt * static_cast<double>(steps);
    qmomkit::RawMomentTensor moments = dense_case.moments;
    double time = 0.0;
    bool diagnostics_successful = true;
    std::size_t retry_count = 0;
    qmomkit::TestDenseCqmomTransportResult result;
    for (std::size_t step = 0; step < steps; ++step)
    {
        result = qmomkit::advance_transport_cqmom(moments, dense_case.hierarchy, time, dt, cqmom_velocity, source, options);
        if (!result.ok())
        {
            append_failure(verification, "CQMOM transport integrator failed: " + result.message);
            return verification;
        }
        if (!cqmom_diagnostics_are_successful(result.diagnostics))
        {
            append_failure(verification, "CQMOM transport diagnostics did not report successful completion of every substage.");
            diagnostics_successful = false;
        }
        retry_count += result.diagnostics.source_retries;
        moments = result.moments;
        time = result.final_time;
    }

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = dt;
    summary.steps = steps;
    summary.final_time = time;
    summary.retry_count = retry_count;
    summary.success = result.ok() && diagnostics_successful;
    summary.output_row = true;
    const std::vector<double> exact = cqmom_exact_values(dense_case.moments, elapsed);
    const qmomkit::CqmomMomentLayout transport_layout = qmomkit::create_cqmom_moment_layout(dense_case.hierarchy);
    const qmomkit::CqmomMomentData result_required = qmomkit::extract_cqmom_moments(result.moments, transport_layout);
    const qmomkit::CqmomMomentData exact_required =
        qmomkit::extract_cqmom_moments(qmomkit::RawMomentTensor{dense_case.moments.extents, exact}, transport_layout);
    update_errors(summary.errors, result_required.values, exact_required.values);
    summary.error_norm = scaled_infinity_error(result_required.values, exact_required.values);
    if (result.quadrature.has_value())
    {
        const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(*result.quadrature, dense_case.moments.extents);
        const qmomkit::CqmomMomentData reconstructed_required = qmomkit::extract_cqmom_moments(reconstructed, transport_layout);
        update_reconstruction_errors(summary.errors, reconstructed_required.values, exact_required.values);
    }
    else
    {
        summary.success = false;
        append_failure(verification, "CQMOM transport result did not retain final quadrature.");
    }
    print_vector_comparison("Final CQMOM transport-layout moments", result_required.values, exact_required.values);
    if (std::abs(time - elapsed) > 1.0e-12)
    {
        summary.success = false;
        append_failure(verification, "CQMOM final time does not match the requested final time.");
    }
    if (summary.errors.moment_scaled > 3.0e-6 || summary.errors.reconstruction_scaled > 3.0e-6)
    {
        summary.success = false;
        append_failure(verification, "CQMOM commuting transport exceeded the analytical tolerance.");
    }
    verification.runs.push_back(summary);
    return verification;
}

VerificationResult run_cqmom_boundary_contribution(const std::filesystem::path& configuration_path)
{
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load CQMOM boundary contribution configuration.");
        return verification;
    }

    const qmomkit::CqmomCase& cqmom_case = configured.value();
    if (cqmom_case.hierarchy().node_counts() != std::vector<std::size_t>{2, 2} ||
        cqmom_case.hierarchy().conditioning_order() != std::vector<std::size_t>{1, 0})
    {
        append_failure(verification, "CQMOM boundary verification requires a fixed 2x2 hierarchy with conditioning_order = {1, 0}.");
        return verification;
    }

    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_case);
    if (!initial.ok())
    {
        append_failure(verification, "initial CQMOM boundary contribution inversion failed: " + initial.message);
        return verification;
    }

    const std::vector<std::size_t> extents{4, 4};
    const qmomkit::RawMomentTensor initial_moments = qmomkit::reconstruct_raw_moments(initial.quadrature, extents);
    const std::vector<double> weight_rates = cqmom_boundary_weight_rates(initial.quadrature);
    CqmomNativeSourceAudit audit;
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto boundary = [&initial, &weight_rates, &audit](const double time, const qmomkit::CqmomSourceContext& context)
    { return cqmom_boundary_contribution_callback(initial.quadrature, weight_rates, audit, time, context); };

    qmomkit::CqmomState state = cqmom_case.initial_state();
    qmomkit::CqmomTransportOptions options = cqmom_state_options();
    options.scheme = qmomkit::TransportScheme::strang;
    constexpr double dt = 0.02;
    const qmomkit::CqmomTransportResult result =
        qmomkit::advance_transport_cqmom(state, 0.0, dt, zero_velocity, cqmom_state_zero_source, boundary, options);

    if (!result.ok())
    {
        append_failure(verification, "CQMOM boundary contribution transport failed: " + result.message);
        return verification;
    }
    if (!audit.passed)
    {
        append_failure(verification, audit.message);
    }
    if (audit.invocation_count == 0 || !audit.checked_initial_context_values)
    {
        append_failure(verification, "CQMOM boundary contribution callback did not inspect the source context.");
    }
    if (!cqmom_diagnostics_are_successful(result.diagnostics))
    {
        append_failure(verification, "CQMOM boundary contribution diagnostics did not report successful completion of every substage.");
    }

    std::vector<double> exact = initial_moments.values;
    for (std::size_t flat = 0; flat < exact.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, extents);
        exact[flat] += dt * cqmom_boundary_weight_rate_moment(initial.quadrature, weight_rates, index);
    }

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = dt;
    summary.steps = 1;
    summary.final_time = result.final_time;
    summary.retry_count = result.diagnostics.source_retries;
    summary.success = result.ok() && audit.passed && audit.invocation_count > 0 && cqmom_diagnostics_are_successful(result.diagnostics);
    summary.output_row = true;
    if (result.quadrature.has_value())
    {
        const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(*result.quadrature, extents);
        update_errors(summary.errors, reconstructed.values, exact);
        update_reconstruction_errors(summary.errors, reconstructed.values, exact);
        summary.error_norm = scaled_infinity_error(reconstructed.values, exact);
        print_vector_comparison("Final CQMOM boundary contribution moments", reconstructed.values, exact);
    }
    else
    {
        summary.success = false;
        append_failure(verification, "CQMOM boundary contribution result did not retain final quadrature.");
    }
    if (std::abs(result.final_time - dt) > 1.0e-12)
    {
        summary.success = false;
        append_failure(verification, "CQMOM boundary contribution final time does not match the requested final time.");
    }
    if (summary.errors.moment_scaled > 5.0e-9 || summary.errors.reconstruction_scaled > 5.0e-9)
    {
        summary.success = false;
        append_failure(verification, "CQMOM boundary contribution moments exceeded the analytical tolerance.");
    }
    verification.runs.push_back(summary);
    return verification;
}

VerificationResult run_cqmom_combined_native_source(const std::filesystem::path& configuration_path)
{
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load combined CQMOM transport/source configuration.");
        return verification;
    }

    const qmomkit::CqmomCase& cqmom_case = configured.value();
    if (cqmom_case.hierarchy().node_counts() != std::vector<std::size_t>{2, 2} ||
        cqmom_case.hierarchy().conditioning_order() != std::vector<std::size_t>{0, 1})
    {
        append_failure(verification, "combined CQMOM verification requires a fixed 2x2 hierarchy with conditioning_order = {0, 1}.");
        return verification;
    }

    // This case exercises the public workflow used by coupled CQMOM physics
    // models: a configured state is continuously advected, exposed to a native
    // CqmomMomentSource callback, then reconstructed again after the final
    // continuous half-step.
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_case);
    if (!initial.ok())
    {
        append_failure(verification, "initial combined CQMOM inversion failed: " + initial.message);
        return verification;
    }

    constexpr double dt = 0.08;
    const std::vector<std::size_t> extents{4, 4};
    const qmomkit::RawMomentTensor initial_moments = qmomkit::reconstruct_raw_moments(initial.quadrature, extents);
    const std::vector<double> expected_source_values = cqmom_combined_expected_source_context_values(initial_moments, 0.5 * dt);
    CqmomNativeSourceAudit audit;
    const auto source = [&expected_source_values, &audit](const double time, const qmomkit::CqmomSourceContext& context)
    { return cqmom_combined_native_source_callback(expected_source_values, audit, time, context); };

    qmomkit::CqmomState state = cqmom_case.initial_state();
    qmomkit::CqmomTransportOptions options = cqmom_state_options();
    options.scheme = qmomkit::TransportScheme::strang;

    const qmomkit::CqmomTransportResult result = qmomkit::advance_transport_cqmom(state, 0.0, dt, cqmom_velocity, source, options);
    if (!result.ok())
    {
        append_failure(verification, "combined CQMOM transport/source workflow failed: " + result.message);
        return verification;
    }
    if (!audit.passed)
    {
        append_failure(verification, audit.message);
    }
    if (audit.invocation_count == 0)
    {
        append_failure(verification, "combined CQMOM source callback was not invoked.");
    }
    if (!audit.checked_split_context_values)
    {
        append_failure(verification, "combined CQMOM source callback did not inspect the post-continuous-half-step context values.");
    }
    if (!cqmom_diagnostics_are_successful(result.diagnostics))
    {
        append_failure(verification, "combined CQMOM transport diagnostics did not report successful completion of every substage.");
    }
    if (result.diagnostics.total_node_advection_steps == 0)
    {
        append_failure(verification, "combined CQMOM continuous transport did not report any node-advection steps.");
    }
    if (result.diagnostics.source_diagnostics.ssp_rk3_status != qmomkit::SspRk3Status::success ||
        std::abs(result.diagnostics.source_diagnostics.time_step_used - dt) > 1.0e-14)
    {
        append_failure(verification, "combined CQMOM source diagnostics did not report the expected successful SSP-RK3 source step.");
    }

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = dt;
    summary.steps = 1;
    summary.final_time = result.final_time;
    summary.retry_count = result.diagnostics.source_retries;
    summary.success = result.ok() && audit.passed && audit.invocation_count > 0 && audit.checked_split_context_values &&
                      cqmom_diagnostics_are_successful(result.diagnostics) && result.diagnostics.total_node_advection_steps > 0 &&
                      result.diagnostics.source_diagnostics.ssp_rk3_status == qmomkit::SspRk3Status::success &&
                      std::abs(result.diagnostics.source_diagnostics.time_step_used - dt) <= 1.0e-14;
    summary.output_row = true;
    const std::vector<double> exact = cqmom_combined_native_source_exact_values(initial_moments, dt);
    if (result.quadrature.has_value())
    {
        const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(*result.quadrature, extents);
        update_errors(summary.errors, reconstructed.values, exact);
        update_reconstruction_errors(summary.errors, reconstructed.values, exact);
        summary.error_norm = scaled_infinity_error(reconstructed.values, exact);
        summary.errors.conservation = std::abs(reconstructed.values.front() - exact.front());
        print_vector_comparison("Final combined CQMOM transport/source moments", reconstructed.values, exact);
    }
    else
    {
        summary.success = false;
        append_failure(verification, "combined CQMOM transport result did not retain final quadrature.");
    }
    if (std::abs(result.final_time - dt) > 1.0e-12)
    {
        summary.success = false;
        append_failure(verification, "combined CQMOM final time does not match the requested final time.");
    }
    if (summary.errors.moment_scaled > 5.0e-9 || summary.errors.reconstruction_scaled > 5.0e-9)
    {
        summary.success = false;
        append_failure(verification, "combined CQMOM transport/source moments exceeded the analytical split-reference tolerance.");
    }
    verification.runs.push_back(summary);
    return verification;
}

VerificationResult run_cqmom_native_source(const std::filesystem::path& configuration_path)
{
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load CQMOM-native source configuration.");
        return verification;
    }

    const qmomkit::CqmomCase& cqmom_case = configured.value();
    if (cqmom_case.hierarchy().node_counts() != std::vector<std::size_t>{2, 2} ||
        cqmom_case.hierarchy().conditioning_order() != std::vector<std::size_t>{0, 1})
    {
        append_failure(verification, "CQMOM-native source verification requires a fixed 2x2 hierarchy with conditioning_order = {0, 1}.");
        return verification;
    }

    // This is the normal user workflow: load_cqmom(...) owns the configured
    // state, compute_cqmom(...) verifies the initial analytical quadrature, and
    // CqmomCase::initial_state() is passed to advance_transport_cqmom(...).
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_case);
    if (!initial.ok())
    {
        append_failure(verification, "initial CQMOM-native source inversion failed: " + initial.message);
        return verification;
    }

    const std::vector<std::size_t> extents{4, 4};
    const qmomkit::RawMomentTensor initial_moments = qmomkit::reconstruct_raw_moments(initial.quadrature, extents);
    const std::vector<double> expected_initial_values = cqmom_native_source_required_values(initial_moments);
    CqmomNativeSourceAudit audit;
    const auto source = [&expected_initial_values, &audit](const double time, const qmomkit::CqmomSourceContext& context)
    { return cqmom_native_source_callback(expected_initial_values, audit, time, context); };
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    qmomkit::CqmomState state = cqmom_case.initial_state();
    qmomkit::CqmomTransportOptions options = cqmom_state_options();
    options.scheme = qmomkit::TransportScheme::strang;

    constexpr double dt = 0.08;
    const qmomkit::CqmomTransportResult result = qmomkit::advance_transport_cqmom(state, 0.0, dt, zero_velocity, source, options);
    if (!result.ok())
    {
        append_failure(verification, "CQMOM-native source transport failed: " + result.message);
        return verification;
    }
    if (!audit.passed)
    {
        append_failure(verification, audit.message);
    }
    if (audit.invocation_count == 0)
    {
        append_failure(verification, "CQMOM-native source callback was not invoked.");
    }
    if (!audit.checked_initial_context_values)
    {
        append_failure(verification, "CQMOM-native source callback did not inspect the initial context values.");
    }
    if (!cqmom_diagnostics_are_successful(result.diagnostics))
    {
        append_failure(verification, "CQMOM-native source transport diagnostics did not report successful completion of every substage.");
    }
    const bool diagnostics_successful = cqmom_diagnostics_are_successful(result.diagnostics);

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = dt;
    summary.steps = 1;
    summary.final_time = result.final_time;
    summary.retry_count = result.diagnostics.source_retries;
    summary.success = result.ok() && audit.passed && audit.invocation_count > 0 && audit.checked_initial_context_values && diagnostics_successful;
    summary.output_row = true;
    const std::vector<double> exact = cqmom_native_source_exact_values(initial_moments, dt);
    if (result.quadrature.has_value())
    {
        const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(*result.quadrature, extents);
        update_errors(summary.errors, reconstructed.values, exact);
        update_reconstruction_errors(summary.errors, reconstructed.values, exact);
        summary.error_norm = scaled_infinity_error(reconstructed.values, exact);
        summary.errors.conservation = std::abs(reconstructed.values.front() - exact.front());
        print_vector_comparison("Final CQMOM-native source moments", reconstructed.values, exact);
    }
    else
    {
        summary.success = false;
        append_failure(verification, "CQMOM-native source transport result did not retain final quadrature.");
    }
    if (std::abs(result.final_time - dt) > 1.0e-12)
    {
        summary.success = false;
        append_failure(verification, "CQMOM-native source final time does not match the requested final time.");
    }
    if (summary.errors.moment_scaled > 5.0e-10 || summary.errors.reconstruction_scaled > 5.0e-10)
    {
        summary.success = false;
        append_failure(verification, "CQMOM-native source moments exceeded the analytical SSP-RK3 tolerance.");
    }
    verification.runs.push_back(summary);
    return verification;
}

void verify_cqmom_continuous_failure_result(const std::string& case_name, const qmomkit::CqmomTransportResult& result,
                                            const std::string& expected_message_fragment, VerificationResult& verification)
{
    // This is the CQMOM state counterpart of the QMOM velocity-callback
    // failure: the public result should identify continuous transport as the
    // failing operation and should not report source or reconstruction success
    // for stages that were never reached.
    if (result.ok())
    {
        append_failure(verification, case_name + " unexpectedly reported CQMOM transport success.");
    }
    if (result.status != qmomkit::TransportStatus::continuous_failed)
    {
        append_failure(verification, case_name + " did not report TransportStatus::continuous_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::first_continuous_half_step)
    {
        append_failure(verification, case_name + " did not identify the first continuous half-step as the failed stage.");
    }
    if (result.diagnostics.first_continuous_status != qmomkit::NodeAdvectionStatus::velocity_failed)
    {
        append_failure(verification, case_name + " did not report NodeAdvectionStatus::velocity_failed.");
    }
    if (!contains_text(result.message, expected_message_fragment) || !contains_text(result.diagnostics.message, expected_message_fragment) ||
        !contains_text(result.diagnostics.first_continuous_diagnostics.message, expected_message_fragment))
    {
        append_failure(verification, case_name + " did not preserve a useful velocity-callback failure message.");
    }
    if (std::abs(result.final_time) > 1.0e-14)
    {
        append_failure(verification, case_name + " advanced final_time despite continuous transport failure.");
    }
    if (result.diagnostics.source_diagnostics.accepted || result.diagnostics.source_diagnostics.rejected ||
        result.diagnostics.first_reconstruction_diagnostics.success || result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, case_name + " reported later transport-stage diagnostics after the failed continuous substep.");
    }
}

void verify_cqmom_rk23_minimum_step_failure_result(const qmomkit::CqmomTransportResult& result, VerificationResult& verification)
{
    // This mirrors the QMOM RK23 minimum-step case through the public CqmomState
    // workflow. The user velocity is valid; the adaptive node integrator is the
    // component that cannot complete the requested continuous half-step.
    if (result.ok())
    {
        append_failure(verification, "CQMOM RK23 minimum-step case unexpectedly reported transport success.");
    }
    if (result.status != qmomkit::TransportStatus::continuous_failed)
    {
        append_failure(verification, "CQMOM RK23 minimum-step case did not report TransportStatus::continuous_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::first_continuous_half_step)
    {
        append_failure(verification, "CQMOM RK23 minimum-step case did not identify the first continuous half-step as the failed stage.");
    }
    if (result.diagnostics.first_continuous_status != qmomkit::NodeAdvectionStatus::integration_failed ||
        result.diagnostics.first_continuous_diagnostics.adaptive_rk23_status != qmomkit::AdaptiveRk23Status::step_size_underflow)
    {
        append_failure(verification, "CQMOM RK23 minimum-step case did not preserve the RK23 step-size-underflow status.");
    }
    if (result.diagnostics.first_continuous_diagnostics.rejected_step_count == 0 || result.diagnostics.total_node_advection_rejected_steps == 0 ||
        result.diagnostics.first_continuous_diagnostics.step_count != 0)
    {
        append_failure(verification, "CQMOM RK23 minimum-step case did not preserve rejected-step diagnostics.");
    }
    if (!contains_text(result.message, "RK23 required a step smaller than min_step.") ||
        !contains_text(result.diagnostics.message, "RK23 required a step smaller than min_step."))
    {
        append_failure(verification, "CQMOM RK23 minimum-step case did not preserve the minimum-step message.");
    }
    if (std::abs(result.final_time) > 1.0e-14 || std::abs(result.diagnostics.first_continuous_diagnostics.final_time) > 1.0e-14)
    {
        append_failure(verification, "CQMOM RK23 minimum-step case advanced time despite rejecting all candidate steps.");
    }
    if (result.diagnostics.source_diagnostics.accepted || result.diagnostics.first_reconstruction_diagnostics.success ||
        result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, "CQMOM RK23 minimum-step case reported misleading successful output after failure.");
    }
}

void verify_cqmom_post_source_inversion_failure_result(const qmomkit::CqmomTransportResult& result, const double expected_failure_time,
                                                       VerificationResult& verification)
{
    // The source callback is finite and correctly sized, so the source update
    // itself succeeds. It deliberately drives M[0,0] negative, which is not a
    // realizable CQMOM moment set, and the following public CQMOM inversion
    // must stop the transport workflow before the second continuous half-step.
    if (result.ok())
    {
        append_failure(verification, "post-source inversion failure unexpectedly reported CQMOM transport success.");
    }
    if (result.status != qmomkit::TransportStatus::inversion_failed)
    {
        append_failure(verification, "post-source inversion failure did not report TransportStatus::inversion_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::post_source_inversion)
    {
        append_failure(verification, "post-source inversion failure did not identify post_source_inversion as the failed stage.");
    }
    if (result.diagnostics.source_status != qmomkit::MomentSourceUpdateStatus::success ||
        result.diagnostics.source_diagnostics.ssp_rk3_status != qmomkit::SspRk3Status::success || !result.diagnostics.source_diagnostics.accepted)
    {
        append_failure(verification, "post-source inversion failure did not preserve the successful source-step diagnostics.");
    }
    if (result.diagnostics.post_source_cqmom_status == qmomkit::CqmomStatus::success)
    {
        append_failure(verification, "post-source inversion failure reported a successful CQMOM post-source status.");
    }
    if (result.diagnostics.post_source_cqmom_diagnostics.empty())
    {
        append_failure(verification, "post-source inversion failure did not provide CQMOM diagnostics.");
    }
    if (result.message.empty() || result.diagnostics.message.empty() ||
        (!contains_text(result.message, "m0") && !contains_text(result.message, "realizable")))
    {
        append_failure(verification, "post-source inversion failure did not provide a meaningful reconstruction message.");
    }
    if (std::abs(result.final_time - expected_failure_time) > 1.0e-14)
    {
        append_failure(verification, "post-source inversion failure did not stop at the source-stage time.");
    }
    if (result.diagnostics.second_continuous_status != qmomkit::NodeAdvectionStatus::success ||
        result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, "post-source inversion failure reported later-stage completion after the failed inversion.");
    }
}

void verify_cqmom_initial_inversion_failure_result(const qmomkit::CqmomTransportResult& result, VerificationResult& verification)
{
    // This represents a configured CQMOM state whose selected moments are
    // finite but not realizable. Transport should fail during the initial
    // public CQMOM reconstruction, before continuous advection or source
    // integration can produce misleading progress diagnostics.
    if (result.ok())
    {
        append_failure(verification, "CQMOM initial inversion failure unexpectedly reported transport success.");
    }
    if (result.status != qmomkit::TransportStatus::inversion_failed)
    {
        append_failure(verification, "CQMOM initial inversion failure did not report TransportStatus::inversion_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::initial_inversion)
    {
        append_failure(verification, "CQMOM initial inversion failure did not identify initial_inversion as the failed stage.");
    }
    if (result.diagnostics.initial_cqmom_status == qmomkit::CqmomStatus::success)
    {
        append_failure(verification, "CQMOM initial inversion failure preserved a successful initial_cqmom_status.");
    }
    if (result.diagnostics.initial_cqmom_diagnostics.empty())
    {
        append_failure(verification, "CQMOM initial inversion failure did not provide CQMOM diagnostics.");
    }
    if (result.message.empty() || result.diagnostics.message.empty() ||
        (!contains_text(result.message, "m0") && !contains_text(result.message, "realizable")))
    {
        append_failure(verification, "CQMOM initial inversion failure did not provide a meaningful realizability message.");
    }
    if (std::abs(result.final_time) > 1.0e-14)
    {
        append_failure(verification, "CQMOM initial inversion failure advanced final_time.");
    }
    if (result.diagnostics.source_diagnostics.accepted || result.diagnostics.first_reconstruction_diagnostics.success ||
        result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, "CQMOM initial inversion failure reported later-stage completion.");
    }
}

void verify_cqmom_dense_continuous_failure_result(const std::string& case_name, const qmomkit::TestDenseCqmomTransportResult& result,
                                                  const qmomkit::RawMomentTensor& initial_moments, const std::string& expected_message_fragment,
                                                  VerificationResult& verification)
{
    // Dense transport is the public workflow for callers that own full raw
    // moment tensors. A velocity callback failure should have the same
    // high-level contract as CqmomState transport: continuous_failed, first
    // continuous half-step, and no source/post-reconstruction progress.
    if (result.ok())
    {
        append_failure(verification, case_name + " unexpectedly reported dense CQMOM transport success.");
    }
    if (result.status != qmomkit::TransportStatus::continuous_failed)
    {
        append_failure(verification, case_name + " did not report TransportStatus::continuous_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::first_continuous_half_step)
    {
        append_failure(verification, case_name + " did not identify the first continuous half-step as the failed stage.");
    }
    if (result.diagnostics.first_continuous_status != qmomkit::NodeAdvectionStatus::velocity_failed)
    {
        append_failure(verification, case_name + " did not report NodeAdvectionStatus::velocity_failed.");
    }
    if (!contains_text(result.message, expected_message_fragment) || !contains_text(result.diagnostics.message, expected_message_fragment) ||
        !contains_text(result.diagnostics.first_continuous_diagnostics.message, expected_message_fragment))
    {
        append_failure(verification, case_name + " did not preserve the velocity-callback failure message.");
    }
    if (std::abs(result.final_time) > 1.0e-14)
    {
        append_failure(verification, case_name + " advanced final_time despite dense continuous transport failure.");
    }
    if (result.moments.extents != initial_moments.extents || !vectors_equal(result.moments.values, initial_moments.values, 0.0))
    {
        append_failure(verification, case_name + " changed dense moments during a failed continuous substep.");
    }
    if (result.diagnostics.source_diagnostics.accepted || result.diagnostics.source_diagnostics.rejected ||
        result.diagnostics.first_reconstruction_diagnostics.success || result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, case_name + " reported later-stage diagnostics after dense continuous failure.");
    }
}

void verify_cqmom_lean_continuous_failure_result(const std::string& case_name, const qmomkit::TestLeanCqmomTransportResult& result,
                                                 const qmomkit::TestLeanCqmomTransportState& initial_state,
                                                 const std::string& expected_message_fragment, VerificationResult& verification)
{
    // Explicit lean transport stores only the required CQMOM moments. The
    // failure contract should still identify the continuous node-advection
    // stage and leave the lean moment state unchanged.
    if (result.ok())
    {
        append_failure(verification, case_name + " unexpectedly reported lean CQMOM transport success.");
    }
    if (result.status != qmomkit::TransportStatus::continuous_failed)
    {
        append_failure(verification, case_name + " did not report TransportStatus::continuous_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::first_continuous_half_step)
    {
        append_failure(verification, case_name + " did not identify the first continuous half-step as the failed stage.");
    }
    if (result.diagnostics.first_continuous_status != qmomkit::NodeAdvectionStatus::velocity_failed)
    {
        append_failure(verification, case_name + " did not report NodeAdvectionStatus::velocity_failed.");
    }
    if (!contains_text(result.message, expected_message_fragment) || !contains_text(result.diagnostics.message, expected_message_fragment) ||
        !contains_text(result.diagnostics.first_continuous_diagnostics.message, expected_message_fragment))
    {
        append_failure(verification, case_name + " did not preserve the velocity-callback failure message.");
    }
    if (std::abs(result.final_time) > 1.0e-14)
    {
        append_failure(verification, case_name + " advanced final_time despite lean continuous transport failure.");
    }
    if (result.state.moments.layout.required_indices != initial_state.moments.layout.required_indices ||
        !vectors_equal(result.state.moments.values, initial_state.moments.values, 0.0))
    {
        append_failure(verification, case_name + " changed lean moments during a failed continuous substep.");
    }
    if (result.diagnostics.source_diagnostics.accepted || result.diagnostics.source_diagnostics.rejected ||
        result.diagnostics.first_reconstruction_diagnostics.success || result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, case_name + " reported later-stage diagnostics after lean continuous failure.");
    }
}

void verify_cqmom_dense_source_failure_result(const std::string& case_name, const qmomkit::TestDenseCqmomTransportResult& result,
                                              const std::string& expected_message_fragment, VerificationResult& verification)
{
    // Dense transport uses the generic raw-moment source callback. Returning a
    // wrong-sized or non-finite rate vector is a user model error and should
    // stop at the discontinuous source step.
    if (result.ok())
    {
        append_failure(verification, case_name + " unexpectedly reported dense source transport success.");
    }
    if (result.status != qmomkit::TransportStatus::source_failed)
    {
        append_failure(verification, case_name + " did not report TransportStatus::source_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::discontinuous_source_step)
    {
        append_failure(verification, case_name + " did not identify the discontinuous source step as the failed stage.");
    }
    if (result.diagnostics.source_status != qmomkit::MomentSourceUpdateStatus::source_failed)
    {
        append_failure(verification, case_name + " did not report MomentSourceUpdateStatus::source_failed.");
    }
    if (!contains_text(result.message, expected_message_fragment) || !contains_text(result.diagnostics.message, expected_message_fragment) ||
        !contains_text(result.diagnostics.source_diagnostics.message, expected_message_fragment))
    {
        append_failure(verification, case_name + " did not preserve the dense source failure message.");
    }
    if (result.diagnostics.source_diagnostics.accepted || result.diagnostics.source_diagnostics.rejected)
    {
        append_failure(verification, case_name + " reported accepted/rejected source integration after callback failure.");
    }
    if (result.diagnostics.post_source_cqmom_status != qmomkit::CqmomStatus::success || !result.diagnostics.post_source_cqmom_diagnostics.empty() ||
        result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, case_name + " reported post-source progress after dense source failure.");
    }
}

void verify_cqmom_lean_source_failure_result(const std::string& case_name, const qmomkit::TestLeanCqmomTransportResult& result,
                                             const std::string& expected_message_fragment, VerificationResult& verification)
{
    // Lean transport uses the CQMOM-native source callback. The high-level
    // failure contract should match CqmomState native-source failure handling.
    if (result.ok())
    {
        append_failure(verification, case_name + " unexpectedly reported lean source transport success.");
    }
    if (result.status != qmomkit::TransportStatus::source_failed)
    {
        append_failure(verification, case_name + " did not report TransportStatus::source_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::discontinuous_source_step)
    {
        append_failure(verification, case_name + " did not identify the discontinuous source step as the failed stage.");
    }
    if (result.diagnostics.source_status != qmomkit::MomentSourceUpdateStatus::source_failed)
    {
        append_failure(verification, case_name + " did not report MomentSourceUpdateStatus::source_failed.");
    }
    if (!contains_text(result.message, expected_message_fragment) || !contains_text(result.diagnostics.message, expected_message_fragment) ||
        !contains_text(result.diagnostics.source_diagnostics.message, expected_message_fragment))
    {
        append_failure(verification, case_name + " did not preserve the lean source failure message.");
    }
    if (result.diagnostics.source_diagnostics.accepted || result.diagnostics.source_diagnostics.rejected)
    {
        append_failure(verification, case_name + " reported accepted/rejected source integration after callback failure.");
    }
    if (result.diagnostics.post_source_cqmom_status != qmomkit::CqmomStatus::success || !result.diagnostics.post_source_cqmom_diagnostics.empty() ||
        result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, case_name + " reported post-source progress after lean source failure.");
    }
}

void verify_cqmom_dense_post_source_inversion_failure_result(const qmomkit::TestDenseCqmomTransportResult& result, const double expected_failure_time,
                                                             VerificationResult& verification)
{
    // Dense full-tensor transport integrates a finite source successfully, but
    // the manufactured source drives total mass negative. The following CQMOM
    // inversion should fail with diagnostics and stop before final advection.
    if (result.ok())
    {
        append_failure(verification, "dense post-source inversion failure unexpectedly reported transport success.");
    }
    if (result.status != qmomkit::TransportStatus::inversion_failed)
    {
        append_failure(verification, "dense post-source inversion failure did not report TransportStatus::inversion_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::post_source_inversion)
    {
        append_failure(verification, "dense post-source inversion failure did not identify post_source_inversion as the failed stage.");
    }
    if (result.diagnostics.source_status != qmomkit::MomentSourceUpdateStatus::success ||
        result.diagnostics.source_diagnostics.ssp_rk3_status != qmomkit::SspRk3Status::success || !result.diagnostics.source_diagnostics.accepted)
    {
        append_failure(verification, "dense post-source inversion failure did not preserve successful source diagnostics.");
    }
    if (result.diagnostics.post_source_cqmom_status == qmomkit::CqmomStatus::success || result.diagnostics.post_source_cqmom_diagnostics.empty())
    {
        append_failure(verification, "dense post-source inversion failure did not preserve failed CQMOM diagnostics.");
    }
    if (result.message.empty() || result.diagnostics.message.empty() ||
        (!contains_text(result.message, "m0") && !contains_text(result.message, "realizable")))
    {
        append_failure(verification, "dense post-source inversion failure did not provide a meaningful reconstruction message.");
    }
    if (std::abs(result.final_time - expected_failure_time) > 1.0e-14)
    {
        append_failure(verification, "dense post-source inversion failure did not stop at the source-stage time.");
    }
    if (result.diagnostics.final_reconstruction_diagnostics.success || result.diagnostics.second_continuous_diagnostics.step_count > 0)
    {
        append_failure(verification, "dense post-source inversion failure reported later-stage completion.");
    }
}

void verify_cqmom_lean_post_source_inversion_failure_result(const qmomkit::TestLeanCqmomTransportResult& result, const double expected_failure_time,
                                                            VerificationResult& verification)
{
    // Lean transport sees the same post-source realizability event through the
    // selected required-moment representation. The public status contract must
    // match dense transport even though the stored moment layout differs.
    if (result.ok())
    {
        append_failure(verification, "lean post-source inversion failure unexpectedly reported transport success.");
    }
    if (result.status != qmomkit::TransportStatus::inversion_failed)
    {
        append_failure(verification, "lean post-source inversion failure did not report TransportStatus::inversion_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::post_source_inversion)
    {
        append_failure(verification, "lean post-source inversion failure did not identify post_source_inversion as the failed stage.");
    }
    if (result.diagnostics.source_status != qmomkit::MomentSourceUpdateStatus::success ||
        result.diagnostics.source_diagnostics.ssp_rk3_status != qmomkit::SspRk3Status::success || !result.diagnostics.source_diagnostics.accepted)
    {
        append_failure(verification, "lean post-source inversion failure did not preserve successful source diagnostics.");
    }
    if (result.diagnostics.post_source_cqmom_status == qmomkit::CqmomStatus::success || result.diagnostics.post_source_cqmom_diagnostics.empty())
    {
        append_failure(verification, "lean post-source inversion failure did not preserve failed CQMOM diagnostics.");
    }
    if (result.message.empty() || result.diagnostics.message.empty() ||
        (!contains_text(result.message, "m0") && !contains_text(result.message, "realizable")))
    {
        append_failure(verification, "lean post-source inversion failure did not provide a meaningful reconstruction message.");
    }
    if (std::abs(result.final_time - expected_failure_time) > 1.0e-14)
    {
        append_failure(verification, "lean post-source inversion failure did not stop at the source-stage time.");
    }
    if (result.diagnostics.final_reconstruction_diagnostics.success || result.diagnostics.second_continuous_diagnostics.step_count > 0)
    {
        append_failure(verification, "lean post-source inversion failure reported later-stage completion.");
    }
}

void verify_cqmom_dense_initial_inversion_failure_result(const qmomkit::TestDenseCqmomTransportResult& result, VerificationResult& verification)
{
    // Dense initial failure represents a full raw tensor that cannot produce a
    // realizable CQMOM quadrature. Transport should stop before continuous or
    // source stages and preserve failed CQMOM diagnostics.
    if (result.ok())
    {
        append_failure(verification, "dense initial inversion failure unexpectedly reported transport success.");
    }
    if (result.status != qmomkit::TransportStatus::inversion_failed)
    {
        append_failure(verification, "dense initial inversion failure did not report TransportStatus::inversion_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::initial_inversion)
    {
        append_failure(verification, "dense initial inversion failure did not identify initial_inversion as the failed stage.");
    }
    if (result.diagnostics.initial_cqmom_status == qmomkit::CqmomStatus::success || result.diagnostics.initial_cqmom_diagnostics.empty())
    {
        append_failure(verification, "dense initial inversion failure did not preserve failed CQMOM diagnostics.");
    }
    if (result.message.empty() || result.diagnostics.message.empty() ||
        (!contains_text(result.message, "m0") && !contains_text(result.message, "realizable")))
    {
        append_failure(verification, "dense initial inversion failure did not provide a meaningful realizability message.");
    }
    if (std::abs(result.final_time) > 1.0e-14 || result.diagnostics.source_diagnostics.accepted ||
        result.diagnostics.first_reconstruction_diagnostics.success || result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, "dense initial inversion failure reported time advance or later-stage completion.");
    }
}

void verify_cqmom_lean_initial_inversion_failure_result(const qmomkit::TestLeanCqmomTransportResult& result, VerificationResult& verification)
{
    // Lean initial failure is the same user-visible event with only required
    // moments stored. It protects the selected-layout path from reporting
    // misleading progress after reconstruction fails.
    if (result.ok())
    {
        append_failure(verification, "lean initial inversion failure unexpectedly reported transport success.");
    }
    if (result.status != qmomkit::TransportStatus::inversion_failed)
    {
        append_failure(verification, "lean initial inversion failure did not report TransportStatus::inversion_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::initial_inversion)
    {
        append_failure(verification, "lean initial inversion failure did not identify initial_inversion as the failed stage.");
    }
    if (result.diagnostics.initial_cqmom_status == qmomkit::CqmomStatus::success || result.diagnostics.initial_cqmom_diagnostics.empty())
    {
        append_failure(verification, "lean initial inversion failure did not preserve failed CQMOM diagnostics.");
    }
    if (result.message.empty() || result.diagnostics.message.empty() ||
        (!contains_text(result.message, "m0") && !contains_text(result.message, "realizable")))
    {
        append_failure(verification, "lean initial inversion failure did not provide a meaningful realizability message.");
    }
    if (std::abs(result.final_time) > 1.0e-14 || result.diagnostics.source_diagnostics.accepted ||
        result.diagnostics.first_reconstruction_diagnostics.success || result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, "lean initial inversion failure reported time advance or later-stage completion.");
    }
}

void verify_dense_helper_based_common_cqmom_transport(const qmomkit::CqmomDenseCase& dense_case, VerificationResult& verification)
{
    qmomkit::CqmomDenseCase policy_case = dense_case;
    policy_case.qmom_options.adaptive = true;
    policy_case.adaptive_branch_counts = true;
    policy_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::basic;
    policy_case.vandermonde_solver = qmomkit::VandermondeSolver::dense;
    qmomkit::TestDenseCqmomTransportOptions options = cqmom_options(policy_case);

    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
    const qmomkit::TestDenseCqmomTransportResult result =
        qmomkit::advance_transport_cqmom(policy_case.moments, policy_case.hierarchy, 0.0, 0.02, zero_velocity, zero_source, options);
    if (!result.ok())
    {
        append_failure(verification, "helper-based dense CQMOM policy probe failed: " + result.message);
    }
    else
    {
        verify_helper_based_common_cqmom_policies("helper-based dense CQMOM", result.diagnostics, policy_case.hierarchy.number_of_dimensions(),
                                                  verification);
    }

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = 0.02;
    summary.steps = 1;
    summary.final_time = result.final_time;
    summary.retry_count = result.diagnostics.source_retries;
    summary.success = result.ok() && result.diagnostics.cqmom_diagnostic_level == qmomkit::CqmomDiagnosticLevel::basic &&
                      result.diagnostics.cqmom_vandermonde_solver == qmomkit::VandermondeSolver::dense &&
                      cqmom_diagnostics_enable_adaptive_branch_counts(result.diagnostics.initial_cqmom_diagnostics) &&
                      cqmom_diagnostics_request_solver(result.diagnostics.initial_cqmom_diagnostics, policy_case.hierarchy.number_of_dimensions(),
                                                       qmomkit::VandermondeSolver::dense);
    verification.runs.push_back(summary);
}

void verify_lean_helper_based_common_cqmom_transport(const qmomkit::CqmomLeanCase& lean_case, VerificationResult& verification)
{
    qmomkit::CqmomLeanCase policy_case = lean_case;
    policy_case.qmom_options.adaptive = true;
    policy_case.adaptive_branch_counts = true;
    policy_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::basic;
    policy_case.vandermonde_solver = qmomkit::VandermondeSolver::dense;
    qmomkit::TestLeanCqmomTransportOptions options = cqmom_lean_options(policy_case);

    const qmomkit::TestLeanCqmomTransportState state{policy_case.moments};
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };
    const qmomkit::TestLeanCqmomTransportResult result = qmomkit::advance_transport_cqmom(state, 0.0, 0.02, zero_velocity, zero_source, options);
    if (!result.ok())
    {
        append_failure(verification, "helper-based lean CQMOM policy probe failed: " + result.message);
    }
    else
    {
        verify_helper_based_common_cqmom_policies("helper-based lean CQMOM", result.diagnostics,
                                                  policy_case.moments.layout.hierarchy.number_of_dimensions(), verification);
    }

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = 0.02;
    summary.steps = 1;
    summary.final_time = result.final_time;
    summary.retry_count = result.diagnostics.source_retries;
    summary.success =
        result.ok() && result.diagnostics.cqmom_diagnostic_level == qmomkit::CqmomDiagnosticLevel::basic &&
        result.diagnostics.cqmom_vandermonde_solver == qmomkit::VandermondeSolver::dense &&
        cqmom_diagnostics_enable_adaptive_branch_counts(result.diagnostics.initial_cqmom_diagnostics) &&
        cqmom_diagnostics_request_solver(result.diagnostics.initial_cqmom_diagnostics, policy_case.moments.layout.hierarchy.number_of_dimensions(),
                                         qmomkit::VandermondeSolver::dense);
    verification.runs.push_back(summary);
}

VerificationResult run_cqmom_transport_failures(const std::filesystem::path& configuration_path)
{
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load CQMOM transport failure configuration.");
        return verification;
    }

    const qmomkit::CqmomCase& cqmom_case = configured.value();
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_case);
    if (!initial.ok())
    {
        append_failure(verification, "initial CQMOM transport failure inversion failed: " + initial.message);
        return verification;
    }

    qmomkit::CqmomTransportOptions options = cqmom_state_options();
    options.scheme = qmomkit::TransportScheme::strang;
    constexpr double dt = 0.08;
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    struct ContinuousFailureCase
    {
        std::string name;
        qmomkit::PhaseSpaceVelocity velocity;
        std::string expected_message_fragment;
    };

    const std::vector<ContinuousFailureCase> continuous_failure_cases{{"CQMOM wrong velocity dimension",
                                                                       [](double, const std::vector<double>&)
                                                                       {
                                                                           // Simulates a CQMOM user callback returning one velocity component
                                                                           // for a two-dimensional phase space.
                                                                           return std::vector<double>{0.0};
                                                                       },
                                                                       "Phase-space velocity callback returned the wrong dimension."},
                                                                      {"CQMOM non-finite velocity",
                                                                       [](double, const std::vector<double>&)
                                                                       {
                                                                           // Simulates a coupled model producing a non-finite node velocity.
                                                                           return std::vector<double>{0.0, std::numeric_limits<double>::infinity()};
                                                                       },
                                                                       "Phase-space velocity callback returned a non-finite value."}};

    for (const ContinuousFailureCase& failure_case : continuous_failure_cases)
    {
        qmomkit::CqmomState state = cqmom_case.initial_state();
        const qmomkit::CqmomTransportResult result = qmomkit::advance_transport_cqmom(state, 0.0, dt, failure_case.velocity, zero_source, options);
        verify_cqmom_continuous_failure_result(failure_case.name, result, failure_case.expected_message_fragment, verification);

        RunSummary summary;
        summary.scheme = options.scheme;
        summary.dt = dt;
        summary.steps = 1;
        summary.final_time = result.final_time;
        summary.retry_count = result.diagnostics.source_retries;
        summary.success = !result.ok() && result.status == qmomkit::TransportStatus::continuous_failed &&
                          result.diagnostics.failed_stage == qmomkit::TransportStage::first_continuous_half_step &&
                          result.diagnostics.first_continuous_status == qmomkit::NodeAdvectionStatus::velocity_failed;
        verification.runs.push_back(summary);
    }

    qmomkit::CqmomTransportOptions rk_failure_options = options;
    rk_failure_options.node_options.scheme = qmomkit::NodeAdvectionScheme::adaptive_rk23;
    rk_failure_options.node_options.rk23_options.initial_step = 0.04;
    rk_failure_options.node_options.rk23_options.max_step = 0.04;
    rk_failure_options.node_options.rk23_options.min_step = 0.02;
    const auto stiff_velocity = [](double, const std::vector<double>& xi)
    {
        // Valid two-dimensional velocity with deliberately stiff diagonal
        // growth/decay. The public RK23 diagnostics should explain that the
        // adaptive node step underflowed, not that the callback itself failed.
        return std::vector<double>{1000.0 * xi[0], -1000.0 * xi[1]};
    };
    qmomkit::CqmomState state = cqmom_case.initial_state();
    const qmomkit::CqmomTransportResult rk_failure =
        qmomkit::advance_transport_cqmom(state, 0.0, dt, stiff_velocity, zero_source, rk_failure_options);
    verify_cqmom_rk23_minimum_step_failure_result(rk_failure, verification);

    RunSummary rk_summary;
    rk_summary.scheme = rk_failure_options.scheme;
    rk_summary.dt = dt;
    rk_summary.steps = 1;
    rk_summary.final_time = rk_failure.final_time;
    rk_summary.retry_count = rk_failure.diagnostics.source_retries;
    rk_summary.success =
        !rk_failure.ok() && rk_failure.status == qmomkit::TransportStatus::continuous_failed &&
        rk_failure.diagnostics.failed_stage == qmomkit::TransportStage::first_continuous_half_step &&
        rk_failure.diagnostics.first_continuous_status == qmomkit::NodeAdvectionStatus::integration_failed &&
        rk_failure.diagnostics.first_continuous_diagnostics.adaptive_rk23_status == qmomkit::AdaptiveRk23Status::step_size_underflow &&
        rk_failure.diagnostics.first_continuous_diagnostics.rejected_step_count > 0 &&
        rk_failure.diagnostics.first_continuous_diagnostics.step_count == 0;
    verification.runs.push_back(rk_summary);

    state = cqmom_case.initial_state();
    const auto nonrealizable_source = [](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rates(context.size(), 0.0);
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            if (context.index(entry) == qmomkit::MultiIndex{0, 0})
            {
                rates[entry] = -20.0;
            }
        }
        return qmomkit::CqmomMomentRate{rates};
    };
    const qmomkit::CqmomTransportResult post_source_inversion =
        qmomkit::advance_transport_cqmom(state, 0.0, dt, zero_velocity, nonrealizable_source, options);
    verify_cqmom_post_source_inversion_failure_result(post_source_inversion, 0.5 * dt, verification);

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = dt;
    summary.steps = 1;
    summary.final_time = post_source_inversion.final_time;
    summary.retry_count = post_source_inversion.diagnostics.source_retries;
    summary.success = !post_source_inversion.ok() && post_source_inversion.status == qmomkit::TransportStatus::inversion_failed &&
                      post_source_inversion.diagnostics.failed_stage == qmomkit::TransportStage::post_source_inversion &&
                      post_source_inversion.diagnostics.post_source_cqmom_status != qmomkit::CqmomStatus::success;
    verification.runs.push_back(summary);
    return verification;
}

VerificationResult run_cqmom_transport_inversion_failures(const std::filesystem::path& configuration_path)
{
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load CQMOM transport inversion-failure configuration.");
        return verification;
    }

    const qmomkit::CqmomCase& cqmom_case = configured.value();
    qmomkit::CqmomState state = cqmom_case.initial_state();
    qmomkit::CqmomTransportOptions options = cqmom_state_options();
    options.scheme = qmomkit::TransportScheme::strang;
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::CqmomTransportResult result = qmomkit::advance_transport_cqmom(state, 0.0, 0.08, zero_velocity, cqmom_state_zero_source, options);
    verify_cqmom_initial_inversion_failure_result(result, verification);

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = 0.08;
    summary.steps = 1;
    summary.final_time = result.final_time;
    summary.retry_count = result.diagnostics.source_retries;
    summary.success = !result.ok() && result.status == qmomkit::TransportStatus::inversion_failed &&
                      result.diagnostics.failed_stage == qmomkit::TransportStage::initial_inversion &&
                      result.diagnostics.initial_cqmom_status != qmomkit::CqmomStatus::success;
    verification.runs.push_back(summary);
    return verification;
}

VerificationResult run_cqmom_dense_transport_failures(const std::filesystem::path& configuration_path)
{
    // Dense-backed CqmomState failure propagation covers callers that start from
    // complete raw moment tensors. The cases mirror CQMOM transport-through-state
    // failures while verifying that dense diagnostic moments remain available on
    // failure.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom_dense(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load dense CQMOM transport failure configuration.");
        return verification;
    }

    const qmomkit::CqmomDenseCase& dense_case = configured.value();
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(dense_case);
    if (!initial.ok())
    {
        append_failure(verification, "initial dense CQMOM transport failure inversion failed: " + initial.message);
        return verification;
    }
    verify_dense_helper_based_common_cqmom_transport(dense_case, verification);

    qmomkit::TestDenseCqmomTransportOptions options = cqmom_options(dense_case);
    options.scheme = qmomkit::TransportScheme::strang;
    constexpr double dt = 0.08;
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    struct ContinuousFailureCase
    {
        std::string name;
        qmomkit::PhaseSpaceVelocity velocity;
        std::string expected_message_fragment;
    };
    const std::vector<ContinuousFailureCase> continuous_failure_cases{{"dense CQMOM wrong velocity dimension",
                                                                       [](double, const std::vector<double>&) { return std::vector<double>{0.0}; },
                                                                       "Phase-space velocity callback returned the wrong dimension."},
                                                                      {"dense CQMOM non-finite velocity", [](double, const std::vector<double>&)
                                                                       { return std::vector<double>{0.0, std::numeric_limits<double>::infinity()}; },
                                                                       "Phase-space velocity callback returned a non-finite value."}};

    for (const ContinuousFailureCase& failure_case : continuous_failure_cases)
    {
        const qmomkit::TestDenseCqmomTransportResult result =
            qmomkit::advance_transport_cqmom(dense_case.moments, dense_case.hierarchy, 0.0, dt, failure_case.velocity, zero_source, options);
        verify_cqmom_dense_continuous_failure_result(failure_case.name, result, dense_case.moments, failure_case.expected_message_fragment,
                                                     verification);

        RunSummary summary;
        summary.scheme = options.scheme;
        summary.dt = dt;
        summary.steps = 1;
        summary.final_time = result.final_time;
        summary.retry_count = result.diagnostics.source_retries;
        summary.success = !result.ok() && result.status == qmomkit::TransportStatus::continuous_failed &&
                          result.diagnostics.failed_stage == qmomkit::TransportStage::first_continuous_half_step &&
                          result.diagnostics.first_continuous_status == qmomkit::NodeAdvectionStatus::velocity_failed;
        verification.runs.push_back(summary);
    }

    const qmomkit::CqmomMomentSource wrong_size_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size() + 1, 0.0)}; };
    const qmomkit::TestDenseCqmomTransportResult source_failure =
        qmomkit::advance_transport_cqmom(dense_case.moments, dense_case.hierarchy, 0.0, dt, zero_velocity, wrong_size_source, options);
    verify_cqmom_dense_source_failure_result("dense CQMOM wrong source rate size", source_failure,
                                             "CQMOM source rate values must match the source context size.", verification);

    RunSummary source_summary;
    source_summary.scheme = options.scheme;
    source_summary.dt = dt;
    source_summary.steps = 1;
    source_summary.final_time = source_failure.final_time;
    source_summary.retry_count = source_failure.diagnostics.source_retries;
    source_summary.success = !source_failure.ok() && source_failure.status == qmomkit::TransportStatus::source_failed &&
                             source_failure.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step;
    verification.runs.push_back(source_summary);

    const qmomkit::CqmomMomentSource nonrealizable_source = [](double, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> rates(moments.size(), 0.0);
        for (std::size_t entry = 0; entry < moments.size(); ++entry)
        {
            if (moments.index(entry) == qmomkit::MultiIndex{0, 0})
            {
                rates[entry] = -20.0;
            }
        }
        return qmomkit::CqmomMomentRate{rates};
    };
    const qmomkit::TestDenseCqmomTransportResult post_source_inversion =
        qmomkit::advance_transport_cqmom(dense_case.moments, dense_case.hierarchy, 0.0, dt, zero_velocity, nonrealizable_source, options);
    verify_cqmom_dense_post_source_inversion_failure_result(post_source_inversion, 0.5 * dt, verification);

    RunSummary inversion_summary;
    inversion_summary.scheme = options.scheme;
    inversion_summary.dt = dt;
    inversion_summary.steps = 1;
    inversion_summary.final_time = post_source_inversion.final_time;
    inversion_summary.retry_count = post_source_inversion.diagnostics.source_retries;
    inversion_summary.success = !post_source_inversion.ok() && post_source_inversion.status == qmomkit::TransportStatus::inversion_failed &&
                                post_source_inversion.diagnostics.failed_stage == qmomkit::TransportStage::post_source_inversion &&
                                post_source_inversion.diagnostics.post_source_cqmom_status != qmomkit::CqmomStatus::success;
    verification.runs.push_back(inversion_summary);
    return verification;
}

VerificationResult run_cqmom_lean_transport_failures(const std::filesystem::path& configuration_path)
{
    // Lean CQMOM failure propagation covers callers that transport only the
    // required recursive moment subset. The same callback and realizability
    // failures should map to the same public statuses as dense transport.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom_lean(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load lean CQMOM transport failure configuration.");
        return verification;
    }

    const qmomkit::CqmomLeanCase& lean_case = configured.value();
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(lean_case);
    if (!initial.ok())
    {
        append_failure(verification, "initial lean CQMOM transport failure inversion failed: " + initial.message);
        return verification;
    }
    verify_lean_helper_based_common_cqmom_transport(lean_case, verification);

    qmomkit::TestLeanCqmomTransportOptions options = cqmom_lean_options(lean_case);
    options.scheme = qmomkit::TransportScheme::strang;
    const qmomkit::TestLeanCqmomTransportState initial_state{lean_case.moments};
    constexpr double dt = 0.08;
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    struct ContinuousFailureCase
    {
        std::string name;
        qmomkit::PhaseSpaceVelocity velocity;
        std::string expected_message_fragment;
    };
    const std::vector<ContinuousFailureCase> continuous_failure_cases{{"lean CQMOM wrong velocity dimension",
                                                                       [](double, const std::vector<double>&) { return std::vector<double>{0.0}; },
                                                                       "Phase-space velocity callback returned the wrong dimension."},
                                                                      {"lean CQMOM non-finite velocity", [](double, const std::vector<double>&)
                                                                       { return std::vector<double>{0.0, std::numeric_limits<double>::infinity()}; },
                                                                       "Phase-space velocity callback returned a non-finite value."}};

    for (const ContinuousFailureCase& failure_case : continuous_failure_cases)
    {
        const qmomkit::TestLeanCqmomTransportResult result =
            qmomkit::advance_transport_cqmom(initial_state, 0.0, dt, failure_case.velocity, zero_source, options);
        verify_cqmom_lean_continuous_failure_result(failure_case.name, result, initial_state, failure_case.expected_message_fragment, verification);

        RunSummary summary;
        summary.scheme = options.scheme;
        summary.dt = dt;
        summary.steps = 1;
        summary.final_time = result.final_time;
        summary.retry_count = result.diagnostics.source_retries;
        summary.success = !result.ok() && result.status == qmomkit::TransportStatus::continuous_failed &&
                          result.diagnostics.failed_stage == qmomkit::TransportStage::first_continuous_half_step &&
                          result.diagnostics.first_continuous_status == qmomkit::NodeAdvectionStatus::velocity_failed;
        verification.runs.push_back(summary);
    }

    const auto wrong_size_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size() + 1, 0.0)}; };
    const qmomkit::TestLeanCqmomTransportResult source_failure =
        qmomkit::advance_transport_cqmom(initial_state, 0.0, dt, zero_velocity, wrong_size_source, options);
    verify_cqmom_lean_source_failure_result("lean CQMOM wrong source rate size", source_failure,
                                            "CQMOM source rate values must match the source context size.", verification);

    RunSummary source_summary;
    source_summary.scheme = options.scheme;
    source_summary.dt = dt;
    source_summary.steps = 1;
    source_summary.final_time = source_failure.final_time;
    source_summary.retry_count = source_failure.diagnostics.source_retries;
    source_summary.success = !source_failure.ok() && source_failure.status == qmomkit::TransportStatus::source_failed &&
                             source_failure.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step;
    verification.runs.push_back(source_summary);

    const auto nonrealizable_source = [](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rates(context.size(), 0.0);
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            if (context.index(entry) == qmomkit::MultiIndex{0, 0})
            {
                rates[entry] = -20.0;
            }
        }
        return qmomkit::CqmomMomentRate{rates};
    };
    const qmomkit::TestLeanCqmomTransportResult post_source_inversion =
        qmomkit::advance_transport_cqmom(initial_state, 0.0, dt, zero_velocity, nonrealizable_source, options);
    verify_cqmom_lean_post_source_inversion_failure_result(post_source_inversion, 0.5 * dt, verification);

    RunSummary inversion_summary;
    inversion_summary.scheme = options.scheme;
    inversion_summary.dt = dt;
    inversion_summary.steps = 1;
    inversion_summary.final_time = post_source_inversion.final_time;
    inversion_summary.retry_count = post_source_inversion.diagnostics.source_retries;
    inversion_summary.success = !post_source_inversion.ok() && post_source_inversion.status == qmomkit::TransportStatus::inversion_failed &&
                                post_source_inversion.diagnostics.failed_stage == qmomkit::TransportStage::post_source_inversion &&
                                post_source_inversion.diagnostics.post_source_cqmom_status != qmomkit::CqmomStatus::success;
    verification.runs.push_back(inversion_summary);
    return verification;
}

VerificationResult run_cqmom_dense_transport_inversion_failures(const std::filesystem::path& configuration_path)
{
    // Dense initial inversion failure starts from a configured full tensor that
    // is invalid for CQMOM. This protects early exit before any user velocity or
    // source callback can alter diagnostics.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom_dense(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load dense CQMOM transport inversion-failure configuration.");
        return verification;
    }

    const qmomkit::CqmomDenseCase& dense_case = configured.value();
    qmomkit::TestDenseCqmomTransportOptions options = cqmom_options(dense_case);
    options.scheme = qmomkit::TransportScheme::strang;
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
    const qmomkit::TestDenseCqmomTransportResult result =
        qmomkit::advance_transport_cqmom(dense_case.moments, dense_case.hierarchy, 0.0, 0.08, zero_velocity, zero_source, options);
    verify_cqmom_dense_initial_inversion_failure_result(result, verification);

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = 0.08;
    summary.steps = 1;
    summary.final_time = result.final_time;
    summary.retry_count = result.diagnostics.source_retries;
    summary.success = !result.ok() && result.status == qmomkit::TransportStatus::inversion_failed &&
                      result.diagnostics.failed_stage == qmomkit::TransportStage::initial_inversion &&
                      result.diagnostics.initial_cqmom_status != qmomkit::CqmomStatus::success;
    verification.runs.push_back(summary);
    return verification;
}

VerificationResult run_cqmom_lean_transport_inversion_failures(const std::filesystem::path& configuration_path)
{
    // Lean initial inversion failure checks the selected-layout path for the
    // same early-exit behavior: failed reconstruction, no time advance, and no
    // later-stage diagnostics.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom_lean(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load lean CQMOM transport inversion-failure configuration.");
        return verification;
    }

    const qmomkit::CqmomLeanCase& lean_case = configured.value();
    qmomkit::TestLeanCqmomTransportOptions options = cqmom_lean_options(lean_case);
    options.scheme = qmomkit::TransportScheme::strang;
    const qmomkit::TestLeanCqmomTransportState state{lean_case.moments};
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::TestLeanCqmomTransportResult result =
        qmomkit::advance_transport_cqmom(state, 0.0, 0.08, zero_velocity, cqmom_state_zero_source, options);
    verify_cqmom_lean_initial_inversion_failure_result(result, verification);

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = 0.08;
    summary.steps = 1;
    summary.final_time = result.final_time;
    summary.retry_count = result.diagnostics.source_retries;
    summary.success = !result.ok() && result.status == qmomkit::TransportStatus::inversion_failed &&
                      result.diagnostics.failed_stage == qmomkit::TransportStage::initial_inversion &&
                      result.diagnostics.initial_cqmom_status != qmomkit::CqmomStatus::success;
    verification.runs.push_back(summary);
    return verification;
}

void verify_cqmom_native_source_failure_result(const std::string& case_name, const qmomkit::CqmomTransportResult& result,
                                               const std::string& expected_message_fragment, VerificationResult& verification)
{
    // User callback failures should stop during the discontinuous source step,
    // before post-source inversion or later transport stages can hide the
    // callback contract violation.
    if (result.ok())
    {
        append_failure(verification, case_name + " unexpectedly reported transport success.");
    }
    if (result.status != qmomkit::TransportStatus::source_failed)
    {
        append_failure(verification, case_name + " did not report TransportStatus::source_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::TransportStage::discontinuous_source_step)
    {
        append_failure(verification, case_name + " did not identify the discontinuous source step as the failed stage.");
    }
    if (result.diagnostics.source_status != qmomkit::MomentSourceUpdateStatus::source_failed)
    {
        append_failure(verification, case_name + " did not report MomentSourceUpdateStatus::source_failed.");
    }
    if (result.message.find(expected_message_fragment) == std::string::npos)
    {
        append_failure(verification, case_name + " result message did not explain the source callback problem.");
    }
    if (result.diagnostics.message.find(expected_message_fragment) == std::string::npos)
    {
        append_failure(verification, case_name + " transport diagnostics did not carry the source callback message.");
    }
    if (result.diagnostics.source_diagnostics.message.find(expected_message_fragment) == std::string::npos)
    {
        append_failure(verification, case_name + " source diagnostics did not carry the source callback message.");
    }
    if (result.diagnostics.source_diagnostics.accepted || result.diagnostics.source_diagnostics.rejected)
    {
        append_failure(verification, case_name + " source diagnostics reported accepted/rejected RK state after callback failure.");
    }
    if (result.diagnostics.post_source_cqmom_status != qmomkit::CqmomStatus::success || !result.diagnostics.post_source_cqmom_diagnostics.empty())
    {
        append_failure(verification, case_name + " unexpectedly reached post-source CQMOM reconstruction.");
    }
    if (result.diagnostics.final_reconstruction_diagnostics.success || result.diagnostics.second_continuous_diagnostics.step_count > 0 ||
        result.final_time <= 0.0)
    {
        append_failure(verification, case_name + " reported later-stage progress after source callback failure.");
    }
}

VerificationResult run_cqmom_native_source_failures(const std::filesystem::path& configuration_path)
{
    // Native CQMOM source failures model user physics callbacks that return a
    // rate vector with the wrong size or non-finite entries. The transport
    // result should fail at the source stage while preserving useful callback
    // diagnostics.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load CQMOM-native source failure configuration.");
        return verification;
    }

    const qmomkit::CqmomCase& cqmom_case = configured.value();
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_case);
    if (!initial.ok())
    {
        append_failure(verification, "initial CQMOM-native source failure inversion failed: " + initial.message);
        return verification;
    }

    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    qmomkit::CqmomTransportOptions options = cqmom_state_options();
    options.scheme = qmomkit::TransportScheme::strang;
    constexpr double dt = 0.08;

    struct FailureCase
    {
        std::string name;
        qmomkit::CqmomMomentSource source;
        std::string expected_message_fragment;
    };

    const std::vector<FailureCase> failure_cases{{"wrong rate size",
                                                  [](double, const qmomkit::CqmomSourceContext& context)
                                                  {
                                                      // Simulates a user callback that ignores the documented
                                                      // requirement to return one derivative per context entry.
                                                      return qmomkit::CqmomMomentRate{std::vector<double>(context.size() + 1, 0.0)};
                                                  },
                                                  "CQMOM source rate values must match the source context size."},
                                                 {"NaN rate",
                                                  [](double, const qmomkit::CqmomSourceContext& context)
                                                  {
                                                      std::vector<double> rates(context.size(), 0.0);
                                                      if (!rates.empty())
                                                      {
                                                          rates.front() = std::numeric_limits<double>::quiet_NaN();
                                                      }
                                                      return qmomkit::CqmomMomentRate{rates};
                                                  },
                                                  "CQMOM source rate values must be finite."},
                                                 {"infinite rate",
                                                  [](double, const qmomkit::CqmomSourceContext& context)
                                                  {
                                                      std::vector<double> rates(context.size(), 0.0);
                                                      if (!rates.empty())
                                                      {
                                                          rates.front() = std::numeric_limits<double>::infinity();
                                                      }
                                                      return qmomkit::CqmomMomentRate{rates};
                                                  },
                                                  "CQMOM source rate values must be finite."}};

    for (const FailureCase& failure_case : failure_cases)
    {
        qmomkit::CqmomState state = cqmom_case.initial_state();
        const qmomkit::CqmomTransportResult result = qmomkit::advance_transport_cqmom(state, 0.0, dt, zero_velocity, failure_case.source, options);
        verify_cqmom_native_source_failure_result(failure_case.name, result, failure_case.expected_message_fragment, verification);

        RunSummary summary;
        summary.scheme = options.scheme;
        summary.dt = dt;
        summary.steps = 1;
        summary.final_time = result.final_time;
        summary.retry_count = result.diagnostics.source_retries;
        summary.success = !result.ok() && result.status == qmomkit::TransportStatus::source_failed &&
                          result.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step &&
                          result.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::source_failed;
        verification.runs.push_back(summary);
    }
    return verification;
}

VerificationResult run_commuting_cqmom_state_adaptive_2d(const std::filesystem::path& configuration_path)
{
    // Adaptive CQMOM transport through CqmomState uses a deterministic affine
    // branch-reduction case with commuting continuous/source operators. It
    // verifies that ragged adaptive quadrature remains transportable and
    // reconstructs the exact analytical moment tensor after splitting.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load adaptive CQMOM state configuration.");
        return verification;
    }

    const qmomkit::CqmomCase& cqmom_case = configured.value();
    if (!cqmom_case.adaptive_branch_counts() || !cqmom_case.qmom_options().adaptive)
    {
        append_failure(verification, "adaptive CQMOM state verification requires adaptive_branch_counts and adaptive QMOM.");
        return verification;
    }

    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_case);
    if (!initial.ok())
    {
        append_failure(verification, "initial adaptive CQMOM state inversion failed: " + initial.message);
        return verification;
    }
    if (!initial.quadrature.adaptive_branch_counts_enabled() || initial.quadrature.has_rectangular_storage())
    {
        append_failure(verification, "initial adaptive CQMOM state quadrature did not use ragged adaptive storage.");
    }
    if (initial.quadrature.total_nodes() >= cqmom_case.hierarchy().total_nodes())
    {
        append_failure(verification, "initial adaptive CQMOM state quadrature did not reduce active node count below maximum capacity.");
    }
    if (initial.quadrature.terminal_paths().size() != initial.quadrature.weights().size())
    {
        append_failure(verification, "initial adaptive CQMOM state quadrature did not provide one terminal path per node.");
    }

    const std::vector<std::size_t> extents{4, 4};
    const qmomkit::RawMomentTensor initial_moments = qmomkit::reconstruct_raw_moments(initial.quadrature, extents);
    qmomkit::CqmomState state = cqmom_case.initial_state();
    const qmomkit::CqmomTransportOptions options = cqmom_state_options();
    const double dt = 0.04;
    const std::size_t steps = 5;
    const double elapsed = dt * static_cast<double>(steps);
    double time = 0.0;
    bool diagnostics_successful = true;
    std::size_t retry_count = 0;
    std::optional<qmomkit::CqmomTransportResult> result;
    for (std::size_t step = 0; step < steps; ++step)
    {
        result = qmomkit::advance_transport_cqmom(state, time, dt, cqmom_velocity, cqmom_state_zero_source, options);
        if (!result->ok())
        {
            append_failure(verification, "adaptive CQMOM transport through CqmomState failed: " + result->message);
            return verification;
        }
        if (!cqmom_diagnostics_are_successful(result->diagnostics))
        {
            append_failure(verification,
                           "adaptive CQMOM transport through CqmomState diagnostics did not report successful completion of every substage.");
            diagnostics_successful = false;
        }
        retry_count += result->diagnostics.source_retries;
        state = result->state;
        time = result->final_time;
    }

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = dt;
    summary.steps = steps;
    summary.final_time = time;
    summary.retry_count = retry_count;
    summary.success = result.has_value() && result->ok() && diagnostics_successful;
    summary.output_row = true;
    const std::vector<double> exact = cqmom_continuous_exact_values(initial_moments, elapsed);
    if (result->quadrature.has_value())
    {
        if (!result->quadrature->adaptive_branch_counts_enabled() || result->quadrature->has_rectangular_storage())
        {
            summary.success = false;
            append_failure(verification, "final adaptive CQMOM state quadrature did not preserve ragged adaptive storage.");
        }
        if (result->quadrature->total_nodes() > cqmom_case.hierarchy().total_nodes())
        {
            summary.success = false;
            append_failure(verification, "final adaptive CQMOM state quadrature exceeded configured maximum node capacity.");
        }
        if (result->quadrature->terminal_paths().size() != result->quadrature->weights().size())
        {
            summary.success = false;
            append_failure(verification, "final adaptive CQMOM state quadrature did not provide one terminal path per node.");
        }
        const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(*result->quadrature, extents);
        update_errors(summary.errors, reconstructed.values, exact);
        update_reconstruction_errors(summary.errors, reconstructed.values, exact);
        summary.error_norm = scaled_infinity_error(reconstructed.values, exact);
        print_vector_comparison("Final adaptive CQMOM state moments", reconstructed.values, exact);
    }
    else
    {
        summary.success = false;
        append_failure(verification, "adaptive CQMOM transport through CqmomState result did not retain final quadrature.");
    }
    if (std::abs(time - elapsed) > 1.0e-12)
    {
        summary.success = false;
        append_failure(verification, "adaptive CQMOM state final time does not match the requested final time.");
    }
    if (summary.errors.moment_scaled > 3.0e-6 || summary.errors.reconstruction_scaled > 3.0e-6)
    {
        summary.success = false;
        append_failure(verification, "adaptive CQMOM state commuting transport exceeded the analytical tolerance.");
    }
    verification.runs.push_back(summary);
    return verification;
}

void print_summary(const std::string& case_name, const VerificationResult& result)
{
    std::cout << std::setprecision(17);
    std::cout << "Case: " << case_name << '\n';
    std::cout << "scheme dt steps final_time error_norm observed_order max_abs max_scaled reconstruction_abs reconstruction_scaled m0_error "
                 "retry_count success\n";
    for (const RunSummary& run : result.runs)
    {
        std::cout << transport_scheme_name(run.scheme) << ' ' << run.dt << ' ' << run.steps << ' ' << run.final_time << ' ' << run.error_norm << ' '
                  << run.observed_order << ' ' << run.errors.moment_absolute << ' ' << run.errors.moment_scaled << ' '
                  << run.errors.reconstruction_absolute << ' ' << run.errors.reconstruction_scaled << ' ' << run.errors.conservation << ' '
                  << run.retry_count << ' ' << run.success << '\n';
    }
}

bool write_output_file(const std::filesystem::path& output_path, const std::string& case_name, const VerificationResult& result)
{
    std::ofstream output(output_path);
    if (!output)
    {
        return false;
    }
    output << std::setprecision(17);
    output << "case scheme dt step_count final_time error_norm observed_order zeroth_moment_error reconstruction_error retry_count success\n";
    for (const RunSummary& run : result.runs)
    {
        if (!run.output_row)
        {
            continue;
        }
        output << case_name << ' ' << transport_scheme_name(run.scheme) << ' ' << run.dt << ' ' << run.steps << ' ' << run.final_time << ' '
               << run.error_norm << ' ' << run.observed_order << ' ' << run.errors.conservation << ' ' << run.errors.reconstruction_scaled << ' '
               << run.retry_count << ' ' << run.success << '\n';
    }
    return output.good();
}

void print_usage()
{
    std::cerr << "Usage: transport_analytical_verification <configuration.qmom> <case_identifier> [--output <path>]\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3 && argc != 5)
    {
        print_usage();
        return 2;
    }

    const std::filesystem::path configuration_path = argv[1];
    const std::string case_name = argv[2];
    std::filesystem::path output_path;
    if (argc == 5)
    {
        if (std::string(argv[3]) != "--output")
        {
            print_usage();
            return 2;
        }
        output_path = argv[4];
    }

    VerificationResult result;
    if (case_name == "commuting_qmom")
    {
        result = run_commuting_qmom(configuration_path);
    }
    else if (case_name == "qmom_boundary_contribution")
    {
        result = run_qmom_boundary_contribution(configuration_path);
    }
    else if (case_name == "qmom_transport_failures")
    {
        result = run_qmom_transport_failures(configuration_path);
    }
    else if (case_name == "qmom_transport_inversion_failures")
    {
        result = run_qmom_transport_inversion_failures(configuration_path);
    }
    else if (case_name == "noncommuting_qmom")
    {
        result = run_noncommuting_qmom(configuration_path);
    }
    else if (case_name == "commuting_cqmom_2d")
    {
        result = run_commuting_cqmom_2d(configuration_path);
    }
    else if (case_name == "cqmom_boundary_contribution")
    {
        result = run_cqmom_boundary_contribution(configuration_path);
    }
    else if (case_name == "cqmom_native_source")
    {
        result = run_cqmom_native_source(configuration_path);
    }
    else if (case_name == "cqmom_combined_native_source")
    {
        result = run_cqmom_combined_native_source(configuration_path);
    }
    else if (case_name == "cqmom_native_source_failures")
    {
        result = run_cqmom_native_source_failures(configuration_path);
    }
    else if (case_name == "cqmom_transport_failures")
    {
        result = run_cqmom_transport_failures(configuration_path);
    }
    else if (case_name == "cqmom_transport_inversion_failures")
    {
        result = run_cqmom_transport_inversion_failures(configuration_path);
    }
    else if (case_name == "cqmom_dense_transport_failures")
    {
        result = run_cqmom_dense_transport_failures(configuration_path);
    }
    else if (case_name == "cqmom_dense_transport_inversion_failures")
    {
        result = run_cqmom_dense_transport_inversion_failures(configuration_path);
    }
    else if (case_name == "cqmom_lean_transport_failures")
    {
        result = run_cqmom_lean_transport_failures(configuration_path);
    }
    else if (case_name == "cqmom_lean_transport_inversion_failures")
    {
        result = run_cqmom_lean_transport_inversion_failures(configuration_path);
    }
    else if (case_name == "commuting_cqmom_state_adaptive_2d")
    {
        result = run_commuting_cqmom_state_adaptive_2d(configuration_path);
    }
    else
    {
        std::cerr << "Unknown transport analytical case '" << case_name << "'.\n";
        return 2;
    }

    print_summary(case_name, result);
    if (!output_path.empty() && !write_output_file(output_path, case_name, result))
    {
        std::cerr << "Unable to write transport analytical output file '" << output_path.string() << "'.\n";
        return 1;
    }
    if (!result.passed)
    {
        std::cerr << result.message << '\n';
        return 1;
    }
    return 0;
}
