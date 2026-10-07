// Phase-space evolution against closed-form solutions. This file holds several
// related studies sharing one set of helpers; the registered case names say
// which is which.
//
// What is verified
//   Commuting cases. The velocity is u(x) = a*x and the source removes each
//     moment order at its own rate, so both operators act on moment order k as
//     a multiple of M_k. They therefore commute, splitting introduces no error
//     at all, and all three orderings must return the same answer. This is an
//     exactness check, not a convergence one.
//
//   Non-commuting cases. The node velocity is a constant, which couples M_k to
//     M_(k-1) rather than to itself, alongside a decay. The operators no longer
//     commute, so splitting error appears and its size distinguishes the
//     schemes. This is where the splitting order is measured.
//
//   Boundary contributions and native sources, checked against the same kind of
//     closed-form reference, plus the failure paths: what CQMOM and QMOM
//     evolution report when a step or an inversion cannot succeed.
//
// How the reference is constructed
//   Written down, not computed. For the commuting case the moments follow
//
//       M_k(t) = M_k(0) * exp((k*a - lambda_k) * t),
//
//   with a = 0.2 and lambda_k = 0.05*k. For the non-commuting case a constant
//   node velocity b = 0.45 with decay gamma = 0.17 gives a triangular linear
//   system in the moments, solved here in closed form. Both expressions are
//   evaluated directly and depend on no part of the evolution machinery. The
//   starting moments come from chosen weighted points, so the initial state is
//   exact as well.
//
// Error metric and pass criteria
//   Maximum moment error relative to max(|expected|, 1). Observed splitting
//   order is the base-2 logarithm of the ratio of errors between successive
//   halvings of the time step.
//
//   Commuting cases are held to 2.0e-6 to 8.0e-6: the splitting is exact
//   here, so the remaining error is the time integrators', dominated by the
//   SSP-RK3 moment-source step, and a looser bound would hide a splitting
//   defect.
//
//   Non-commuting cases are held to 8.0e-4 for Strang and 3.0e-2 for the
//   first-order orderings, because at a fixed step those schemes really are
//   that far from the exact answer; the meaningful statement is the rate at
//   which the gap closes. The observed order must exceed 1.85 for Strang and
//   0.85 for the first-order orderings on the finest pair, five to eight
//   percent below the nominal 2 and 1; the CQMOM case requires at least 1.8
//   and 0.8 on its two finest pairs. A separate convergence summary requires
//   an order above 2.4 where a higher-order combination is expected.
//
//   Failure cases require a specific non-success status and leave the state
//   unchanged, since a partially advanced state that looks plausible is the
//   dangerous outcome.
//
// Reproduce
//   ctest -R verification_phase_space_evolution
//
// Build modes
//   Both. Failures are reported through an exit code, not assert().

#include "../../../src/internal/cqmom_case.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace qmomkit
{

struct TestCqmomPhaseSpaceEvolutionOptions : CqmomPhaseSpaceEvolutionOptions
{
    ConditionalMomentOptions conditional_options;
    QmomOptions qmom_options;
    bool adaptive_branch_counts = false;
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::none;
    VandermondeSolver vandermonde_solver = VandermondeSolver::automatic;
};

struct TestCqmomPhaseSpaceEvolutionResult
{
    CqmomMomentData moments;
    std::optional<CqmomQuadrature> quadrature;
    PhaseSpaceEvolutionStatus status = PhaseSpaceEvolutionStatus::success;
    PhaseSpaceEvolutionDiagnostics diagnostics;
    double final_time = 0.0;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == PhaseSpaceEvolutionStatus::success;
    }
};

// Advances a state built directly from a CqmomCase made of these moments and
// the case options in `options`. The state has no accepted quadrature, so
// phase-space evolution performs the initial CQMOM inversion itself.
TestCqmomPhaseSpaceEvolutionResult advance_cqmom_phase_space(const CqmomMomentData& moments, const double t, const double time_step,
                                                             const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                             const TestCqmomPhaseSpaceEvolutionOptions& options)
{
    const CqmomMomentLayout& layout = moments.layout;
    CqmomCase cqmom_case{moments,
                         options.conditional_options,
                         options.qmom_options,
                         options.adaptive_branch_counts,
                         options.diagnostic_level,
                         options.vandermonde_solver};
    // CQMOM requires the QMOM adaptive option to agree with the branch-count flag.
    cqmom_case.qmom_options.adaptive = cqmom_case.adaptive_branch_counts;
    CqmomState state = internal::CqmomConstructionAccess::create_state(cqmom_case);
    const CqmomPhaseSpaceEvolutionResult result =
        advance_cqmom_phase_space(state, t, time_step, velocity, source, static_cast<const CqmomPhaseSpaceEvolutionOptions&>(options));
    const CqmomState& result_state = result.ok() ? result.state : state;
    return TestCqmomPhaseSpaceEvolutionResult{CqmomMomentData{layout, result_state.moment_values()},
                                              result.ok() ? std::optional<CqmomQuadrature>{result.state.quadrature()} : result.diagnostic_quadrature,
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
    qmomkit::PhaseSpaceEvolutionScheme scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    double dt = 0.0;
    std::size_t steps = 0;
    double final_time = 0.0;
    double error_norm = 0.0;
    double observed_order = std::numeric_limits<double>::quiet_NaN();
    ErrorMetrics errors;
    std::size_t retry_count = 0;
    bool success = false;
};

struct VerificationResult
{
    bool passed = true;
    std::string message;
    std::vector<RunSummary> runs;
};

struct QmomStepRun
{
    qmomkit::QmomPhaseSpaceEvolutionResult result;
    bool diagnostics_successful = true;
    std::size_t retry_count = 0;
};

qmomkit::QmomState qmom_state_from_moments(const std::vector<double>& moments, const std::size_t node_count, const qmomkit::QmomOptions& options)
{
    return qmomkit::QmomCase{qmomkit::MomentSet(moments), node_count, options}.initial_state();
}

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

bool values_are_finite(const std::vector<double>& values)
{
    return std::all_of(values.begin(), values.end(), [](const double value) { return std::isfinite(value); });
}

bool qmom_quadrature_values_are_finite(const qmomkit::QmomResult& quadrature)
{
    return values_are_finite(quadrature.weights) && values_are_finite(quadrature.abscissas);
}

bool cqmom_quadrature_values_are_finite(const qmomkit::CqmomQuadrature& quadrature)
{
    if (!values_are_finite(quadrature.weights()))
    {
        return false;
    }
    for (const std::vector<double>& dimension : quadrature.abscissas())
    {
        if (!values_are_finite(dimension))
        {
            return false;
        }
    }
    return true;
}

bool cqmom_quadrature_topology_is_consistent(const qmomkit::CqmomQuadrature& quadrature)
{
    const std::size_t active_nodes = quadrature.weights().size();
    if (quadrature.total_nodes() != active_nodes || quadrature.abscissas().size() != quadrature.number_of_dimensions())
    {
        return false;
    }
    for (const std::vector<double>& dimension : quadrature.abscissas())
    {
        if (dimension.size() != active_nodes)
        {
            return false;
        }
    }
    return quadrature.terminal_paths().empty() || quadrature.terminal_paths().size() == active_nodes;
}

bool active_branch_counts_are_valid(const std::vector<qmomkit::CqmomActiveBranchCount>& counts)
{
    for (const qmomkit::CqmomActiveBranchCount& count : counts)
    {
        if (count.active_node_count == 0 || count.active_node_count > count.configured_node_count ||
            count.parent_path.size() != count.conditioning_level)
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

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::QmomCase>& configured)
{
    std::cerr << configured.error() << '\n';
}

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::CqmomCase>& configured)
{
    std::cerr << configured.error() << '\n';
}

void print_configuration_errors(const qmomkit::config::ConfigurationResult<qmomkit::CqmomPreparedCase>& configured)
{
    std::cerr << configured.error() << '\n';
}

qmomkit::QmomPhaseSpaceEvolutionOptions qmom_options(const double rk_relative_tolerance, const double rk_absolute_tolerance, const double max_step,
                                                     const qmomkit::PhaseSpaceEvolutionScheme scheme = qmomkit::PhaseSpaceEvolutionScheme::strang)
{
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    options.node_options.adaptive_options.relative_tolerance = rk_relative_tolerance;
    options.node_options.adaptive_options.absolute_tolerance = rk_absolute_tolerance;
    options.node_options.adaptive_options.initial_step = max_step;
    options.node_options.adaptive_options.max_step = max_step;
    options.source_options.ssp_rk3_options.estimate_error = true;
    options.scheme = scheme;
    return options;
}

const char* phase_space_evolution_scheme_name(const qmomkit::PhaseSpaceEvolutionScheme scheme)
{
    switch (scheme)
    {
    case qmomkit::PhaseSpaceEvolutionScheme::strang:
        return "strang";
    case qmomkit::PhaseSpaceEvolutionScheme::advsrc:
        return "advsrc";
    case qmomkit::PhaseSpaceEvolutionScheme::srcadv:
        return "srcadv";
    }
    return "unknown";
}

qmomkit::TestCqmomPhaseSpaceEvolutionOptions cqmom_options(const qmomkit::CqmomCase& configured_cqmom_case)
{
    qmomkit::TestCqmomPhaseSpaceEvolutionOptions options;
    options.conditional_options = configured_cqmom_case.conditional_options;
    options.qmom_options = configured_cqmom_case.qmom_options;
    options.adaptive_branch_counts = configured_cqmom_case.adaptive_branch_counts;
    options.diagnostic_level = configured_cqmom_case.diagnostic_level;
    options.vandermonde_solver = configured_cqmom_case.vandermonde_solver;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    options.node_options.adaptive_options.relative_tolerance = 1.0e-12;
    options.node_options.adaptive_options.absolute_tolerance = 1.0e-14;
    options.node_options.adaptive_options.initial_step = 0.01;
    options.node_options.adaptive_options.max_step = 0.01;
    options.source_options.ssp_rk3_options.estimate_error = true;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    return options;
}

qmomkit::CqmomPhaseSpaceEvolutionOptions cqmom_state_options()
{
    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    options.node_options.adaptive_options.relative_tolerance = 1.0e-12;
    options.node_options.adaptive_options.absolute_tolerance = 1.0e-14;
    options.node_options.adaptive_options.initial_step = 0.01;
    options.node_options.adaptive_options.max_step = 0.01;
    options.source_options.ssp_rk3_options.estimate_error = true;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    return options;
}

bool diagnostics_are_successful(const qmomkit::PhaseSpaceEvolutionDiagnostics& diagnostics)
{
    bool advection_success = false;
    bool reconstruction_success = false;
    switch (diagnostics.scheme)
    {
    case qmomkit::PhaseSpaceEvolutionScheme::strang:
        advection_success = diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::success &&
                            diagnostics.second_advection_status == qmomkit::NodeAdvectionStatus::success;
        reconstruction_success = diagnostics.first_reconstruction_diagnostics.success && diagnostics.final_reconstruction_diagnostics.success;
        break;
    case qmomkit::PhaseSpaceEvolutionScheme::advsrc:
        advection_success = diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::success;
        reconstruction_success = diagnostics.first_reconstruction_diagnostics.success;
        break;
    case qmomkit::PhaseSpaceEvolutionScheme::srcadv:
        advection_success = diagnostics.second_advection_status == qmomkit::NodeAdvectionStatus::success;
        reconstruction_success = diagnostics.final_reconstruction_diagnostics.success;
        break;
    }
    return diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::none &&
           diagnostics.post_source_qmom_status == qmomkit::QmomStatus::success && advection_success &&
           diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success && diagnostics.moment_source_diagnostics.accepted &&
           !diagnostics.moment_source_diagnostics.rejected && diagnostics.moment_source_retries == 0 && reconstruction_success;
}

bool cqmom_diagnostics_are_successful(const qmomkit::PhaseSpaceEvolutionDiagnostics& diagnostics)
{
    return diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::none && diagnostics.initial_cqmom_status == qmomkit::CqmomStatus::success &&
           diagnostics.post_source_cqmom_status == qmomkit::CqmomStatus::success &&
           diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::success &&
           diagnostics.second_advection_status == qmomkit::NodeAdvectionStatus::success &&
           diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success && diagnostics.moment_source_diagnostics.accepted &&
           !diagnostics.moment_source_diagnostics.rejected && diagnostics.moment_source_retries == 0 &&
           diagnostics.first_reconstruction_diagnostics.success && diagnostics.final_reconstruction_diagnostics.success;
}

bool cqmom_split_diagnostics_are_successful(const qmomkit::PhaseSpaceEvolutionDiagnostics& diagnostics,
                                            const qmomkit::PhaseSpaceEvolutionScheme scheme)
{
    bool advection_success = false;
    bool reconstruction_success = false;
    switch (scheme)
    {
    case qmomkit::PhaseSpaceEvolutionScheme::strang:
        advection_success = diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::success &&
                            diagnostics.second_advection_status == qmomkit::NodeAdvectionStatus::success;
        reconstruction_success = diagnostics.first_reconstruction_diagnostics.success && diagnostics.final_reconstruction_diagnostics.success;
        break;
    case qmomkit::PhaseSpaceEvolutionScheme::advsrc:
        advection_success = diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::success;
        reconstruction_success = diagnostics.first_reconstruction_diagnostics.success;
        break;
    case qmomkit::PhaseSpaceEvolutionScheme::srcadv:
        advection_success = diagnostics.second_advection_status == qmomkit::NodeAdvectionStatus::success;
        reconstruction_success = diagnostics.final_reconstruction_diagnostics.success;
        break;
    }
    return diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::none && diagnostics.initial_cqmom_status == qmomkit::CqmomStatus::success &&
           diagnostics.post_source_cqmom_status == qmomkit::CqmomStatus::success && advection_success &&
           diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success && diagnostics.moment_source_diagnostics.accepted &&
           !diagnostics.moment_source_diagnostics.rejected && diagnostics.moment_source_retries == 0 && reconstruction_success;
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

void verify_helper_based_common_cqmom_policies(const std::string& case_name, const qmomkit::PhaseSpaceEvolutionDiagnostics& diagnostics,
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

void commuting_velocity(double, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
{
    if (dimensions != 1)
    {
        throw std::invalid_argument("This verification velocity model is one-dimensional.");
    }
    for (std::size_t entry = 0; entry < coordinates.size(); ++entry)
    {
        velocities[entry] = commuting_a * coordinates[entry];
    }
}

qmomkit::QmomMomentRate commuting_source(double, const qmomkit::QmomSourceContext& context)
{
    std::vector<double> rhs(context.size(), 0.0);
    for (std::size_t order = 1; order < context.size(); ++order)
    {
        rhs[order] = -commuting_lambda_scale * static_cast<double>(order) * context.value(order);
    }
    return qmomkit::QmomMomentRate{std::move(rhs)};
}

qmomkit::QmomMomentRate zero_qmom_source(double, const qmomkit::QmomSourceContext& context)
{
    return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)};
}

qmomkit::QmomMomentRate qmom_boundary_injection_removal(double, const qmomkit::QmomSourceContext& moments)
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
    return qmomkit::QmomMomentRate{rate};
}

std::vector<double> qmom_boundary_rate_values(const std::size_t count)
{
    // The context is a non-owning view, so the moment vector must outlive the call.
    const std::vector<double> moments(count, 0.0);
    return qmom_boundary_injection_removal(0.0, qmomkit::QmomSourceContext{0.0, moments}).values;
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

qmomkit::QmomMomentRate qmom_constant_source(double, const qmomkit::QmomSourceContext& context)
{
    return qmomkit::QmomMomentRate{qmom_constant_source_values(context.size())};
}

std::vector<double> noncommuting_exact(const std::vector<double>& initial, const double elapsed)
{
    // Noncommuting QMOM reference: constant node velocity couples M_k to
    // M_{k-1}, while the source damps each order diagonally. The closed-form
    // triangular exponential below is independent of the phase-space evolution splitting
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

void noncommuting_velocity(double, const std::vector<double>&, std::vector<double>& velocities, const std::size_t dimensions)
{
    if (dimensions != 1)
    {
        throw std::invalid_argument("This verification velocity model is one-dimensional.");
    }
    std::fill(velocities.begin(), velocities.end(), noncommuting_b);
}

qmomkit::QmomMomentRate noncommuting_source(double, const qmomkit::QmomSourceContext& context)
{
    std::vector<double> rhs(context.size(), 0.0);
    for (std::size_t order = 1; order < context.size(); ++order)
    {
        rhs[order] = -noncommuting_gamma * static_cast<double>(order) * context.value(order);
    }
    return qmomkit::QmomMomentRate{std::move(rhs)};
}

QmomStepRun run_qmom_steps(const std::vector<double>& initial, const std::size_t node_count, const double t0, const double dt,
                           const std::size_t steps, const qmomkit::PhaseSpaceVelocity& velocity, const qmomkit::QmomMomentSource& source,
                           const qmomkit::QmomOptions& reconstruction_options, const qmomkit::QmomPhaseSpaceEvolutionOptions& options)
{
    qmomkit::QmomState state = qmom_state_from_moments(initial, node_count, reconstruction_options);
    double time = t0;
    std::optional<qmomkit::QmomPhaseSpaceEvolutionResult> result;
    bool diagnostics_successful = true;
    std::size_t retry_count = 0;
    if (steps == 0)
    {
        throw std::invalid_argument("QMOM phase-space evolution verification requires at least one step.");
    }
    for (std::size_t step = 0; step < steps; ++step)
    {
        result = qmomkit::advance_qmom_phase_space(state, time, dt, velocity, source, options);
        if (!result->ok())
        {
            retry_count += result->diagnostics.moment_source_retries;
            return {std::move(*result), diagnostics_successful, retry_count};
        }
        diagnostics_successful = diagnostics_successful && diagnostics_are_successful(result->diagnostics);
        retry_count += result->diagnostics.moment_source_retries;
        state = result->state;
        time = result->final_time;
    }
    return {std::move(*result), diagnostics_successful, retry_count};
}

RunSummary validate_qmom_run(VerificationResult& verification, const std::vector<double>& initial, const std::size_t node_count, const double t0,
                             const double dt, const std::size_t steps, const qmomkit::PhaseSpaceVelocity& velocity,
                             const qmomkit::QmomMomentSource& source, const qmomkit::QmomPhaseSpaceEvolutionOptions& options,
                             const qmomkit::QmomOptions& reconstruction_options, const std::vector<double>& exact, const double moment_tolerance,
                             const double reconstruction_tolerance)
{
    const QmomStepRun step_run = run_qmom_steps(initial, node_count, t0, dt, steps, velocity, source, reconstruction_options, options);
    const qmomkit::QmomPhaseSpaceEvolutionResult& result = step_run.result;
    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = dt;
    summary.steps = steps;
    summary.retry_count = step_run.retry_count;
    if (!result.ok())
    {
        summary.final_time = result.final_time;
        append_failure(verification, "phase-space evolution integrator failed: " + result.message);
        return summary;
    }
    summary.final_time = result.final_time;
    summary.success = true;
    if (!step_run.diagnostics_successful)
    {
        summary.success = false;
        append_failure(verification, "phase-space evolution diagnostics did not report successful completion of every substage.");
    }
    if (std::abs(result.final_time - (t0 + dt * static_cast<double>(steps))) > 1.0e-12)
    {
        summary.success = false;
        append_failure(verification, "phase-space evolution final time does not match the requested final time.");
    }

    update_errors(summary.errors, result.state.moment_values(), exact);
    summary.error_norm = scaled_infinity_error(result.state.moment_values(), exact);
    const qmomkit::Quadrature1D final_quadrature{result.state.quadrature().abscissas, result.state.quadrature().weights,
                                                 qmomkit::QuadratureStatus::success, ""};
    const std::vector<double> reconstructed = qmomkit::reconstruct_raw_moments(final_quadrature, exact.size());
    update_reconstruction_errors(summary.errors, reconstructed, exact);
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

void validate_commuting_component_checks(VerificationResult& verification, const qmomkit::QmomCase& qmom_case)
{
    const qmomkit::QmomResult initial = qmomkit::compute_qmom(qmom_case);
    if (!initial.ok())
    {
        append_failure(verification, "initial QMOM inversion failed during component checks.");
        return;
    }

    qmomkit::NodeAdvectionOptions node_options;
    node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    node_options.adaptive_options.relative_tolerance = 1.0e-12;
    node_options.adaptive_options.absolute_tolerance = 1.0e-14;
    node_options.adaptive_options.initial_step = 0.005;
    node_options.adaptive_options.max_step = 0.005;
    constexpr double elapsed = 0.03;
    const qmomkit::QmomNodeAdvectionResult advanced = qmomkit::advance_nodes(initial, 0.0, elapsed, commuting_velocity, node_options);
    if (!advanced.ok())
    {
        append_failure(verification, "advection operator component check failed.");
        return;
    }
    const double scale = std::exp(commuting_a * elapsed);
    for (std::size_t node = 0; node < initial.weights.size(); ++node)
    {
        if (std::abs(advanced.quadrature.weights[node] - initial.weights[node]) > 1.0e-13 ||
            std::abs(advanced.quadrature.abscissas[node] - scale * initial.abscissas[node]) > 1.0e-10)
        {
            append_failure(verification, "advection operator did not preserve weights while scaling abscissas.");
            break;
        }
    }

    constexpr double source_probe_step = 0.002;
    const qmomkit::QmomMomentSourceUpdateResult sourced =
        qmomkit::advance_qmom_moments_with_source(qmom_case.moments.values(), 0.0, source_probe_step, commuting_source, {});
    if (!sourced.ok())
    {
        append_failure(verification, "moment-source operator component check failed.");
        return;
    }
    const std::vector<double> exact = commuting_source_exact(qmom_case.moments.values(), source_probe_step);
    if (scaled_infinity_error(sourced.moments, exact) > 2.0e-9)
    {
        append_failure(verification, "moment-source operator diagonal-law check exceeds tolerance.");
    }
}

VerificationResult run_commuting_qmom(const std::filesystem::path& configuration_path)
{
    // Verify QMOM phase-space evolution when phase-space advection and moment-source evolution
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
    const auto options = qmom_options(1.0e-12, 1.0e-14, 0.01);

    verification.runs.push_back(validate_qmom_run(verification, initial, qmom_case.node_count, t0, 0.08, 1, commuting_velocity, commuting_source,
                                                  options, qmom_case.qmom_options, commuting_exact(initial, 0.08), 2.0e-6, 2.0e-6));
    verification.runs.push_back(validate_qmom_run(verification, initial, qmom_case.node_count, t0, 0.04, 6, commuting_velocity, commuting_source,
                                                  options, qmom_case.qmom_options, commuting_exact(initial, 0.24), 2.0e-6, 2.0e-6));

    const std::vector<double> refinement_steps{0.12, 0.06, 0.03};
    double previous_error = 0.0;
    for (const double dt : refinement_steps)
    {
        const std::size_t steps = static_cast<std::size_t>(std::lround(0.24 / dt));
        RunSummary summary = validate_qmom_run(verification, initial, qmom_case.node_count, t0, dt, steps, commuting_velocity, commuting_source,
                                               options, qmom_case.qmom_options, commuting_exact(initial, 0.24), 8.0e-6, 8.0e-6);
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
    const qmomkit::QmomState initial_state = qmom_state_from_moments(initial, qmom_case.node_count, qmom_case.qmom_options);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    qmomkit::QmomPhaseSpaceEvolutionOptions options = qmom_options(1.0e-12, 1.0e-14, 0.01, qmomkit::PhaseSpaceEvolutionScheme::srcadv);

    constexpr double boundary_dt = 0.04;
    const qmomkit::QmomPhaseSpaceEvolutionResult boundary_result =
        qmomkit::advance_qmom_phase_space(initial_state, 0.0, boundary_dt, zero_velocity, zero_qmom_source, qmom_boundary_injection_removal, options);
    RunSummary boundary_summary;
    boundary_summary.scheme = options.scheme;
    boundary_summary.dt = boundary_dt;
    boundary_summary.steps = 1;
    boundary_summary.final_time = boundary_result.final_time;
    boundary_summary.retry_count = boundary_result.diagnostics.moment_source_retries;
    boundary_summary.success = boundary_result.ok() && diagnostics_are_successful(boundary_result.diagnostics);
    std::vector<double> boundary_exact = initial;
    const std::vector<double> boundary_rate = qmom_boundary_rate_values(initial.size());
    for (std::size_t i = 0; i < boundary_exact.size(); ++i)
    {
        boundary_exact[i] += boundary_dt * boundary_rate[i];
    }
    if (!boundary_result.ok())
    {
        append_failure(verification, "QMOM boundary injection/removal phase-space evolution failed: " + boundary_result.message);
    }
    else
    {
        update_errors(boundary_summary.errors, boundary_result.state.moment_values(), boundary_exact);
        const qmomkit::Quadrature1D final_quadrature{boundary_result.state.quadrature().abscissas, boundary_result.state.quadrature().weights,
                                                     qmomkit::QuadratureStatus::success, ""};
        update_reconstruction_errors(boundary_summary.errors, qmomkit::reconstruct_raw_moments(final_quadrature, boundary_exact.size()),
                                     boundary_exact);
        boundary_summary.error_norm = scaled_infinity_error(boundary_result.state.moment_values(), boundary_exact);
        if (boundary_summary.errors.moment_scaled > 5.0e-11 || boundary_summary.errors.reconstruction_scaled > 5.0e-11)
        {
            boundary_summary.success = false;
            append_failure(verification, "QMOM boundary injection/removal moments exceeded the analytical tolerance.");
        }
    }
    verification.runs.push_back(boundary_summary);

    constexpr double combined_dt = 0.03;
    const qmomkit::QmomPhaseSpaceEvolutionResult combined_result = qmomkit::advance_qmom_phase_space(
        initial_state, 0.0, combined_dt, zero_velocity, qmom_constant_source, qmom_boundary_injection_removal, options);
    RunSummary combined_summary;
    combined_summary.scheme = options.scheme;
    combined_summary.dt = combined_dt;
    combined_summary.steps = 1;
    combined_summary.final_time = combined_result.final_time;
    combined_summary.retry_count = combined_result.diagnostics.moment_source_retries;
    combined_summary.success = combined_result.ok() && diagnostics_are_successful(combined_result.diagnostics);
    std::vector<double> combined_exact = initial;
    const std::vector<double> source_rate = qmom_constant_source_values(initial.size());
    for (std::size_t i = 0; i < combined_exact.size(); ++i)
    {
        combined_exact[i] += combined_dt * (source_rate[i] + boundary_rate[i]);
    }
    if (!combined_result.ok())
    {
        append_failure(verification, "QMOM combined source plus boundary contribution phase-space evolution failed: " + combined_result.message);
    }
    else
    {
        update_errors(combined_summary.errors, combined_result.state.moment_values(), combined_exact);
        const qmomkit::Quadrature1D final_quadrature{combined_result.state.quadrature().abscissas, combined_result.state.quadrature().weights,
                                                     qmomkit::QuadratureStatus::success, ""};
        update_reconstruction_errors(combined_summary.errors, qmomkit::reconstruct_raw_moments(final_quadrature, combined_exact.size()),
                                     combined_exact);
        combined_summary.error_norm = scaled_infinity_error(combined_result.state.moment_values(), combined_exact);
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
        qmomkit::QmomMomentBoundaryContribution boundary;
        std::string expected_message_fragment;
    };
    const std::vector<BoundaryFailureCase> failure_cases{{"wrong size", [](double, const qmomkit::QmomSourceContext& moments)
                                                          { return qmomkit::QmomMomentRate{std::vector<double>(moments.size() + 1, 0.0)}; },
                                                          "boundary contribution"},
                                                         {"NaN",
                                                          [](double, const qmomkit::QmomSourceContext& moments)
                                                          {
                                                              std::vector<double> rate(moments.size(), 0.0);
                                                              rate.front() = std::numeric_limits<double>::quiet_NaN();
                                                              return qmomkit::QmomMomentRate{rate};
                                                          },
                                                          "boundary contribution"},
                                                         {"infinity",
                                                          [](double, const qmomkit::QmomSourceContext& moments)
                                                          {
                                                              std::vector<double> rate(moments.size(), 0.0);
                                                              rate.front() = std::numeric_limits<double>::infinity();
                                                              return qmomkit::QmomMomentRate{rate};
                                                          },
                                                          "boundary contribution"}};

    for (const BoundaryFailureCase& failure_case : failure_cases)
    {
        const qmomkit::QmomPhaseSpaceEvolutionResult failure_result =
            qmomkit::advance_qmom_phase_space(initial_state, 0.0, 0.01, zero_velocity, zero_qmom_source, failure_case.boundary, options);
        RunSummary failure_summary;
        failure_summary.scheme = options.scheme;
        failure_summary.dt = 0.01;
        failure_summary.steps = 1;
        failure_summary.final_time = failure_result.final_time;
        failure_summary.retry_count = failure_result.diagnostics.moment_source_retries;
        failure_summary.success = !failure_result.ok() && failure_result.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed &&
                                  failure_result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step &&
                                  failure_result.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::source_failed &&
                                  contains_text(failure_result.message, failure_case.expected_message_fragment);
        if (!failure_summary.success)
        {
            append_failure(verification,
                           "QMOM boundary failure case '" + failure_case.name + "' did not report a boundary failure in the moment-source stage.");
        }
        verification.runs.push_back(failure_summary);
    }

    const qmomkit::QmomMomentBoundaryContribution nonrealizable_boundary = [](double, const qmomkit::QmomSourceContext& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        rate.front() = -100.0;
        return qmomkit::QmomMomentRate{rate};
    };
    const qmomkit::QmomPhaseSpaceEvolutionResult nonrealizable_result =
        qmomkit::advance_qmom_phase_space(initial_state, 0.0, 0.02, zero_velocity, zero_qmom_source, nonrealizable_boundary, options);
    RunSummary nonrealizable_summary;
    nonrealizable_summary.scheme = options.scheme;
    nonrealizable_summary.dt = 0.02;
    nonrealizable_summary.steps = 1;
    nonrealizable_summary.final_time = nonrealizable_result.final_time;
    nonrealizable_summary.retry_count = nonrealizable_result.diagnostics.moment_source_retries;
    nonrealizable_summary.success = !nonrealizable_result.ok() &&
                                    nonrealizable_result.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed &&
                                    nonrealizable_result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::post_source_inversion &&
                                    nonrealizable_result.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success &&
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
    // advsrc and srcadv compositions should approach first order and differ at coarse
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

    const std::vector<qmomkit::PhaseSpaceEvolutionScheme> schemes{
        qmomkit::PhaseSpaceEvolutionScheme::strang, qmomkit::PhaseSpaceEvolutionScheme::advsrc, qmomkit::PhaseSpaceEvolutionScheme::srcadv};
    for (const qmomkit::PhaseSpaceEvolutionScheme scheme : schemes)
    {
        const auto options = qmom_options(1.0e-12, 1.0e-14, 0.01, scheme);
        const std::size_t first_run_index = verification.runs.size();
        double previous_error = 0.0;
        for (const double dt : refinement_steps)
        {
            const std::size_t steps = static_cast<std::size_t>(std::lround(final_time / dt));
            const double tolerance = (scheme == qmomkit::PhaseSpaceEvolutionScheme::strang) ? 8.0e-4 : 3.0e-2;
            RunSummary summary = validate_qmom_run(verification, initial, qmom_case.node_count, t0, dt, steps, noncommuting_velocity,
                                                   noncommuting_source, options, qmom_case.qmom_options, exact, tolerance, tolerance);
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
            if (scheme == qmomkit::PhaseSpaceEvolutionScheme::strang)
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

    const auto advsrc_options = qmom_options(1.0e-12, 1.0e-14, 0.01, qmomkit::PhaseSpaceEvolutionScheme::advsrc);
    const auto srcadv_options = qmom_options(1.0e-12, 1.0e-14, 0.01, qmomkit::PhaseSpaceEvolutionScheme::srcadv);
    const QmomStepRun cd_coarse = run_qmom_steps(initial, qmom_case.node_count, t0, refinement_steps.front(),
                                                 static_cast<std::size_t>(std::lround(final_time / refinement_steps.front())), noncommuting_velocity,
                                                 noncommuting_source, qmom_case.qmom_options, advsrc_options);
    const QmomStepRun dc_coarse = run_qmom_steps(initial, qmom_case.node_count, t0, refinement_steps.front(),
                                                 static_cast<std::size_t>(std::lround(final_time / refinement_steps.front())), noncommuting_velocity,
                                                 noncommuting_source, qmom_case.qmom_options, srcadv_options);
    if (cd_coarse.result.ok() && dc_coarse.result.ok() &&
        scaled_infinity_error(cd_coarse.result.state.moment_values(), dc_coarse.result.state.moment_values()) < 1.0e-4)
    {
        append_failure(verification, "noncommuting advsrc and srcadv schemes did not differ at finite timestep.");
    }
    return verification;
}

void verify_qmom_advection_failure_result(const std::string& case_name, const qmomkit::QmomPhaseSpaceEvolutionResult& result,
                                          const std::vector<double>& initial_moments, const std::string& expected_message_fragment,
                                          VerificationResult& verification)
{
    // A velocity callback contract violation should stop in the advection
    // node-advection stage. The returned state is the accepted state from before
    // the failed substep, not a completed final state.
    if (result.ok())
    {
        append_failure(verification, case_name + " unexpectedly reported QMOM phase-space evolution success.");
    }
    if (result.status != qmomkit::PhaseSpaceEvolutionStatus::advection_failed)
    {
        append_failure(verification, case_name + " did not report PhaseSpaceEvolutionStatus::advection_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::PhaseSpaceEvolutionStage::first_advection_half_step)
    {
        append_failure(verification, case_name + " did not identify the first advection half-step as the failed stage.");
    }
    if (result.diagnostics.first_advection_status != qmomkit::NodeAdvectionStatus::velocity_failed)
    {
        append_failure(verification, case_name + " did not report NodeAdvectionStatus::velocity_failed.");
    }
    if (!contains_text(result.message, expected_message_fragment) || !contains_text(result.diagnostics.message, expected_message_fragment) ||
        !contains_text(result.diagnostics.first_advection_diagnostics.message, expected_message_fragment))
    {
        append_failure(verification, case_name + " did not preserve a useful velocity-callback failure message.");
    }
    if (std::abs(result.final_time) > 1.0e-14)
    {
        append_failure(verification, case_name + " advanced final_time despite phase-space advection failure.");
    }
    if (!vectors_equal(result.state.moment_values(), initial_moments, 0.0))
    {
        append_failure(verification, case_name + " returned moments changed by a failed advection substep.");
    }
    if (result.diagnostics.moment_source_diagnostics.accepted || result.diagnostics.moment_source_diagnostics.rejected ||
        result.diagnostics.first_reconstruction_diagnostics.success || result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, case_name + " reported later phase-space evolution stage diagnostics after the failed advection substep.");
    }
}

void verify_qmom_rk23_minimum_step_acceptance_result(const qmomkit::QmomPhaseSpaceEvolutionResult& result, const std::vector<double>& exact,
                                                     VerificationResult& verification)
{
    // The velocity callback is finite and dimensionally valid. Adaptive RK23
    // rejects the first requested step, then advances at the effective minimum
    // step instead of failing because the controller would prefer a smaller one.
    if (!result.ok())
    {
        append_failure(verification, "QMOM RK23 minimum-step case did not report phase-space evolution success.");
    }
    if (result.status != qmomkit::PhaseSpaceEvolutionStatus::success)
    {
        append_failure(verification, "QMOM RK23 minimum-step case did not report PhaseSpaceEvolutionStatus::success.");
    }
    if (result.diagnostics.failed_stage != qmomkit::PhaseSpaceEvolutionStage::none)
    {
        append_failure(verification, "QMOM RK23 minimum-step case reported a failed phase-space evolution stage.");
    }
    if (result.diagnostics.first_advection_status != qmomkit::NodeAdvectionStatus::success ||
        result.diagnostics.second_advection_status != qmomkit::NodeAdvectionStatus::success ||
        result.diagnostics.first_advection_diagnostics.adaptive_status != qmomkit::AdaptiveRkStatus::success ||
        result.diagnostics.second_advection_diagnostics.adaptive_status != qmomkit::AdaptiveRkStatus::success)
    {
        append_failure(verification, "QMOM RK23 minimum-step case did not preserve successful RK23 diagnostics.");
    }
    if (result.diagnostics.first_advection_diagnostics.rejected_step_count == 0 || result.diagnostics.total_node_advection_rejected_steps == 0 ||
        result.diagnostics.first_advection_diagnostics.step_count == 0 || result.diagnostics.second_advection_diagnostics.step_count == 0)
    {
        append_failure(verification, "QMOM RK23 minimum-step case did not preserve accepted and rejected RK23 diagnostics.");
    }
    if (std::abs(result.final_time - 0.08) > 1.0e-14 || std::abs(result.diagnostics.first_advection_diagnostics.final_time - 0.04) > 1.0e-14 ||
        std::abs(result.diagnostics.second_advection_diagnostics.final_time - 0.08) > 1.0e-14)
    {
        append_failure(verification, "QMOM RK23 minimum-step case did not reach the expected endpoint.");
    }
    if (!result.diagnostics.moment_source_diagnostics.accepted || !result.diagnostics.first_reconstruction_diagnostics.success ||
        !result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, "QMOM RK23 minimum-step case did not complete the later phase-space evolution stages.");
    }
    if (!values_are_finite(result.state.moment_values()))
    {
        append_failure(verification, "QMOM RK23 minimum-step case returned a non-finite final moment.");
    }
    if (!qmom_quadrature_values_are_finite(result.state.quadrature()))
    {
        append_failure(verification, "QMOM RK23 minimum-step case returned a non-finite final quadrature value.");
    }
    if (result.state.moment_values().size() != exact.size())
    {
        append_failure(verification, "QMOM RK23 minimum-step case returned an unexpected final moment count.");
    }
    const qmomkit::Quadrature1D final_quadrature{result.state.quadrature().abscissas, result.state.quadrature().weights,
                                                 qmomkit::QuadratureStatus::success, ""};
    const std::vector<double> reconstructed = qmomkit::reconstruct_raw_moments(final_quadrature, exact.size());
    if (result.state.moment_values().size() == exact.size() &&
        (scaled_infinity_error(result.state.moment_values(), exact) > 2.0e-6 || scaled_infinity_error(reconstructed, exact) > 2.0e-6))
    {
        append_failure(verification, "QMOM RK23 minimum-step case final state exceeded the analytical advection-reference tolerance.");
    }
}

void verify_qmom_post_source_inversion_failure_result(const qmomkit::QmomPhaseSpaceEvolutionResult& result, const double expected_failure_time,
                                                      VerificationResult& verification)
{
    // The source callback is finite, correctly sized, and successfully
    // integrated, but it drives m0 negative. The following QMOM inversion must
    // reject the evolved moments instead of returning a successful final
    // quadrature.
    if (result.ok())
    {
        append_failure(verification, "QMOM post-source inversion failure unexpectedly reported phase-space evolution success.");
    }
    if (result.status != qmomkit::PhaseSpaceEvolutionStatus::inversion_failed)
    {
        append_failure(verification, "QMOM post-source inversion failure did not report PhaseSpaceEvolutionStatus::inversion_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::PhaseSpaceEvolutionStage::post_source_inversion)
    {
        append_failure(verification, "QMOM post-source inversion failure did not identify post_source_inversion as the failed stage.");
    }
    if (result.diagnostics.moment_source_status != qmomkit::MomentSourceUpdateStatus::success ||
        result.diagnostics.moment_source_diagnostics.ssp_rk3_status != qmomkit::SspRk3Status::success ||
        !result.diagnostics.moment_source_diagnostics.accepted)
    {
        append_failure(verification, "QMOM post-source inversion failure did not preserve the successful source-step diagnostics.");
    }
    if (result.diagnostics.post_source_qmom_status == qmomkit::QmomStatus::success)
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
        append_failure(verification, "QMOM post-source inversion failure did not stop at the time of the moment-source stage.");
    }
    if (result.diagnostics.second_advection_diagnostics.step_count > 0 || result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, "QMOM post-source inversion failure reported later-stage completion after the failed inversion.");
    }
}

VerificationResult run_qmom_phase_space_evolution_failures(const std::filesystem::path& configuration_path)
{
    // Exercise QMOM phase-space evolution failure propagation after a valid initial
    // inversion: velocity callback errors and a moment-source stage that creates a
    // non-realizable post-source moment sequence.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_qmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load QMOM phase-space evolution failure configuration.");
        return verification;
    }

    const qmomkit::QmomCase& qmom_case = configured.value();
    const std::vector<double> initial = qmom_case.moments.values();
    const qmomkit::QmomState initial_state = qmom_state_from_moments(initial, qmom_case.node_count, qmom_case.qmom_options);
    const auto options = qmom_options(1.0e-12, 1.0e-14, 0.01);
    constexpr double dt = 0.08;

    struct FailureCase
    {
        std::string name;
        qmomkit::PhaseSpaceVelocity velocity;
        std::string expected_message_fragment;
    };

    const std::vector<FailureCase> failure_cases{{"wrong velocity dimension",
                                                  [](double, const std::vector<double>&, std::vector<double>& velocities, std::size_t)
                                                  {
                                                      // Simulates a user callback writing more velocity components than
                                                      // the configured one-dimensional phase space has.
                                                      velocities.assign(velocities.size() + 1, 0.0);
                                                  },
                                                  "Phase-space velocity callback returned the wrong dimension:"},
                                                 {"non-finite velocity",
                                                  [](double, const std::vector<double>&, std::vector<double>& velocities, std::size_t)
                                                  {
                                                      // Simulates a physical model producing an invalid node velocity.
                                                      std::fill(velocities.begin(), velocities.end(), std::numeric_limits<double>::quiet_NaN());
                                                  },
                                                  "Phase-space velocity callback returned a non-finite value:"}};

    for (const FailureCase& failure_case : failure_cases)
    {
        const qmomkit::QmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_qmom_phase_space(initial_state, 0.0, dt, failure_case.velocity, zero_qmom_source, options);
        verify_qmom_advection_failure_result(failure_case.name, result, initial_state.moment_values(), failure_case.expected_message_fragment,
                                             verification);

        RunSummary summary;
        summary.scheme = options.scheme;
        summary.dt = dt;
        summary.steps = 1;
        summary.final_time = result.final_time;
        summary.retry_count = result.diagnostics.moment_source_retries;
        summary.success = !result.ok() && result.status == qmomkit::PhaseSpaceEvolutionStatus::advection_failed &&
                          result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::first_advection_half_step &&
                          result.diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::velocity_failed &&
                          vectors_equal(result.state.moment_values(), initial_state.moment_values(), 0.0);
        verification.runs.push_back(summary);
    }

    qmomkit::QmomPhaseSpaceEvolutionOptions rk_minimum_options = options;
    rk_minimum_options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    rk_minimum_options.node_options.adaptive_options.initial_step = 0.04;
    rk_minimum_options.node_options.adaptive_options.max_step = 0.04;
    rk_minimum_options.node_options.adaptive_options.min_step = 0.02;
    const qmomkit::PhaseSpaceVelocity stiff_velocity = qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>& xi)
        {
            // This is a valid velocity callback, but the tight RK23 tolerances create
            // an oversized local-error estimate for the first half-step. With
            // min_step = dt/4, the next accepted update
            // advances at the configured minimum instead of failing.
            return std::vector<double>{xi[0]};
        });
    const qmomkit::QmomPhaseSpaceEvolutionResult rk_minimum =
        qmomkit::advance_qmom_phase_space(initial_state, 0.0, dt, stiff_velocity, zero_qmom_source, rk_minimum_options);
    std::vector<double> rk_minimum_exact(initial_state.moment_values().size(), 0.0);
    for (std::size_t order = 0; order < rk_minimum_exact.size(); ++order)
    {
        rk_minimum_exact[order] = initial_state.moment_values()[order] * std::exp(static_cast<double>(order) * dt);
    }
    verify_qmom_rk23_minimum_step_acceptance_result(rk_minimum, rk_minimum_exact, verification);

    RunSummary rk_summary;
    rk_summary.scheme = rk_minimum_options.scheme;
    rk_summary.dt = dt;
    rk_summary.steps = 1;
    rk_summary.final_time = rk_minimum.final_time;
    rk_summary.retry_count = rk_minimum.diagnostics.moment_source_retries;
    rk_summary.success = rk_minimum.ok() && rk_minimum.status == qmomkit::PhaseSpaceEvolutionStatus::success &&
                         rk_minimum.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::none &&
                         rk_minimum.diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::success &&
                         rk_minimum.diagnostics.second_advection_status == qmomkit::NodeAdvectionStatus::success &&
                         rk_minimum.diagnostics.total_node_advection_rejected_steps > 0 && rk_minimum.diagnostics.total_node_advection_steps > 0;
    verification.runs.push_back(rk_summary);

    const auto nonrealizable_source = [](double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> rates(context.size(), 0.0);
        if (!rates.empty())
        {
            rates.front() = -20.0;
        }
        return qmomkit::QmomMomentRate{std::move(rates)};
    };
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const qmomkit::QmomPhaseSpaceEvolutionResult post_source_inversion =
        qmomkit::advance_qmom_phase_space(initial_state, 0.0, dt, zero_velocity, nonrealizable_source, options);
    verify_qmom_post_source_inversion_failure_result(post_source_inversion, 0.0, verification);

    RunSummary reconstruction_summary;
    reconstruction_summary.scheme = options.scheme;
    reconstruction_summary.dt = dt;
    reconstruction_summary.steps = 1;
    reconstruction_summary.final_time = post_source_inversion.final_time;
    reconstruction_summary.retry_count = post_source_inversion.diagnostics.moment_source_retries;
    reconstruction_summary.success = !post_source_inversion.ok() &&
                                     post_source_inversion.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed &&
                                     post_source_inversion.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::post_source_inversion &&
                                     post_source_inversion.diagnostics.post_source_qmom_status != qmomkit::QmomStatus::success;
    verification.runs.push_back(reconstruction_summary);
    return verification;
}

VerificationResult run_qmom_phase_space_evolution_inversion_failures(const std::filesystem::path& configuration_path)
{
    // Construction of an accepted state owns initial inversion, so an invalid
    // configured moment sequence fails before a state can enter evolution.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_qmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load QMOM phase-space evolution inversion-failure configuration.");
        return verification;
    }

    const qmomkit::QmomCase& qmom_case = configured.value();
    bool initial_state_failed = false;
    try
    {
        (void)qmom_case.initial_state();
    }
    catch (const std::exception&)
    {
        initial_state_failed = true;
    }
    if (!initial_state_failed)
    {
        append_failure(verification, "QMOM invalid initial moments unexpectedly produced an accepted state.");
    }

    RunSummary summary;
    summary.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    summary.dt = 0.08;
    summary.steps = 1;
    summary.final_time = 0.0;
    summary.retry_count = 0;
    summary.success = initial_state_failed;
    verification.runs.push_back(summary);
    return verification;
}

std::vector<double> cqmom_exact_values(const qmomkit::CqmomMomentData& initial, const double elapsed)
{
    // CQMOM commuting reference: diagonal phase-space scaling and diagonal
    // source damping give each tensor moment M_ij an exact scalar exponential
    // factor based on its multi-index.
    std::vector<double> exact(initial.values.size(), 0.0);
    const std::vector<qmomkit::MultiIndex> indices = initial.layout.required_indices.to_multi_indices();
    for (std::size_t flat = 0; flat < initial.values.size(); ++flat)
    {
        const qmomkit::MultiIndex& index = indices[flat];
        const double lambda = cqmom_lambda_scale * static_cast<double>(index[0] + index[1]);
        const double exponent = static_cast<double>(index[0]) * cqmom_a + static_cast<double>(index[1]) * cqmom_c - lambda;
        exact[flat] = initial.values[flat] * std::exp(exponent * elapsed);
    }
    return exact;
}

std::vector<double> cqmom_advection_exact_values(const qmomkit::CqmomMomentData& initial, const double elapsed)
{
    std::vector<double> exact(initial.values.size(), 0.0);
    const std::vector<qmomkit::MultiIndex> indices = initial.layout.required_indices.to_multi_indices();
    for (std::size_t flat = 0; flat < initial.values.size(); ++flat)
    {
        const qmomkit::MultiIndex& index = indices[flat];
        const double exponent = static_cast<double>(index[0]) * cqmom_a + static_cast<double>(index[1]) * cqmom_c;
        exact[flat] = initial.values[flat] * std::exp(exponent * elapsed);
    }
    return exact;
}

std::vector<double> cqmom_noncommuting_exact_values(const qmomkit::CqmomQuadrature& initial, const std::vector<std::size_t>& extents,
                                                    const double elapsed)
{
    const double scale = std::exp(elapsed);
    std::vector<double> exact(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat = 0; flat < exact.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, extents);
        for (std::size_t node = 0; node < initial.total_nodes(); ++node)
        {
            const double x = scale * initial.abscissas()[0][node];
            const double y = initial.abscissas()[1][node] + (scale - 1.0) * initial.abscissas()[0][node];
            exact[flat] += initial.weights()[node] * integer_power(x, index[0]) * integer_power(y, index[1]);
        }
    }
    return exact;
}

void cqmom_velocity(double, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
{
    if (dimensions != 2)
    {
        throw std::invalid_argument("This verification velocity model is two-dimensional.");
    }
    for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
    {
        velocities[node * dimensions] = cqmom_a * coordinates[node * dimensions];
        velocities[node * dimensions + 1] = cqmom_c * coordinates[node * dimensions + 1];
    }
}

void cqmom_noncommuting_velocity(double, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
{
    if (dimensions != 2)
    {
        throw std::invalid_argument("This verification velocity model is two-dimensional.");
    }
    for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
    {
        velocities[node * dimensions] = coordinates[node * dimensions];
        velocities[node * dimensions + 1] = 0.0;
    }
}

qmomkit::CqmomMomentRate cqmom_source_for_context(double, const qmomkit::CqmomSourceContext& moments)
{
    std::vector<double> rhs(moments.size(), 0.0);
    for (std::size_t entry = 0; entry < moments.size(); ++entry)
    {
        const qmomkit::MultiIndexView index = moments.index(entry);
        const double lambda = cqmom_lambda_scale * static_cast<double>(index[0] + index[1]);
        rhs[entry] = -lambda * moments.value(entry);
    }
    return qmomkit::CqmomMomentRate{std::move(rhs)};
}

qmomkit::CqmomMomentRate cqmom_state_zero_source(double, const qmomkit::CqmomSourceContext& moments)
{
    return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)};
}

qmomkit::CqmomMomentSource make_cqmom_noncommuting_source(const double time_step)
{
    std::optional<std::vector<double>> base;
    std::optional<std::vector<double>> first_stage_rate;
    std::optional<std::vector<double>> second_stage_rate;
    std::size_t rhs_call = 0;

    return [time_step, base = std::move(base), first_stage_rate = std::move(first_stage_rate), second_stage_rate = std::move(second_stage_rate),
            rhs_call = std::move(rhs_call)](double, const qmomkit::CqmomSourceContext& context) mutable
    {
        const std::vector<std::size_t> extents{5, 4};
        const auto full_rates_from_full = [&extents](const std::vector<double>& full)
        {
            std::vector<double> rates(full.size(), 0.0);
            for (std::size_t flat = 0; flat < full.size(); ++flat)
            {
                const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, extents);
                if (index[1] == 0 || index[0] + 1 >= extents[0])
                {
                    continue;
                }
                const qmomkit::MultiIndex lower_index{index[0] + 1, index[1] - 1};
                rates[flat] = static_cast<double>(index[1]) * full[qmomkit::flatten_multi_index(lower_index, extents)];
            }
            return rates;
        };
        const auto required_rates_from_full = [&context, &extents](const std::vector<double>& full_rates)
        {
            std::vector<double> rates(context.size(), 0.0);
            for (std::size_t entry = 0; entry < context.size(); ++entry)
            {
                rates[entry] = full_rates[qmomkit::flatten_multi_index(context.index(entry), extents)];
            }
            return rates;
        };

        if (rhs_call == 0)
        {
            const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(context.hierarchy());
            const qmomkit::CqmomCase current_cqmom_case{qmomkit::CqmomMomentData{layout, context.values()},
                                                        {},
                                                        qmomkit::QmomOptions{false},
                                                        false,
                                                        qmomkit::CqmomDiagnosticLevel::none,
                                                        qmomkit::VandermondeSolver::automatic};
            const qmomkit::CqmomResult current = qmomkit::compute_cqmom(current_cqmom_case);
            if (!current.ok())
            {
                throw std::runtime_error("Manufactured CQMOM source could not reconstruct its current required moments.");
            }
            base = qmomkit::reconstruct_raw_moments(current.quadrature, extents).values;
            first_stage_rate = full_rates_from_full(*base);
            ++rhs_call;
            return qmomkit::CqmomMomentRate{required_rates_from_full(*first_stage_rate)};
        }

        if (rhs_call == 1)
        {
            std::vector<double> stage_one = *base;
            for (std::size_t entry = 0; entry < stage_one.size(); ++entry)
            {
                stage_one[entry] += time_step * (*first_stage_rate)[entry];
            }
            second_stage_rate = full_rates_from_full(stage_one);
            ++rhs_call;
            return qmomkit::CqmomMomentRate{required_rates_from_full(*second_stage_rate)};
        }

        std::vector<double> stage_two = *base;
        for (std::size_t entry = 0; entry < stage_two.size(); ++entry)
        {
            const double stage_one = (*base)[entry] + time_step * (*first_stage_rate)[entry];
            stage_two[entry] = 0.75 * (*base)[entry] + 0.25 * (stage_one + time_step * (*second_stage_rate)[entry]);
        }
        ++rhs_call;
        return qmomkit::CqmomMomentRate{required_rates_from_full(full_rates_from_full(stage_two))};
    };
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

std::vector<double> cqmom_native_source_weight_rates()
{
    return {0.02, -0.01, 0.03, -0.015};
}

double cqmom_native_source_weight_rate_moment(const qmomkit::MultiIndexView index)
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
    // advection half-step. For dx/dt = a x and dy/dt = c y, raw moment M_ij
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
    // The manufactured source is dM_alpha/dt = gamma M_alpha. The phase-space evolution
    // framework advances that source with one SSP-RK3 step, whose stability
    // polynomial is P(z) = 1 + z + z^2/2 + z^3/6. The split reference is
    // therefore advection(dt/2), source(dt), advection(dt/2).
    std::vector<double> exact(initial.values.size(), 0.0);
    const double source_factor = ssp_rk3_linear_affine_factor(cqmom_combined_source_gamma * time_step);
    for (std::size_t flat = 0; flat < initial.values.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, initial.extents);
        const double advection_exponent = static_cast<double>(index[0]) * cqmom_a + static_cast<double>(index[1]) * cqmom_c;
        exact[flat] = initial.values[flat] * std::exp(advection_exponent * time_step) * source_factor;
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
        record_audit_failure(audit, "CQMOM source context size did not match the expected phase_space_evolved moment count.");
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
        const qmomkit::MultiIndexView index = context.index(entry);
        if (index != expected_indices[entry])
        {
            record_audit_failure(audit, "CQMOM source context index order did not match the expected public moment order.");
            continue;
        }
        if (audit.invocation_count == 1 && entry < expected_initial_values.size() &&
            std::abs(context.value(entry) - expected_initial_values[entry]) > 1.0e-12)
        {
            record_audit_failure(audit, "CQMOM source context initial value did not match the corresponding phase_space_evolved moment.");
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
    // first advection half-step. Verifying those context values protects the
    // public timing of native CQMOM physics callbacks.
    ++audit.invocation_count;
    const std::vector<qmomkit::MultiIndex> expected_indices = cqmom_native_source_context_indices();
    if (context.size() != expected_indices.size())
    {
        record_audit_failure(audit, "combined CQMOM source context size did not match the expected phase_space_evolved moment count.");
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
        const qmomkit::MultiIndexView index = context.index(entry);
        if (index != expected_indices[entry])
        {
            record_audit_failure(audit, "combined CQMOM source context index order did not match the expected public moment order.");
            continue;
        }
        if (audit.invocation_count == 1 && entry < expected_source_values.size() &&
            std::abs(context.value(entry) - expected_source_values[entry]) > 5.0e-11)
        {
            record_audit_failure(audit, "combined CQMOM source context value did not match the post-advection-half-step moment.");
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
                                         const qmomkit::MultiIndexView index)
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
        record_audit_failure(audit, "CQMOM boundary context size did not match the expected phase_space_evolved moment count.");
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
        const qmomkit::MultiIndexView index = context.index(entry);
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
    // Verify CQMOM phase-space evolution on a two-dimensional commuting analytical
    // problem. CQMOM source callbacks act on the phase-space evolution moment layout, so the
    // final phase_space_evolved moments and quadrature reconstruction are compared with
    // the exact tensor exponential reference on that layout.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load commuting CQMOM configuration.");
        return verification;
    }

    const qmomkit::CqmomCase& configured_cqmom_case = configured.value();
    const qmomkit::CqmomMomentSource source = cqmom_source_for_context;
    const auto options = cqmom_options(configured_cqmom_case);
    const double dt = 0.04;
    const std::size_t steps = 5;
    const double elapsed = dt * static_cast<double>(steps);
    qmomkit::CqmomMomentData moments = configured_cqmom_case.moments;
    double time = 0.0;
    bool diagnostics_successful = true;
    std::size_t retry_count = 0;
    qmomkit::TestCqmomPhaseSpaceEvolutionResult result{
        configured_cqmom_case.moments, std::nullopt, qmomkit::PhaseSpaceEvolutionStatus::success, {}, 0.0, ""};
    for (std::size_t step = 0; step < steps; ++step)
    {
        result = qmomkit::advance_cqmom_phase_space(moments, time, dt, cqmom_velocity, source, options);
        if (!result.ok())
        {
            append_failure(verification, "CQMOM phase-space evolution integrator failed: " + result.message);
            return verification;
        }
        if (!cqmom_diagnostics_are_successful(result.diagnostics))
        {
            append_failure(verification, "CQMOM phase-space evolution diagnostics did not report successful completion of every substage.");
            diagnostics_successful = false;
        }
        retry_count += result.diagnostics.moment_source_retries;
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
    const std::vector<double> exact = cqmom_exact_values(configured_cqmom_case.moments, elapsed);
    update_errors(summary.errors, result.moments.values, exact);
    summary.error_norm = scaled_infinity_error(result.moments.values, exact);
    if (result.quadrature.has_value())
    {
        const qmomkit::CqmomMomentData reconstructed = qmomkit::reconstruct_cqmom_moments(*result.quadrature, configured_cqmom_case.moments.layout);
        update_reconstruction_errors(summary.errors, reconstructed.values, exact);
    }
    else
    {
        summary.success = false;
        append_failure(verification, "CQMOM phase-space evolution result did not retain final quadrature.");
    }
    if (std::abs(time - elapsed) > 1.0e-12)
    {
        summary.success = false;
        append_failure(verification, "CQMOM final time does not match the requested final time.");
    }
    if (summary.errors.moment_scaled > 3.0e-6 || summary.errors.reconstruction_scaled > 3.0e-6)
    {
        summary.success = false;
        append_failure(verification, "CQMOM commuting phase-space evolution exceeded the analytical tolerance.");
    }
    verification.runs.push_back(summary);
    return verification;
}

VerificationResult run_noncommuting_cqmom(const std::filesystem::path& configuration_path)
{
    // Public full CQMOM phase-space evolution verification for the noncommuting operators
    // C: dx/dt = x, dy/dt = 0 and D: dM_(i,j)/dt = j M_(i+1,j-1).
    // The manufactured discrete quadrature therefore has
    // x(t) = exp(t) x(0) and y(t) = y(0) + (exp(t) - 1) x(0).
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom_prepared(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load noncommuting CQMOM configuration.");
        return verification;
    }

    const qmomkit::CqmomPreparedCase& cqmom_prepared_case = configured.value();
    const qmomkit::CqmomResult initial_result = qmomkit::compute_cqmom(cqmom_prepared_case);
    if (!initial_result.ok())
    {
        append_failure(verification, "noncommuting CQMOM initial reconstruction failed: " + initial_result.message);
        return verification;
    }

    const std::vector<std::size_t> extents = cqmom_prepared_case.hierarchy().moment_extents();
    const double final_time = 1.2;
    const std::vector<double> refinement_steps{0.3, 0.15, 0.075, 0.0375, 0.01875};
    const std::vector<qmomkit::PhaseSpaceEvolutionScheme> schemes{
        qmomkit::PhaseSpaceEvolutionScheme::strang, qmomkit::PhaseSpaceEvolutionScheme::advsrc, qmomkit::PhaseSpaceEvolutionScheme::srcadv};
    std::optional<std::vector<double>> cd_coarse;
    std::optional<std::vector<double>> dc_coarse;

    for (const qmomkit::PhaseSpaceEvolutionScheme scheme : schemes)
    {
        qmomkit::CqmomPhaseSpaceEvolutionOptions options = cqmom_state_options();
        options.scheme = scheme;
        const std::size_t first_run_index = verification.runs.size();
        double previous_error = 0.0;
        for (const double dt : refinement_steps)
        {
            const std::size_t steps = static_cast<std::size_t>(std::lround(final_time / dt));
            qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
            double time = 0.0;
            std::optional<qmomkit::CqmomPhaseSpaceEvolutionResult> result;
            bool run_success = true;
            bool diagnostics_successful = true;
            for (std::size_t step = 0; step < steps; ++step)
            {
                result =
                    qmomkit::advance_cqmom_phase_space(state, time, dt, cqmom_noncommuting_velocity, make_cqmom_noncommuting_source(dt), options);
                if (!result->ok())
                {
                    append_failure(verification, "noncommuting CQMOM phase-space evolution failed: " + result->message);
                    run_success = false;
                    break;
                }
                diagnostics_successful = diagnostics_successful && cqmom_split_diagnostics_are_successful(result->diagnostics, scheme);
                state = result->state;
                time = result->final_time;
            }

            RunSummary summary;
            summary.scheme = scheme;
            summary.dt = dt;
            summary.steps = steps;
            summary.final_time = time;
            if (!run_success || !result.has_value())
            {
                summary.success = false;
                summary.error_norm = std::numeric_limits<double>::infinity();
                verification.runs.push_back(summary);
                continue;
            }

            const std::vector<double> exact = cqmom_noncommuting_exact_values(initial_result.quadrature, extents, final_time);
            const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(cqmom_prepared_case.hierarchy());
            const qmomkit::CqmomMomentData exact_required = qmomkit::extract_cqmom_moments(qmomkit::RawMomentTensor{extents, exact}, layout);
            update_errors(summary.errors, state.moment_values(), exact_required.values);
            bool quadrature_successful = false;
            if (result->ok())
            {
                quadrature_successful = true;
                const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(result->state.quadrature(), extents);
                update_errors(summary.errors, reconstructed.values, exact);
                update_reconstruction_errors(summary.errors, reconstructed.values, exact);
                summary.error_norm =
                    std::max(scaled_infinity_error(state.moment_values(), exact_required.values), scaled_infinity_error(reconstructed.values, exact));
                if (scheme == qmomkit::PhaseSpaceEvolutionScheme::advsrc && dt == refinement_steps.front())
                {
                    cd_coarse = reconstructed.values;
                }
                if (scheme == qmomkit::PhaseSpaceEvolutionScheme::srcadv && dt == refinement_steps.front())
                {
                    dc_coarse = reconstructed.values;
                }
            }
            else
            {
                summary.success = false;
                append_failure(verification, "noncommuting CQMOM phase-space evolution did not retain its final quadrature.");
            }
            summary.retry_count = result->diagnostics.moment_source_retries;
            summary.success = run_success && diagnostics_successful && result->ok() && quadrature_successful;
            if (previous_error > 0.0 && summary.error_norm > 0.0)
            {
                summary.observed_order = std::log(previous_error / summary.error_norm) / std::log(2.0);
            }
            previous_error = summary.error_norm;
            verification.runs.push_back(summary);
        }

        if (verification.runs.size() >= first_run_index + refinement_steps.size())
        {
            const double penultimate_order = verification.runs[first_run_index + refinement_steps.size() - 2].observed_order;
            const double finest_order = verification.runs[first_run_index + refinement_steps.size() - 1].observed_order;
            const double expected_order = scheme == qmomkit::PhaseSpaceEvolutionScheme::strang ? 2.0 : 1.0;
            const double order_tolerance = scheme == qmomkit::PhaseSpaceEvolutionScheme::strang ? 0.2 : 0.2;
            if (penultimate_order < expected_order - order_tolerance || finest_order < expected_order - order_tolerance)
            {
                append_failure(verification, "noncommuting CQMOM " + std::string(phase_space_evolution_scheme_name(scheme)) +
                                                 " refinement did not show the expected splitting order.");
            }
        }
    }

    if (!cd_coarse.has_value() || !dc_coarse.has_value() || scaled_infinity_error(*cd_coarse, *dc_coarse) < 1.0e-5)
    {
        append_failure(verification, "noncommuting CQMOM advsrc and srcadv schemes did not differ at finite timestep.");
    }
    return verification;
}

VerificationResult run_cqmom_boundary_contribution(const std::filesystem::path& configuration_path)
{
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom_prepared(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load CQMOM boundary contribution configuration.");
        return verification;
    }

    const qmomkit::CqmomPreparedCase& cqmom_prepared_case = configured.value();
    if (cqmom_prepared_case.hierarchy().node_counts() != std::vector<std::size_t>{2, 2} ||
        cqmom_prepared_case.hierarchy().conditioning_order() != std::vector<std::size_t>{1, 0})
    {
        append_failure(verification, "CQMOM boundary verification requires a fixed 2x2 hierarchy with conditioning_order = {1, 0}.");
        return verification;
    }

    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_prepared_case);
    if (!initial.ok())
    {
        append_failure(verification, "initial CQMOM boundary contribution inversion failed: " + initial.message);
        return verification;
    }

    const std::vector<std::size_t> extents{4, 4};
    const qmomkit::RawMomentTensor initial_moments = qmomkit::reconstruct_raw_moments(initial.quadrature, extents);
    const std::vector<double> weight_rates = cqmom_boundary_weight_rates(initial.quadrature);
    CqmomNativeSourceAudit audit;
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto boundary = [&initial, &weight_rates, &audit](const double time, const qmomkit::CqmomSourceContext& context)
    { return cqmom_boundary_contribution_callback(initial.quadrature, weight_rates, audit, time, context); };

    qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
    qmomkit::CqmomPhaseSpaceEvolutionOptions options = cqmom_state_options();
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    constexpr double dt = 0.02;
    const qmomkit::CqmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_cqmom_phase_space(state, 0.0, dt, zero_velocity, cqmom_state_zero_source, boundary, options);

    if (!result.ok())
    {
        append_failure(verification, "CQMOM boundary contribution phase-space evolution failed: " + result.message);
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
    summary.retry_count = result.diagnostics.moment_source_retries;
    summary.success = result.ok() && audit.passed && audit.invocation_count > 0 && cqmom_diagnostics_are_successful(result.diagnostics);
    if (result.ok())
    {
        const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(result.state.quadrature(), extents);
        update_errors(summary.errors, reconstructed.values, exact);
        update_reconstruction_errors(summary.errors, reconstructed.values, exact);
        summary.error_norm = scaled_infinity_error(reconstructed.values, exact);
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
    const auto configured = qmomkit::config::load_cqmom_prepared(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load combined CQMOM phase-space evolution and moment-source configuration.");
        return verification;
    }

    const qmomkit::CqmomPreparedCase& cqmom_prepared_case = configured.value();
    if (cqmom_prepared_case.hierarchy().node_counts() != std::vector<std::size_t>{2, 2} ||
        cqmom_prepared_case.hierarchy().conditioning_order() != std::vector<std::size_t>{0, 1})
    {
        append_failure(verification, "combined CQMOM verification requires a fixed 2x2 hierarchy with conditioning_order = {0, 1}.");
        return verification;
    }

    // This case exercises the public workflow used by coupled CQMOM physics
    // models: a configured state undergoes phase-space advection, is exposed to
    // a CQMOM-native moment-source callback, then is reconstructed after the final
    // advection half-step.
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_prepared_case);
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

    qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
    qmomkit::CqmomPhaseSpaceEvolutionOptions options = cqmom_state_options();
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;

    const qmomkit::CqmomPhaseSpaceEvolutionResult result = qmomkit::advance_cqmom_phase_space(state, 0.0, dt, cqmom_velocity, source, options);
    if (!result.ok())
    {
        append_failure(verification, "combined CQMOM phase-space evolution and moment-source workflow failed: " + result.message);
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
        append_failure(verification, "combined CQMOM source callback did not inspect the post-advection-half-step context values.");
    }
    if (!cqmom_diagnostics_are_successful(result.diagnostics))
    {
        append_failure(verification, "combined CQMOM phase-space evolution diagnostics did not report successful completion of every substage.");
    }
    if (result.diagnostics.total_node_advection_steps == 0)
    {
        append_failure(verification, "combined CQMOM phase-space advection did not report any node-advection steps.");
    }
    if (result.diagnostics.moment_source_diagnostics.ssp_rk3_status != qmomkit::SspRk3Status::success ||
        std::abs(result.diagnostics.moment_source_diagnostics.time_step_used - dt) > 1.0e-14)
    {
        append_failure(verification, "combined CQMOM moment-source diagnostics did not report the expected successful SSP-RK3 stage.");
    }

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = dt;
    summary.steps = 1;
    summary.final_time = result.final_time;
    summary.retry_count = result.diagnostics.moment_source_retries;
    summary.success = result.ok() && audit.passed && audit.invocation_count > 0 && audit.checked_split_context_values &&
                      cqmom_diagnostics_are_successful(result.diagnostics) && result.diagnostics.total_node_advection_steps > 0 &&
                      result.diagnostics.moment_source_diagnostics.ssp_rk3_status == qmomkit::SspRk3Status::success &&
                      std::abs(result.diagnostics.moment_source_diagnostics.time_step_used - dt) <= 1.0e-14;
    const std::vector<double> exact = cqmom_combined_native_source_exact_values(initial_moments, dt);
    if (result.ok())
    {
        const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(result.state.quadrature(), extents);
        update_errors(summary.errors, reconstructed.values, exact);
        update_reconstruction_errors(summary.errors, reconstructed.values, exact);
        summary.error_norm = scaled_infinity_error(reconstructed.values, exact);
        summary.errors.conservation = std::abs(reconstructed.values.front() - exact.front());
    }
    else
    {
        summary.success = false;
        append_failure(verification, "combined CQMOM phase-space evolution result did not retain final quadrature.");
    }
    if (std::abs(result.final_time - dt) > 1.0e-12)
    {
        summary.success = false;
        append_failure(verification, "combined CQMOM final time does not match the requested final time.");
    }
    if (summary.errors.moment_scaled > 5.0e-9 || summary.errors.reconstruction_scaled > 5.0e-9)
    {
        summary.success = false;
        append_failure(verification, "combined CQMOM moments exceeded the analytical split-reference tolerance.");
    }
    verification.runs.push_back(summary);
    return verification;
}

VerificationResult run_cqmom_native_source(const std::filesystem::path& configuration_path)
{
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom_prepared(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load CQMOM-native source configuration.");
        return verification;
    }

    const qmomkit::CqmomPreparedCase& cqmom_prepared_case = configured.value();
    if (cqmom_prepared_case.hierarchy().node_counts() != std::vector<std::size_t>{2, 2} ||
        cqmom_prepared_case.hierarchy().conditioning_order() != std::vector<std::size_t>{0, 1})
    {
        append_failure(verification, "CQMOM-native source verification requires a fixed 2x2 hierarchy with conditioning_order = {0, 1}.");
        return verification;
    }

    // This is the normal user workflow: load_cqmom_prepared(...) owns the configured
    // state, compute_cqmom(...) verifies the initial analytical quadrature, and
    // CqmomPreparedCase::initial_state() is passed to advance_cqmom_phase_space(...).
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_prepared_case);
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
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });

    qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
    qmomkit::CqmomPhaseSpaceEvolutionOptions options = cqmom_state_options();
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;

    constexpr double dt = 0.08;
    const qmomkit::CqmomPhaseSpaceEvolutionResult result = qmomkit::advance_cqmom_phase_space(state, 0.0, dt, zero_velocity, source, options);
    if (!result.ok())
    {
        append_failure(verification, "CQMOM-native moment-source evolution failed: " + result.message);
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
        append_failure(verification, "CQMOM-native moment-source evolution diagnostics did not report successful completion of every substage.");
    }
    const bool diagnostics_successful = cqmom_diagnostics_are_successful(result.diagnostics);

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = dt;
    summary.steps = 1;
    summary.final_time = result.final_time;
    summary.retry_count = result.diagnostics.moment_source_retries;
    summary.success = result.ok() && audit.passed && audit.invocation_count > 0 && audit.checked_initial_context_values && diagnostics_successful;
    const std::vector<double> exact = cqmom_native_source_exact_values(initial_moments, dt);
    if (result.ok())
    {
        const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(result.state.quadrature(), extents);
        update_errors(summary.errors, reconstructed.values, exact);
        update_reconstruction_errors(summary.errors, reconstructed.values, exact);
        summary.error_norm = scaled_infinity_error(reconstructed.values, exact);
        summary.errors.conservation = std::abs(reconstructed.values.front() - exact.front());
    }
    else
    {
        summary.success = false;
        append_failure(verification, "CQMOM-native moment-source evolution result did not retain final quadrature.");
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

void verify_cqmom_advection_failure_result(const std::string& case_name, const qmomkit::CqmomPhaseSpaceEvolutionResult& result,
                                           const std::string& expected_message_fragment, VerificationResult& verification)
{
    // This is the CQMOM state counterpart of the QMOM velocity-callback
    // failure: the public result should identify phase-space advection as the
    // failing operation and should not report source or reconstruction success
    // for stages that were never reached.
    if (result.ok())
    {
        append_failure(verification, case_name + " unexpectedly reported CQMOM phase-space evolution success.");
    }
    if (result.status != qmomkit::PhaseSpaceEvolutionStatus::advection_failed)
    {
        append_failure(verification, case_name + " did not report PhaseSpaceEvolutionStatus::advection_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::PhaseSpaceEvolutionStage::first_advection_half_step)
    {
        append_failure(verification, case_name + " did not identify the first advection half-step as the failed stage.");
    }
    if (result.diagnostics.first_advection_status != qmomkit::NodeAdvectionStatus::velocity_failed)
    {
        append_failure(verification, case_name + " did not report NodeAdvectionStatus::velocity_failed.");
    }
    if (!contains_text(result.message, expected_message_fragment) || !contains_text(result.diagnostics.message, expected_message_fragment) ||
        !contains_text(result.diagnostics.first_advection_diagnostics.message, expected_message_fragment))
    {
        append_failure(verification, case_name + " did not preserve a useful velocity-callback failure message.");
    }
    if (std::abs(result.final_time) > 1.0e-14)
    {
        append_failure(verification, case_name + " advanced final_time despite phase-space advection failure.");
    }
    if (result.diagnostics.moment_source_diagnostics.accepted || result.diagnostics.moment_source_diagnostics.rejected ||
        result.diagnostics.first_reconstruction_diagnostics.success || result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, case_name + " reported later phase-space evolution stage diagnostics after the failed advection substep.");
    }
}

void verify_cqmom_rk23_minimum_step_acceptance_result(const qmomkit::CqmomPhaseSpaceEvolutionResult& result, const qmomkit::RawMomentTensor& exact,
                                                      VerificationResult& verification)
{
    // This mirrors the QMOM RK23 minimum-step case through the public CqmomState
    // workflow. The user velocity is valid; the adaptive node integrator should
    // advance at the effective minimum step instead of treating it as a callback
    // or integration failure.
    if (!result.ok())
    {
        append_failure(verification, "CQMOM RK23 minimum-step case did not report phase-space evolution success.");
    }
    if (result.status != qmomkit::PhaseSpaceEvolutionStatus::success)
    {
        append_failure(verification, "CQMOM RK23 minimum-step case did not report PhaseSpaceEvolutionStatus::success.");
    }
    if (result.diagnostics.failed_stage != qmomkit::PhaseSpaceEvolutionStage::none)
    {
        append_failure(verification, "CQMOM RK23 minimum-step case reported a failed phase-space evolution stage.");
    }
    if (result.diagnostics.first_advection_status != qmomkit::NodeAdvectionStatus::success ||
        result.diagnostics.second_advection_status != qmomkit::NodeAdvectionStatus::success ||
        result.diagnostics.first_advection_diagnostics.adaptive_status != qmomkit::AdaptiveRkStatus::success ||
        result.diagnostics.second_advection_diagnostics.adaptive_status != qmomkit::AdaptiveRkStatus::success)
    {
        append_failure(verification, "CQMOM RK23 minimum-step case did not preserve successful RK23 diagnostics.");
    }
    if (result.diagnostics.first_advection_diagnostics.rejected_step_count == 0 || result.diagnostics.total_node_advection_rejected_steps == 0 ||
        result.diagnostics.first_advection_diagnostics.step_count == 0 || result.diagnostics.second_advection_diagnostics.step_count == 0)
    {
        append_failure(verification, "CQMOM RK23 minimum-step case did not preserve accepted and rejected RK23 diagnostics.");
    }
    if (std::abs(result.final_time - 0.08) > 1.0e-14 || std::abs(result.diagnostics.first_advection_diagnostics.final_time - 0.04) > 1.0e-14 ||
        std::abs(result.diagnostics.second_advection_diagnostics.final_time - 0.08) > 1.0e-14)
    {
        append_failure(verification, "CQMOM RK23 minimum-step case did not reach the expected endpoint.");
    }
    if (!result.diagnostics.moment_source_diagnostics.accepted || !result.diagnostics.first_reconstruction_diagnostics.success ||
        !result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, "CQMOM RK23 minimum-step case did not complete the later phase-space evolution stages.");
    }
    if (!values_are_finite(result.state.moment_values()))
    {
        append_failure(verification, "CQMOM RK23 minimum-step case returned a non-finite final stored moment.");
    }
    const qmomkit::CqmomMomentLayout phase_space_evolution_layout = qmomkit::create_cqmom_moment_layout(result.state.hierarchy());
    const qmomkit::CqmomMomentData exact_required = qmomkit::extract_cqmom_moments(exact, phase_space_evolution_layout);
    if (result.state.moment_values().size() != exact_required.values.size())
    {
        append_failure(verification, "CQMOM RK23 minimum-step case returned an unexpected stored moment count.");
    }
    if (result.state.moment_values().size() == exact_required.values.size() &&
        scaled_infinity_error(result.state.moment_values(), exact_required.values) > 3.0e-6)
    {
        append_failure(verification, "CQMOM RK23 minimum-step case final stored moments exceeded the analytical advection-reference tolerance.");
    }
    if (!active_branch_counts_are_valid(result.state.active_branch_counts()))
    {
        append_failure(verification, "CQMOM RK23 minimum-step case returned invalid active branch counts.");
    }
    if (!result.ok())
    {
        append_failure(verification, "CQMOM RK23 minimum-step case did not retain final quadrature.");
    }
    else
    {
        if (!cqmom_quadrature_values_are_finite(result.state.quadrature()))
        {
            append_failure(verification, "CQMOM RK23 minimum-step case returned a non-finite final terminal node or weight.");
        }
        if (!cqmom_quadrature_topology_is_consistent(result.state.quadrature()))
        {
            append_failure(verification, "CQMOM RK23 minimum-step case returned inconsistent final quadrature topology.");
        }
        const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(result.state.quadrature(), exact.extents);
        if (scaled_infinity_error(reconstructed.values, exact.values) > 3.0e-6)
        {
            append_failure(verification, "CQMOM RK23 minimum-step case final quadrature exceeded the analytical advection-reference tolerance.");
        }
    }
}

void verify_cqmom_post_source_inversion_failure_result(const qmomkit::CqmomPhaseSpaceEvolutionResult& result, const double expected_failure_time,
                                                       VerificationResult& verification)
{
    // The source callback is finite and correctly sized, so the moment-source update
    // itself succeeds. It deliberately drives M[0,0] negative, which is not a
    // realizable CQMOM moment set, and the following public CQMOM inversion
    // must stop the phase-space evolution workflow before the second advection half-step.
    if (result.ok())
    {
        append_failure(verification, "post-source inversion failure unexpectedly reported CQMOM phase-space evolution success.");
    }
    if (result.status != qmomkit::PhaseSpaceEvolutionStatus::inversion_failed)
    {
        append_failure(verification, "post-source inversion failure did not report PhaseSpaceEvolutionStatus::inversion_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::PhaseSpaceEvolutionStage::post_source_inversion)
    {
        append_failure(verification, "post-source inversion failure did not identify post_source_inversion as the failed stage.");
    }
    if (result.diagnostics.moment_source_status != qmomkit::MomentSourceUpdateStatus::success ||
        result.diagnostics.moment_source_diagnostics.ssp_rk3_status != qmomkit::SspRk3Status::success ||
        !result.diagnostics.moment_source_diagnostics.accepted)
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
        append_failure(verification, "post-source inversion failure did not stop at the time of the moment-source stage.");
    }
    if (result.diagnostics.second_advection_status != qmomkit::NodeAdvectionStatus::success ||
        result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, "post-source inversion failure reported later-stage completion after the failed inversion.");
    }
}

void verify_cqmom_initial_inversion_failure_result(const qmomkit::CqmomPhaseSpaceEvolutionResult& result, VerificationResult& verification)
{
    // This represents a configured CQMOM state whose selected moments are
    // finite but not realizable. Phase-space evolution should fail during the initial
    // public CQMOM reconstruction, before phase-space advection or moment-source
    // integration can produce misleading progress diagnostics.
    if (result.ok())
    {
        append_failure(verification, "CQMOM initial inversion failure unexpectedly reported phase-space evolution success.");
    }
    if (result.status != qmomkit::PhaseSpaceEvolutionStatus::inversion_failed)
    {
        append_failure(verification, "CQMOM initial inversion failure did not report PhaseSpaceEvolutionStatus::inversion_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::PhaseSpaceEvolutionStage::initial_inversion)
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
    if (result.diagnostics.moment_source_diagnostics.accepted || result.diagnostics.first_reconstruction_diagnostics.success ||
        result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, "CQMOM initial inversion failure reported later-stage completion.");
    }
}

void verify_cqmom_case_advection_failure_result(const std::string& case_name, const qmomkit::TestCqmomPhaseSpaceEvolutionResult& result,
                                                const qmomkit::CqmomMomentData& initial_moments, const std::string& expected_message_fragment,
                                                VerificationResult& verification)
{
    // A velocity-callback failure has the same contract here as for a state
    // created by a prepared case: advection_failed at the first advection
    // half-step and no source or reconstruction progress. This helper also
    // checks that the failed step leaves the moments unchanged.
    if (result.ok())
    {
        append_failure(verification, case_name + " unexpectedly reported CQMOM phase-space evolution success.");
    }
    if (result.status != qmomkit::PhaseSpaceEvolutionStatus::advection_failed)
    {
        append_failure(verification, case_name + " did not report PhaseSpaceEvolutionStatus::advection_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::PhaseSpaceEvolutionStage::first_advection_half_step)
    {
        append_failure(verification, case_name + " did not identify the first advection half-step as the failed stage.");
    }
    if (result.diagnostics.first_advection_status != qmomkit::NodeAdvectionStatus::velocity_failed)
    {
        append_failure(verification, case_name + " did not report NodeAdvectionStatus::velocity_failed.");
    }
    if (!contains_text(result.message, expected_message_fragment) || !contains_text(result.diagnostics.message, expected_message_fragment) ||
        !contains_text(result.diagnostics.first_advection_diagnostics.message, expected_message_fragment))
    {
        append_failure(verification, case_name + " did not preserve the velocity-callback failure message.");
    }
    if (std::abs(result.final_time) > 1.0e-14)
    {
        append_failure(verification, case_name + " advanced final_time despite phase-space advection failure.");
    }
    if (result.moments.layout.required_indices != initial_moments.layout.required_indices ||
        !vectors_equal(result.moments.values, initial_moments.values, 0.0))
    {
        append_failure(verification, case_name + " changed its moments during a failed advection substep.");
    }
    if (result.diagnostics.moment_source_diagnostics.accepted || result.diagnostics.moment_source_diagnostics.rejected ||
        result.diagnostics.first_reconstruction_diagnostics.success || result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, case_name + " reported later-stage diagnostics after advection failure.");
    }
}

void verify_cqmom_case_source_failure_result(const std::string& case_name, const qmomkit::TestCqmomPhaseSpaceEvolutionResult& result,
                                             const std::string& expected_message_fragment, VerificationResult& verification)
{
    // A CqmomMomentSource callback that returns a wrong-sized rate vector is a
    // user model error, so phase-space evolution must stop at the moment-source
    // stage.
    if (result.ok())
    {
        append_failure(verification, case_name + " unexpectedly reported moment-source evolution success.");
    }
    if (result.status != qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed)
    {
        append_failure(verification, case_name + " did not report PhaseSpaceEvolutionStatus::moment_source_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::PhaseSpaceEvolutionStage::moment_source_step)
    {
        append_failure(verification, case_name + " did not identify the moment-source stage as the failed stage.");
    }
    if (result.diagnostics.moment_source_status != qmomkit::MomentSourceUpdateStatus::source_failed)
    {
        append_failure(verification, case_name + " did not report MomentSourceUpdateStatus::source_failed.");
    }
    if (!contains_text(result.message, expected_message_fragment) || !contains_text(result.diagnostics.message, expected_message_fragment) ||
        !contains_text(result.diagnostics.moment_source_diagnostics.message, expected_message_fragment))
    {
        append_failure(verification, case_name + " did not preserve the source failure message.");
    }
    if (result.diagnostics.moment_source_diagnostics.accepted || result.diagnostics.moment_source_diagnostics.rejected)
    {
        append_failure(verification, case_name + " reported accepted/rejected source integration after callback failure.");
    }
    if (result.diagnostics.post_source_cqmom_status != qmomkit::CqmomStatus::success || !result.diagnostics.post_source_cqmom_diagnostics.empty() ||
        result.diagnostics.final_reconstruction_diagnostics.success)
    {
        append_failure(verification, case_name + " reported post-source progress after source failure.");
    }
}

void verify_cqmom_case_post_source_inversion_failure_result(const qmomkit::TestCqmomPhaseSpaceEvolutionResult& result,
                                                            const double expected_failure_time, VerificationResult& verification)
{
    // The source is finite and correctly sized, so the moment-source update
    // succeeds, but it drives M[0,0] negative. The following CQMOM inversion
    // must fail with diagnostics and stop before the second advection half-step.
    if (result.ok())
    {
        append_failure(verification, "post-source inversion failure unexpectedly reported phase-space evolution success.");
    }
    if (result.status != qmomkit::PhaseSpaceEvolutionStatus::inversion_failed)
    {
        append_failure(verification, "post-source inversion failure did not report PhaseSpaceEvolutionStatus::inversion_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::PhaseSpaceEvolutionStage::post_source_inversion)
    {
        append_failure(verification, "post-source inversion failure did not identify post_source_inversion as the failed stage.");
    }
    if (result.diagnostics.moment_source_status != qmomkit::MomentSourceUpdateStatus::success ||
        result.diagnostics.moment_source_diagnostics.ssp_rk3_status != qmomkit::SspRk3Status::success ||
        !result.diagnostics.moment_source_diagnostics.accepted)
    {
        append_failure(verification, "post-source inversion failure did not preserve successful moment-source diagnostics.");
    }
    if (result.diagnostics.post_source_cqmom_status == qmomkit::CqmomStatus::success || result.diagnostics.post_source_cqmom_diagnostics.empty())
    {
        append_failure(verification, "post-source inversion failure did not preserve failed CQMOM diagnostics.");
    }
    if (result.message.empty() || result.diagnostics.message.empty() ||
        (!contains_text(result.message, "m0") && !contains_text(result.message, "realizable")))
    {
        append_failure(verification, "post-source inversion failure did not provide a meaningful reconstruction message.");
    }
    if (std::abs(result.final_time - expected_failure_time) > 1.0e-14)
    {
        append_failure(verification, "post-source inversion failure did not stop at the time of the moment-source stage.");
    }
    if (result.diagnostics.final_reconstruction_diagnostics.success || result.diagnostics.second_advection_diagnostics.step_count > 0)
    {
        append_failure(verification, "post-source inversion failure reported later-stage completion.");
    }
}

void verify_helper_based_common_cqmom_phase_space_evolution(const qmomkit::CqmomCase& configured_cqmom_case, VerificationResult& verification)
{
    qmomkit::CqmomCase policy_cqmom_case = configured_cqmom_case;
    policy_cqmom_case.qmom_options.adaptive = true;
    policy_cqmom_case.adaptive_branch_counts = true;
    policy_cqmom_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::basic;
    policy_cqmom_case.vandermonde_solver = qmomkit::VandermondeSolver::dense;
    qmomkit::TestCqmomPhaseSpaceEvolutionOptions options = cqmom_options(policy_cqmom_case);

    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
    const qmomkit::TestCqmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_cqmom_phase_space(policy_cqmom_case.moments, 0.0, 0.02, zero_velocity, zero_source, options);
    if (!result.ok())
    {
        append_failure(verification, "helper-based CQMOM policy probe failed: " + result.message);
    }
    else
    {
        verify_helper_based_common_cqmom_policies("helper-based CQMOM", result.diagnostics,
                                                  policy_cqmom_case.moments.layout.hierarchy.number_of_dimensions(), verification);
    }

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = 0.02;
    summary.steps = 1;
    summary.final_time = result.final_time;
    summary.retry_count = result.diagnostics.moment_source_retries;
    summary.success =
        result.ok() && result.diagnostics.cqmom_diagnostic_level == qmomkit::CqmomDiagnosticLevel::basic &&
        result.diagnostics.cqmom_vandermonde_solver == qmomkit::VandermondeSolver::dense &&
        cqmom_diagnostics_enable_adaptive_branch_counts(result.diagnostics.initial_cqmom_diagnostics) &&
        cqmom_diagnostics_request_solver(result.diagnostics.initial_cqmom_diagnostics,
                                         policy_cqmom_case.moments.layout.hierarchy.number_of_dimensions(), qmomkit::VandermondeSolver::dense);
    verification.runs.push_back(summary);
}

VerificationResult run_cqmom_phase_space_evolution_failures(const std::filesystem::path& configuration_path)
{
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom_prepared(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load CQMOM phase-space evolution failure configuration.");
        return verification;
    }

    const qmomkit::CqmomPreparedCase& cqmom_prepared_case = configured.value();
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_prepared_case);
    if (!initial.ok())
    {
        append_failure(verification, "initial CQMOM phase-space evolution failure inversion failed: " + initial.message);
        return verification;
    }

    qmomkit::CqmomPhaseSpaceEvolutionOptions options = cqmom_state_options();
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    constexpr double dt = 0.08;
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    struct AdvectionFailureCase
    {
        std::string name;
        qmomkit::PhaseSpaceVelocity velocity;
        std::string expected_message_fragment;
    };

    const std::vector<AdvectionFailureCase> advection_failure_cases{
        {"CQMOM wrong velocity dimension",
         [](double, const std::vector<double>&, std::vector<double>& velocities, const std::size_t dimensions)
         {
             // Simulates a CQMOM user callback that fills one node's worth of
             // components instead of one per node.
             velocities.assign(dimensions, 0.0);
         },
         "Phase-space velocity callback returned the wrong dimension:"},
        {"CQMOM non-finite velocity",
         [](double, const std::vector<double>&, std::vector<double>& velocities, std::size_t)
         {
             // Simulates a coupled model producing a non-finite node velocity.
             std::fill(velocities.begin(), velocities.end(), 0.0);
             velocities[1] = std::numeric_limits<double>::infinity();
         },
         "Phase-space velocity callback returned a non-finite value:"}};

    for (const AdvectionFailureCase& failure_case : advection_failure_cases)
    {
        qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
        const qmomkit::CqmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_cqmom_phase_space(state, 0.0, dt, failure_case.velocity, zero_source, options);
        verify_cqmom_advection_failure_result(failure_case.name, result, failure_case.expected_message_fragment, verification);

        RunSummary summary;
        summary.scheme = options.scheme;
        summary.dt = dt;
        summary.steps = 1;
        summary.final_time = result.final_time;
        summary.retry_count = result.diagnostics.moment_source_retries;
        summary.success = !result.ok() && result.status == qmomkit::PhaseSpaceEvolutionStatus::advection_failed &&
                          result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::first_advection_half_step &&
                          result.diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::velocity_failed;
        verification.runs.push_back(summary);
    }

    qmomkit::CqmomPhaseSpaceEvolutionOptions rk_minimum_options = options;
    rk_minimum_options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    rk_minimum_options.node_options.adaptive_options.initial_step = 0.04;
    rk_minimum_options.node_options.adaptive_options.max_step = 0.04;
    rk_minimum_options.node_options.adaptive_options.min_step = 0.02;
    const qmomkit::PhaseSpaceVelocity stiff_velocity = qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>& xi)
        {
            // Valid two-dimensional velocity with tight RK23 tolerances. RK23
            // should accept finite updates at the effective minimum step instead of
            // reporting callback or integration failure.
            return std::vector<double>{xi[0], -xi[1]};
        });
    qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
    const qmomkit::CqmomPhaseSpaceEvolutionResult rk_minimum =
        qmomkit::advance_cqmom_phase_space(state, 0.0, dt, stiff_velocity, zero_source, rk_minimum_options);
    const std::vector<std::size_t> rk_minimum_extents{4, 4};
    const qmomkit::RawMomentTensor rk_minimum_initial = qmomkit::reconstruct_raw_moments(initial.quadrature, rk_minimum_extents);
    qmomkit::RawMomentTensor rk_minimum_exact{rk_minimum_extents, std::vector<double>(qmomkit::product_extents(rk_minimum_extents), 0.0)};
    for (std::size_t flat = 0; flat < rk_minimum_exact.values.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, rk_minimum_extents);
        const double exponent = static_cast<double>(index[0]) - static_cast<double>(index[1]);
        rk_minimum_exact.values[flat] = rk_minimum_initial.values[flat] * std::exp(exponent * dt);
    }
    verify_cqmom_rk23_minimum_step_acceptance_result(rk_minimum, rk_minimum_exact, verification);

    RunSummary rk_summary;
    rk_summary.scheme = rk_minimum_options.scheme;
    rk_summary.dt = dt;
    rk_summary.steps = 1;
    rk_summary.final_time = rk_minimum.final_time;
    rk_summary.retry_count = rk_minimum.diagnostics.moment_source_retries;
    rk_summary.success = rk_minimum.ok() && rk_minimum.status == qmomkit::PhaseSpaceEvolutionStatus::success &&
                         rk_minimum.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::none &&
                         rk_minimum.diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::success &&
                         rk_minimum.diagnostics.second_advection_status == qmomkit::NodeAdvectionStatus::success &&
                         rk_minimum.diagnostics.total_node_advection_rejected_steps > 0 && rk_minimum.diagnostics.total_node_advection_steps > 0;
    verification.runs.push_back(rk_summary);

    state = cqmom_prepared_case.initial_state();
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
    const qmomkit::CqmomPhaseSpaceEvolutionResult post_source_inversion =
        qmomkit::advance_cqmom_phase_space(state, 0.0, dt, zero_velocity, nonrealizable_source, options);
    verify_cqmom_post_source_inversion_failure_result(post_source_inversion, 0.0, verification);

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = dt;
    summary.steps = 1;
    summary.final_time = post_source_inversion.final_time;
    summary.retry_count = post_source_inversion.diagnostics.moment_source_retries;
    summary.success = !post_source_inversion.ok() && post_source_inversion.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed &&
                      post_source_inversion.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::post_source_inversion &&
                      post_source_inversion.diagnostics.post_source_cqmom_status != qmomkit::CqmomStatus::success;
    verification.runs.push_back(summary);
    return verification;
}

VerificationResult run_cqmom_phase_space_evolution_inversion_failures(const std::filesystem::path& configuration_path)
{
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom_prepared(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load CQMOM phase-space evolution inversion-failure configuration.");
        return verification;
    }

    const auto configured_cqmom_case_result = qmomkit::config::load_cqmom(configuration_path);
    if (!configured_cqmom_case_result.ok())
    {
        print_configuration_errors(configured_cqmom_case_result);
        append_failure(verification, "unable to load CQMOM phase-space evolution inversion-failure configuration as a CqmomCase.");
        return verification;
    }

    qmomkit::CqmomState state = qmomkit::internal::CqmomConstructionAccess::create_state(configured_cqmom_case_result.value());
    qmomkit::CqmomPhaseSpaceEvolutionOptions options = cqmom_state_options();
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const qmomkit::CqmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 0.08, zero_velocity, cqmom_state_zero_source, options);
    verify_cqmom_initial_inversion_failure_result(result, verification);

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = 0.08;
    summary.steps = 1;
    summary.final_time = result.final_time;
    summary.retry_count = result.diagnostics.moment_source_retries;
    summary.success = !result.ok() && result.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed &&
                      result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::initial_inversion &&
                      result.diagnostics.initial_cqmom_status != qmomkit::CqmomStatus::success;
    verification.runs.push_back(summary);
    return verification;
}

VerificationResult run_cqmom_case_phase_space_evolution_failures(const std::filesystem::path& configuration_path)
{
    // Failure propagation for states built directly from the configured
    // CqmomCase, with no accepted quadrature, so phase-space evolution performs
    // the initial inversion: failing velocity and source callbacks, a source
    // that makes the post-source inversion fail, and the case policies reaching
    // the CQMOM reconstructions. run_cqmom_phase_space_evolution_failures covers
    // states created by a prepared case.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load CQMOM phase-space evolution failure configuration.");
        return verification;
    }

    const qmomkit::CqmomCase& configured_cqmom_case = configured.value();
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(configured_cqmom_case);
    if (!initial.ok())
    {
        append_failure(verification, "initial CQMOM phase-space evolution failure inversion failed: " + initial.message);
        return verification;
    }
    verify_helper_based_common_cqmom_phase_space_evolution(configured_cqmom_case, verification);

    qmomkit::TestCqmomPhaseSpaceEvolutionOptions options = cqmom_options(configured_cqmom_case);
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    constexpr double dt = 0.08;
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    struct AdvectionFailureCase
    {
        std::string name;
        qmomkit::PhaseSpaceVelocity velocity;
        std::string expected_message_fragment;
    };
    const std::vector<AdvectionFailureCase> advection_failure_cases{
        {"CQMOM wrong velocity dimension", [](double, const std::vector<double>&, std::vector<double>& velocities, const std::size_t dimensions)
         { velocities.assign(dimensions, 0.0); }, "Phase-space velocity callback returned the wrong dimension:"},
        {"CQMOM non-finite velocity",
         [](double, const std::vector<double>&, std::vector<double>& velocities, std::size_t)
         {
             std::fill(velocities.begin(), velocities.end(), 0.0);
             velocities[1] = std::numeric_limits<double>::infinity();
         },
         "Phase-space velocity callback returned a non-finite value:"}};

    for (const AdvectionFailureCase& failure_case : advection_failure_cases)
    {
        const qmomkit::TestCqmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_cqmom_phase_space(configured_cqmom_case.moments, 0.0, dt, failure_case.velocity, zero_source, options);
        verify_cqmom_case_advection_failure_result(failure_case.name, result, configured_cqmom_case.moments, failure_case.expected_message_fragment,
                                                   verification);

        RunSummary summary;
        summary.scheme = options.scheme;
        summary.dt = dt;
        summary.steps = 1;
        summary.final_time = result.final_time;
        summary.retry_count = result.diagnostics.moment_source_retries;
        summary.success = !result.ok() && result.status == qmomkit::PhaseSpaceEvolutionStatus::advection_failed &&
                          result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::first_advection_half_step &&
                          result.diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::velocity_failed;
        verification.runs.push_back(summary);
    }

    const qmomkit::CqmomMomentSource wrong_size_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size() + 1, 0.0)}; };
    const qmomkit::TestCqmomPhaseSpaceEvolutionResult source_failure =
        qmomkit::advance_cqmom_phase_space(configured_cqmom_case.moments, 0.0, dt, zero_velocity, wrong_size_source, options);
    verify_cqmom_case_source_failure_result("CQMOM wrong moment-source rate size", source_failure,
                                            "CQMOM moment-source rate values must match the source context size.", verification);

    RunSummary source_summary;
    source_summary.scheme = options.scheme;
    source_summary.dt = dt;
    source_summary.steps = 1;
    source_summary.final_time = source_failure.final_time;
    source_summary.retry_count = source_failure.diagnostics.moment_source_retries;
    source_summary.success = !source_failure.ok() && source_failure.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed &&
                             source_failure.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step;
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
    const qmomkit::TestCqmomPhaseSpaceEvolutionResult post_source_inversion =
        qmomkit::advance_cqmom_phase_space(configured_cqmom_case.moments, 0.0, dt, zero_velocity, nonrealizable_source, options);
    verify_cqmom_case_post_source_inversion_failure_result(post_source_inversion, 0.0, verification);

    RunSummary inversion_summary;
    inversion_summary.scheme = options.scheme;
    inversion_summary.dt = dt;
    inversion_summary.steps = 1;
    inversion_summary.final_time = post_source_inversion.final_time;
    inversion_summary.retry_count = post_source_inversion.diagnostics.moment_source_retries;
    inversion_summary.success = !post_source_inversion.ok() && post_source_inversion.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed &&
                                post_source_inversion.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::post_source_inversion &&
                                post_source_inversion.diagnostics.post_source_cqmom_status != qmomkit::CqmomStatus::success;
    verification.runs.push_back(inversion_summary);
    return verification;
}

void verify_cqmom_native_source_failure_result(const std::string& case_name, const qmomkit::CqmomPhaseSpaceEvolutionResult& result,
                                               const std::string& expected_message_fragment, VerificationResult& verification)
{
    // User callback failures should stop during the moment-source stage,
    // before post-source inversion or later phase-space evolution stages can hide the
    // callback contract violation.
    if (result.ok())
    {
        append_failure(verification, case_name + " unexpectedly reported phase-space evolution success.");
    }
    if (result.status != qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed)
    {
        append_failure(verification, case_name + " did not report PhaseSpaceEvolutionStatus::moment_source_failed.");
    }
    if (result.diagnostics.failed_stage != qmomkit::PhaseSpaceEvolutionStage::moment_source_step)
    {
        append_failure(verification, case_name + " did not identify the moment-source stage as the failed stage.");
    }
    if (result.diagnostics.moment_source_status != qmomkit::MomentSourceUpdateStatus::source_failed)
    {
        append_failure(verification, case_name + " did not report MomentSourceUpdateStatus::source_failed.");
    }
    if (result.message.find(expected_message_fragment) == std::string::npos)
    {
        append_failure(verification, case_name + " result message did not explain the source callback problem.");
    }
    if (result.diagnostics.message.find(expected_message_fragment) == std::string::npos)
    {
        append_failure(verification, case_name + " phase-space evolution diagnostics did not carry the source callback message.");
    }
    if (result.diagnostics.moment_source_diagnostics.message.find(expected_message_fragment) == std::string::npos)
    {
        append_failure(verification, case_name + " moment-source diagnostics did not carry the source callback message.");
    }
    if (result.diagnostics.moment_source_diagnostics.accepted || result.diagnostics.moment_source_diagnostics.rejected)
    {
        append_failure(verification, case_name + " moment-source diagnostics reported accepted/rejected RK state after callback failure.");
    }
    if (result.diagnostics.post_source_cqmom_status != qmomkit::CqmomStatus::success || !result.diagnostics.post_source_cqmom_diagnostics.empty())
    {
        append_failure(verification, case_name + " unexpectedly reached post-source CQMOM reconstruction.");
    }
    if (result.diagnostics.final_reconstruction_diagnostics.success || result.diagnostics.second_advection_diagnostics.step_count > 0 ||
        std::abs(result.final_time) > 1.0e-14)
    {
        append_failure(verification, case_name + " advanced final_time past the source callback failure.");
    }
}

VerificationResult run_cqmom_native_source_failures(const std::filesystem::path& configuration_path)
{
    // Native CQMOM source failures model user physics callbacks that return a
    // rate vector with the wrong size or non-finite entries. The phase-space evolution
    // result should fail at the source stage while preserving useful callback
    // diagnostics.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom_prepared(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load CQMOM-native source failure configuration.");
        return verification;
    }

    const qmomkit::CqmomPreparedCase& cqmom_prepared_case = configured.value();
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_prepared_case);
    if (!initial.ok())
    {
        append_failure(verification, "initial CQMOM-native source failure inversion failed: " + initial.message);
        return verification;
    }

    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    qmomkit::CqmomPhaseSpaceEvolutionOptions options = cqmom_state_options();
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
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
                                                  "CQMOM moment-source rate values must match the source context size."},
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
                                                  "CQMOM moment-source rate values must be finite."},
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
                                                  "CQMOM moment-source rate values must be finite."}};

    for (const FailureCase& failure_case : failure_cases)
    {
        qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
        const qmomkit::CqmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_cqmom_phase_space(state, 0.0, dt, zero_velocity, failure_case.source, options);
        verify_cqmom_native_source_failure_result(failure_case.name, result, failure_case.expected_message_fragment, verification);

        RunSummary summary;
        summary.scheme = options.scheme;
        summary.dt = dt;
        summary.steps = 1;
        summary.final_time = result.final_time;
        summary.retry_count = result.diagnostics.moment_source_retries;
        summary.success = !result.ok() && result.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed &&
                          result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step &&
                          result.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::source_failed;
        verification.runs.push_back(summary);
    }
    return verification;
}

VerificationResult run_commuting_cqmom_state_adaptive_2d(const std::filesystem::path& configuration_path)
{
    // Adaptive CQMOM phase-space evolution through CqmomState uses a deterministic affine
    // branch-reduction case with commuting advection/source operators. It
    // verifies that ragged adaptive quadrature remains valid for phase-space evolution and
    // reconstructs the exact analytical moment tensor after splitting.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom_prepared(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load adaptive CQMOM state configuration.");
        return verification;
    }

    const qmomkit::CqmomPreparedCase& cqmom_prepared_case = configured.value();
    if (!cqmom_prepared_case.adaptive_branch_counts() || !cqmom_prepared_case.qmom_options().adaptive)
    {
        append_failure(verification, "adaptive CQMOM state verification requires adaptive_branch_counts and adaptive QMOM.");
        return verification;
    }

    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_prepared_case);
    if (!initial.ok())
    {
        append_failure(verification, "initial adaptive CQMOM state inversion failed: " + initial.message);
        return verification;
    }
    if (!initial.quadrature.adaptive_branch_counts_enabled() || initial.quadrature.has_rectangular_storage())
    {
        append_failure(verification, "initial adaptive CQMOM state quadrature did not use ragged adaptive storage.");
    }
    if (initial.quadrature.total_nodes() >= cqmom_prepared_case.hierarchy().total_nodes())
    {
        append_failure(verification, "initial adaptive CQMOM state quadrature did not reduce active node count below maximum capacity.");
    }
    if (initial.quadrature.terminal_paths().size() != initial.quadrature.weights().size())
    {
        append_failure(verification, "initial adaptive CQMOM state quadrature did not provide one terminal path per node.");
    }

    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(cqmom_prepared_case.hierarchy());
    const qmomkit::CqmomMomentData initial_moments = qmomkit::reconstruct_cqmom_moments(initial.quadrature, layout);
    qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
    const qmomkit::CqmomPhaseSpaceEvolutionOptions options = cqmom_state_options();
    const double dt = 0.04;
    const std::size_t steps = 5;
    const double elapsed = dt * static_cast<double>(steps);
    double time = 0.0;
    bool diagnostics_successful = true;
    std::size_t retry_count = 0;
    std::optional<qmomkit::CqmomPhaseSpaceEvolutionResult> result;
    for (std::size_t step = 0; step < steps; ++step)
    {
        result = qmomkit::advance_cqmom_phase_space(state, time, dt, cqmom_velocity, cqmom_state_zero_source, options);
        if (!result->ok())
        {
            append_failure(verification, "adaptive CQMOM phase-space evolution through CqmomState failed: " + result->message);
            return verification;
        }
        if (!cqmom_diagnostics_are_successful(result->diagnostics))
        {
            append_failure(
                verification,
                "adaptive CQMOM phase-space evolution through CqmomState diagnostics did not report successful completion of every substage.");
            diagnostics_successful = false;
        }
        retry_count += result->diagnostics.moment_source_retries;
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
    const std::vector<double> exact = cqmom_advection_exact_values(initial_moments, elapsed);
    if (result->ok())
    {
        if (!result->state.quadrature().adaptive_branch_counts_enabled() || result->state.quadrature().has_rectangular_storage())
        {
            summary.success = false;
            append_failure(verification, "final adaptive CQMOM state quadrature did not preserve ragged adaptive storage.");
        }
        if (result->state.quadrature().total_nodes() > cqmom_prepared_case.hierarchy().total_nodes())
        {
            summary.success = false;
            append_failure(verification, "final adaptive CQMOM state quadrature exceeded configured maximum node capacity.");
        }
        if (result->state.quadrature().terminal_paths().size() != result->state.quadrature().weights().size())
        {
            summary.success = false;
            append_failure(verification, "final adaptive CQMOM state quadrature did not provide one terminal path per node.");
        }
        const qmomkit::CqmomMomentData reconstructed = qmomkit::reconstruct_cqmom_moments(result->state.quadrature(), layout);
        update_errors(summary.errors, reconstructed.values, exact);
        update_reconstruction_errors(summary.errors, reconstructed.values, exact);
        summary.error_norm = scaled_infinity_error(reconstructed.values, exact);
    }
    else
    {
        summary.success = false;
        append_failure(verification, "adaptive CQMOM phase-space evolution through CqmomState result did not retain final quadrature.");
    }
    if (std::abs(time - elapsed) > 1.0e-12)
    {
        summary.success = false;
        append_failure(verification, "adaptive CQMOM state final time does not match the requested final time.");
    }
    if (summary.errors.moment_scaled > 3.0e-6 || summary.errors.reconstruction_scaled > 3.0e-6)
    {
        summary.success = false;
        append_failure(verification, "adaptive CQMOM state commuting phase-space evolution exceeded the analytical tolerance.");
    }
    verification.runs.push_back(summary);
    return verification;
}

VerificationResult run_automatic_scaling_cqmom_state_adaptive_2d(const std::filesystem::path& configuration_path, const std::vector<double>& units)
{
    // This manufactured CQMOM state is the ragged rule
    // (-1,-1,0.10), (1,-1,0.30), (0.5,2,0.60), conditioned as (xi1,xi0).
    // The two fixtures express that same physical rule in units (1,1) and
    // (1e6,1e-6). Constant physical advection is exact under SSP-RK3, while
    // dM_alpha/dt = gamma M_alpha has the independent SSP-RK3 polynomial
    // oracle P(gamma dt). The source callback visits public terminal nodes,
    // proving source-stage reconstruction presents physical coordinates.
    VerificationResult verification;
    const auto configured = qmomkit::config::load_cqmom_prepared(configuration_path);
    if (!configured.ok())
    {
        print_configuration_errors(configured);
        append_failure(verification, "unable to load automatic-scaling CQMOM evolution configuration.");
        return verification;
    }

    const qmomkit::CqmomPreparedCase& cqmom_prepared_case = configured.value();
    if (!cqmom_prepared_case.qmom_options().automatic_scaling || !cqmom_prepared_case.qmom_options().adaptive ||
        !cqmom_prepared_case.adaptive_branch_counts() || cqmom_prepared_case.hierarchy().conditioning_order() != std::vector<std::size_t>({1, 0}))
    {
        append_failure(verification, "automatic-scaling evolution configuration did not retain its adaptive non-identity CQMOM policy.");
        return verification;
    }

    constexpr double dt = 0.2;
    constexpr double gamma = 0.1;
    const std::vector<std::vector<double>> initial_canonical{{-1.0, -1.0}, {1.0, -1.0}, {0.5, 2.0}};
    const std::vector<double> initial_weights{0.10, 0.30, 0.60};
    const std::vector<double> canonical_velocity{0.4, -0.25};
    const double source_factor = ssp_rk3_linear_affine_factor(gamma * dt);
    bool callback_coordinates_physical = true;
    std::size_t callback_visits = 0;

    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity([&units, &canonical_velocity](double, const std::vector<double>&)
                                           { return std::vector<double>{canonical_velocity[0] * units[0], canonical_velocity[1] * units[1]}; });
    const auto source = [&units, &initial_canonical, &canonical_velocity, &callback_coordinates_physical,
                         &callback_visits](double, const qmomkit::CqmomSourceContext& context)
    {
        context.for_each_terminal_node(
            [&](double, const std::vector<double>& coordinates)
            {
                ++callback_visits;
                bool matched = false;
                for (const std::vector<double>& initial : initial_canonical)
                {
                    const double x = (initial[0] + 0.5 * dt * canonical_velocity[0]) * units[0];
                    const double y = (initial[1] + 0.5 * dt * canonical_velocity[1]) * units[1];
                    if (std::abs(coordinates[0] - x) <= 2.0e-9 * std::max(1.0, std::abs(x)) &&
                        std::abs(coordinates[1] - y) <= 2.0e-9 * std::max(1.0, std::abs(y)))
                    {
                        matched = true;
                    }
                }
                callback_coordinates_physical = callback_coordinates_physical && matched;
            });
        std::vector<double> rates(context.size(), 0.0);
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            rates[entry] = gamma * context.value(entry);
        }
        return qmomkit::CqmomMomentRate{std::move(rates)};
    };

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_cqmom_phase_space(cqmom_prepared_case.initial_state(), 0.0, dt, velocity, source, options);

    RunSummary summary;
    summary.scheme = options.scheme;
    summary.dt = dt;
    summary.steps = 1;
    summary.final_time = result.final_time;
    summary.retry_count = result.diagnostics.moment_source_retries;
    summary.success = result.ok() && cqmom_diagnostics_are_successful(result.diagnostics) && callback_coordinates_physical && callback_visits != 0;
    if (!summary.success)
    {
        append_failure(verification, "automatic-scaling CQMOM evolution did not complete all physical reconstruction stages.");
    }
    if (!result.ok())
    {
        append_failure(verification, "automatic-scaling CQMOM evolution did not retain a final quadrature.");
        summary.success = false;
        verification.runs.push_back(summary);
        return verification;
    }

    const qmomkit::CqmomQuadrature& quadrature = result.state.quadrature();
    if (quadrature.has_rectangular_storage() || !quadrature.adaptive_branch_counts_enabled() || quadrature.total_nodes() != initial_weights.size())
    {
        append_failure(verification, "automatic-scaling CQMOM evolution did not preserve the expected ragged topology.");
        summary.success = false;
    }
    std::vector<std::size_t> actual_paths;
    for (const qmomkit::ConditioningPathView path : quadrature.terminal_paths())
    {
        actual_paths.push_back(path.node_index_at_level(0) * 2 + path.node_index_at_level(1));
    }
    std::sort(actual_paths.begin(), actual_paths.end());
    if (actual_paths != std::vector<std::size_t>({0, 1, 2}))
    {
        append_failure(verification, "automatic-scaling CQMOM evolution changed terminal conditioning paths.");
        summary.success = false;
    }

    std::vector<std::vector<double>> expected_coordinates;
    std::vector<double> expected_weights;
    for (std::size_t node = 0; node < initial_weights.size(); ++node)
    {
        expected_coordinates.push_back({(initial_canonical[node][0] + dt * canonical_velocity[0]) * units[0],
                                        (initial_canonical[node][1] + dt * canonical_velocity[1]) * units[1]});
        expected_weights.push_back(initial_weights[node] * source_factor);
    }
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        const qmomkit::ConditioningPathView path = quadrature.terminal_paths()[node];
        const std::size_t expected_node = path.node_index_at_level(0) == 0 ? path.node_index_at_level(1) : 2;
        update_errors(summary.errors, {quadrature.weights()[node]}, {expected_weights[expected_node]});
        for (std::size_t dimension = 0; dimension < units.size(); ++dimension)
        {
            update_errors(summary.errors, {quadrature.abscissas()[dimension][node]}, {expected_coordinates[expected_node][dimension]});
        }
    }
    std::vector<double> analytical(16, 0.0);
    for (std::size_t flat = 0; flat < analytical.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, {4, 4});
        for (std::size_t node = 0; node < expected_weights.size(); ++node)
        {
            analytical[flat] += expected_weights[node] * integer_power(expected_coordinates[node][0], index[0]) *
                                integer_power(expected_coordinates[node][1], index[1]);
        }
    }
    const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(quadrature, {4, 4});
    update_reconstruction_errors(summary.errors, reconstructed.values, analytical);
    summary.error_norm = scaled_infinity_error(reconstructed.values, analytical);
    if (summary.errors.moment_scaled > 2.0e-9 || summary.errors.reconstruction_scaled > 2.0e-9)
    {
        append_failure(verification, "automatic-scaling CQMOM evolution exceeded the manufactured physical tolerance.");
        summary.success = false;
    }
    verification.runs.push_back(summary);
    return verification;
}

// Silent on success; CTest surfaces stderr with --output-on-failure.
void report_failure(const std::string& case_name, const VerificationResult& result)
{
    std::cerr << std::setprecision(17);
    std::cerr << "Case: " << case_name << '\n';
    std::cerr << "scheme dt steps final_time error_norm observed_order max_abs max_scaled reconstruction_abs reconstruction_scaled m0_error "
                 "retry_count success\n";
    for (const RunSummary& run : result.runs)
    {
        std::cerr << phase_space_evolution_scheme_name(run.scheme) << ' ' << run.dt << ' ' << run.steps << ' ' << run.final_time << ' '
                  << run.error_norm << ' ' << run.observed_order << ' ' << run.errors.moment_absolute << ' ' << run.errors.moment_scaled << ' '
                  << run.errors.reconstruction_absolute << ' ' << run.errors.reconstruction_scaled << ' ' << run.errors.conservation << ' '
                  << run.retry_count << ' ' << run.success << '\n';
    }
}

int run_cqmom_automatic_scaling_is_unit_invariant()
{
    const std::filesystem::path case_directory = std::filesystem::path(__FILE__).parent_path();
    const VerificationResult reference =
        run_automatic_scaling_cqmom_state_adaptive_2d(case_directory / "automatic_scaling_cqmom_state_adaptive_reference.qmom", {1.0, 1.0});
    const VerificationResult rescaled =
        run_automatic_scaling_cqmom_state_adaptive_2d(case_directory / "automatic_scaling_cqmom_state_adaptive_rescaled.qmom", {1.0e6, 1.0e-6});

    if (!reference.passed || !rescaled.passed)
    {
        report_failure("automatic_scaling_cqmom_state_adaptive_reference", reference);
        report_failure("automatic_scaling_cqmom_state_adaptive_rescaled", rescaled);
    }
    if (!reference.passed)
    {
        std::cerr << reference.message << '\n';
    }
    if (!rescaled.passed)
    {
        std::cerr << rescaled.message << '\n';
    }
    return reference.passed && rescaled.passed ? 0 : 1;
}

void print_usage()
{
    std::cerr << "Usage: phase_space_evolution_analytical_verification <configuration.qmom> <case_identifier> | "
                 "cqmom_automatic_scaling_is_unit_invariant\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc == 2 && std::string(argv[1]) == "cqmom_automatic_scaling_is_unit_invariant")
    {
        return run_cqmom_automatic_scaling_is_unit_invariant();
    }
    if (argc != 3)
    {
        print_usage();
        return 2;
    }

    const std::filesystem::path configuration_path = argv[1];
    const std::string case_name = argv[2];

    VerificationResult result;
    if (case_name == "commuting_qmom")
    {
        result = run_commuting_qmom(configuration_path);
    }
    else if (case_name == "qmom_boundary_contribution")
    {
        result = run_qmom_boundary_contribution(configuration_path);
    }
    else if (case_name == "qmom_phase_space_evolution_failures")
    {
        result = run_qmom_phase_space_evolution_failures(configuration_path);
    }
    else if (case_name == "qmom_phase_space_evolution_inversion_failures")
    {
        result = run_qmom_phase_space_evolution_inversion_failures(configuration_path);
    }
    else if (case_name == "noncommuting_qmom")
    {
        result = run_noncommuting_qmom(configuration_path);
    }
    else if (case_name == "noncommuting_cqmom")
    {
        result = run_noncommuting_cqmom(configuration_path);
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
    else if (case_name == "cqmom_phase_space_evolution_failures")
    {
        result = run_cqmom_phase_space_evolution_failures(configuration_path);
    }
    else if (case_name == "cqmom_phase_space_evolution_inversion_failures")
    {
        result = run_cqmom_phase_space_evolution_inversion_failures(configuration_path);
    }
    else if (case_name == "cqmom_case_phase_space_evolution_failures")
    {
        result = run_cqmom_case_phase_space_evolution_failures(configuration_path);
    }
    else if (case_name == "commuting_cqmom_state_adaptive_2d")
    {
        result = run_commuting_cqmom_state_adaptive_2d(configuration_path);
    }
    else
    {
        std::cerr << "Unknown phase-space evolution analytical case '" << case_name << "'.\n";
        return 2;
    }

    if (!result.passed)
    {
        report_failure(case_name, result);
    }
    if (!result.passed)
    {
        std::cerr << result.message << '\n';
        return 1;
    }
    return 0;
}
