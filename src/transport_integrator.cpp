#include "qmomkit/transport_integrator.hpp"

#include "internal/cqmom_case.hpp"
#include "internal/cqmom_validation.hpp"
#include "internal/moment_reconstruction.hpp"
#include "internal/numeric_utils.hpp"
#include "internal/performance.hpp"
#include "qmomkit/quadrature.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace qmomkit
{

// Hybrid transport splits each macro step into two numerical operators:
// C: continuous phase-space transport, which advects quadrature nodes and
//    reconstructs raw moments from the advected quadrature.
// D: discontinuous moment evolution, which combines MomentSource with any
//    MomentBoundaryContribution, advances moments with SSP-RK3, then
//    reconstructs quadrature to check realizability before more node motion.
//
// Supported operator orderings:
//   strang: C(dt/2) -> D(dt) -> C(dt/2)
//   cd:     C(dt)   -> D(dt)
//   dc:     D(dt)   -> C(dt)

namespace
{

enum class CqmomMetadataTrust
{
    public_input,
    trusted_case_state
};

void validate_time_input(const double t, const double time_step)
{
    if (!std::isfinite(t) || !std::isfinite(time_step) || !std::isfinite(t + time_step))
    {
        throw std::invalid_argument("Transport time input contains non-finite values.");
    }
    if (time_step <= 0.0)
    {
        throw std::invalid_argument("Transport timestep must be positive.");
    }
}

void validate_transport_scheme(const TransportScheme scheme)
{
    switch (scheme)
    {
    case TransportScheme::strang:
    case TransportScheme::cd:
    case TransportScheme::dc:
        return;
    }
    throw std::invalid_argument("Unsupported transport scheme.");
}

template <typename ApplyContinuous, typename ApplySource>
[[nodiscard]] bool execute_transport_sequence(const TransportScheme scheme, const double t, const double time_step,
                                              ApplyContinuous&& apply_continuous, ApplySource&& apply_source)
{
    const double half_step = 0.5 * time_step;
    switch (scheme)
    {
    // strang: C(dt/2) -> D(dt) -> C(dt/2)
    case TransportScheme::strang:
        if (!apply_continuous(t, half_step, false))
        {
            return false;
        }
        if (!apply_source(t + half_step))
        {
            return false;
        }
        return apply_continuous(t + half_step, half_step, true);
    // cd: C(dt) -> D(dt)
    case TransportScheme::cd:
        if (!apply_continuous(t, time_step, false))
        {
            return false;
        }
        return apply_source(t);
    // dc: D(dt) -> C(dt)
    case TransportScheme::dc:
        if (!apply_source(t))
        {
            return false;
        }
        return apply_continuous(t, time_step, true);
    }
    throw std::invalid_argument("Unsupported transport scheme.");
}

void validate_qmom_input(const std::vector<double>& moments, const std::size_t node_count, const double t, const double time_step,
                         const PhaseSpaceVelocity& velocity, const MomentSource& source, const TransportScheme scheme)
{
    validate_time_input(t, time_step);
    validate_transport_scheme(scheme);
    if (node_count == 0)
    {
        throw std::invalid_argument("QMOM transport requires at least one node.");
    }
    if (moments.size() < 2 * node_count)
    {
        throw std::invalid_argument("QMOM transport input has too few moments for inversion.");
    }
    if (!internal::all_finite(moments))
    {
        throw std::invalid_argument("QMOM transport moments must be finite.");
    }
    if (!velocity)
    {
        throw std::invalid_argument("Transport phase-space velocity callback must be provided.");
    }
    if (!source)
    {
        throw std::invalid_argument("Transport moment-source callback must be provided.");
    }
}

void validate_boundary(const MomentBoundaryContribution* const boundary)
{
    if (boundary != nullptr && !*boundary)
    {
        throw std::invalid_argument("Moment-boundary contribution callback must be provided.");
    }
}

void validate_cqmom_input(const RawMomentTensor& moments, const PhaseSpaceHierarchy& hierarchy, const double t, const double time_step,
                          const PhaseSpaceVelocity& velocity, const MomentSource& source, const TransportScheme scheme)
{
    validate_time_input(t, time_step);
    validate_transport_scheme(scheme);
    internal::validate_dense_cqmom_moment_tensor(moments, hierarchy, "CQMOM transport moment tensor rank must match the phase-space hierarchy.",
                                                 "CQMOM transport moment tensor size does not match its extents.",
                                                 "CQMOM transport moments must be finite.");
    if (!velocity)
    {
        throw std::invalid_argument("Transport phase-space velocity callback must be provided.");
    }
    if (!source)
    {
        throw std::invalid_argument("Transport moment-source callback must be provided.");
    }
}

void validate_cqmom_boundary(const CqmomMomentBoundaryContribution* const boundary)
{
    if (boundary != nullptr && !*boundary)
    {
        throw std::invalid_argument("CQMOM boundary contribution callback must be provided.");
    }
}

void validate_lean_cqmom_input(const CqmomLeanTransportState& state, const double t, const double time_step, const PhaseSpaceVelocity& velocity,
                               const TransportScheme scheme)
{
    validate_time_input(t, time_step);
    validate_transport_scheme(scheme);
    internal::validate_fixed_cqmom_moment_data(state.moments, {"Lean CQMOM transport requires the fixed CQMOM moment layout.",
                                                               "Lean CQMOM transport moments must match the fixed layout size.",
                                                               "Lean CQMOM transport moments must be finite."});
    if (!velocity)
    {
        throw std::invalid_argument("Transport phase-space velocity callback must be provided.");
    }
}

void validate_trusted_lean_cqmom_input(const CqmomLeanTransportState& state, const double t, const double time_step,
                                       const PhaseSpaceVelocity& velocity, const TransportScheme scheme)
{
    validate_time_input(t, time_step);
    validate_transport_scheme(scheme);
    // CqmomState owns already validated layout metadata. The trusted path still
    // checks the evolving values that transport mutates between reconstructions.
    internal::validate_cqmom_moment_values_match_layout(
        state.moments, {"", "Lean CQMOM transport moments must match the fixed layout size.", "Lean CQMOM transport moments must be finite."});
    if (!velocity)
    {
        throw std::invalid_argument("Transport phase-space velocity callback must be provided.");
    }
}

void validate_lean_cqmom_input(const CqmomLeanTransportState& state, const double t, const double time_step, const PhaseSpaceVelocity& velocity,
                               const CqmomMomentSource& source, const TransportScheme scheme)
{
    validate_lean_cqmom_input(state, t, time_step, velocity, scheme);
    if (!source)
    {
        throw std::invalid_argument("Lean CQMOM transport moment-source callback must be provided.");
    }
}

void validate_trusted_lean_cqmom_input(const CqmomLeanTransportState& state, const double t, const double time_step,
                                       const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source, const TransportScheme scheme)
{
    validate_trusted_lean_cqmom_input(state, t, time_step, velocity, scheme);
    if (!source)
    {
        throw std::invalid_argument("Lean CQMOM transport moment-source callback must be provided.");
    }
}

[[nodiscard]] MomentReconstructionDiagnostics successful_reconstruction(const std::size_t moment_count)
{
    return MomentReconstructionDiagnostics{true, moment_count, ""};
}

void record_continuous_diagnostics(TransportDiagnostics& diagnostics, const bool final_continuous, const NodeAdvectionStatus status,
                                   const NodeAdvectionDiagnostics& continuous_diagnostics)
{
    if (final_continuous)
    {
        diagnostics.second_continuous_status = status;
        diagnostics.second_continuous_diagnostics = continuous_diagnostics;
    }
    else
    {
        diagnostics.first_continuous_status = status;
        diagnostics.first_continuous_diagnostics = continuous_diagnostics;
    }
}

template <typename SourceUpdateResult>
void record_source_update_diagnostics(TransportDiagnostics& diagnostics, const SourceUpdateResult& source_update)
{
    diagnostics.source_status = source_update.status;
    diagnostics.source_diagnostics = source_update.diagnostics;
    diagnostics.source_retries = source_update.diagnostics.retries;
}

[[nodiscard]] QmomResult qmom_from_quadrature(const Quadrature1D& quadrature)
{
    QmomStatus status = QmomStatus::success;
    switch (quadrature.status)
    {
    case QuadratureStatus::success:
        status = QmomStatus::success;
        break;
    case QuadratureStatus::invalid_input:
        status = QmomStatus::invalid_input;
        break;
    case QuadratureStatus::not_realizable:
        status = QmomStatus::not_realizable;
        break;
    case QuadratureStatus::no_adaptive_quadrature:
        status = QmomStatus::no_adaptive_quadrature;
        break;
    }
    return QmomResult{quadrature.abscissas, quadrature.weights, status, quadrature.message};
}

void add_continuous_totals(TransportDiagnostics& diagnostics)
{
    diagnostics.total_node_advection_steps =
        diagnostics.first_continuous_diagnostics.step_count + diagnostics.second_continuous_diagnostics.step_count;
    diagnostics.total_node_advection_rejected_steps =
        diagnostics.first_continuous_diagnostics.rejected_step_count + diagnostics.second_continuous_diagnostics.rejected_step_count;
    diagnostics.total_node_advection_rhs_evaluations =
        diagnostics.first_continuous_diagnostics.rhs_evaluations + diagnostics.second_continuous_diagnostics.rhs_evaluations;
}

void finalize_diagnostics(TransportDiagnostics& diagnostics)
{
    add_continuous_totals(diagnostics);
    diagnostics.source_retries = diagnostics.source_diagnostics.retries;
}

void record_failure_diagnostics(TransportDiagnostics& diagnostics, const TransportStage stage, const std::string& message)
{
    diagnostics.failed_stage = stage;
    diagnostics.message = message;
    finalize_diagnostics(diagnostics);
}

[[nodiscard]] QmomTransportResult qmom_failure(std::vector<double> moments, QmomResult quadrature, const TransportStatus status,
                                               const TransportStage stage, TransportDiagnostics diagnostics, const double final_time,
                                               std::string message)
{
    record_failure_diagnostics(diagnostics, stage, message);
    return QmomTransportResult{std::move(moments), std::move(quadrature), status, std::move(diagnostics), final_time, std::move(message)};
}

[[nodiscard]] CqmomDenseTransportResult cqmom_failure(RawMomentTensor moments, std::optional<CqmomQuadrature> quadrature,
                                                      const TransportStatus status, const TransportStage stage, TransportDiagnostics diagnostics,
                                                      const double final_time, std::string message)
{
    record_failure_diagnostics(diagnostics, stage, message);
    return CqmomDenseTransportResult{std::move(moments), std::move(quadrature), status, std::move(diagnostics), final_time, std::move(message)};
}

[[nodiscard]] CqmomLeanTransportResult lean_cqmom_failure(CqmomLeanTransportState state, std::optional<CqmomQuadrature> quadrature,
                                                          const TransportStatus status, const TransportStage stage, TransportDiagnostics diagnostics,
                                                          const double final_time, std::string message)
{
    record_failure_diagnostics(diagnostics, stage, message);
    return CqmomLeanTransportResult{std::move(state), std::move(quadrature), status, std::move(diagnostics), final_time, std::move(message)};
}

[[nodiscard]] std::string neutral_cqmom_message(std::string message)
{
    constexpr std::string_view lean_prefix = "Lean CQMOM ";
    constexpr std::string_view neutral_prefix = "CQMOM ";
    std::size_t position = message.find(lean_prefix);
    while (position != std::string::npos)
    {
        message.replace(position, lean_prefix.size(), neutral_prefix);
        position = message.find(lean_prefix, position + neutral_prefix.size());
    }
    return message;
}

[[nodiscard]] TransportDiagnostics neutral_cqmom_diagnostics(TransportDiagnostics diagnostics)
{
    diagnostics.message = neutral_cqmom_message(std::move(diagnostics.message));
    diagnostics.source_diagnostics.message = neutral_cqmom_message(std::move(diagnostics.source_diagnostics.message));
    diagnostics.first_reconstruction_diagnostics.message = neutral_cqmom_message(std::move(diagnostics.first_reconstruction_diagnostics.message));
    diagnostics.final_reconstruction_diagnostics.message = neutral_cqmom_message(std::move(diagnostics.final_reconstruction_diagnostics.message));
    return diagnostics;
}

} // namespace

// Public QMOM overloads dispatch here; this helper applies the selected
// C/D numerical transport sequence to the one-dimensional moment vector.
static QmomTransportResult advance_transport_qmom_impl(const std::vector<double>& moments, const std::size_t node_count, const double t,
                                                       const double time_step, const PhaseSpaceVelocity& velocity, const MomentSource& source,
                                                       const MomentBoundaryContribution* const boundary, const QmomTransportOptions& options)
{
    internal::performance::record_transport_macro_step();
    TransportDiagnostics diagnostics;
    diagnostics.scheme = options.scheme;
    try
    {
        validate_qmom_input(moments, node_count, t, time_step, velocity, source, options.scheme);
        validate_boundary(boundary);
    }
    catch (const std::exception& error)
    {
        return qmom_failure(moments, {}, TransportStatus::invalid_input, TransportStage::none, diagnostics, t, error.what());
    }

    QmomResult initial_quadrature = compute_qmom(MomentSet(moments), node_count, options.qmom_options);
    diagnostics.initial_qmom_status = initial_quadrature.status;
    if (!initial_quadrature.ok())
    {
        return qmom_failure(moments, initial_quadrature, TransportStatus::inversion_failed, TransportStage::initial_inversion, diagnostics, t,
                            initial_quadrature.message);
    }

    std::vector<double> current_moments = moments;
    QmomResult current_quadrature = initial_quadrature;
    std::optional<QmomTransportResult> failure;

    const auto apply_continuous = [&](const double start_time, const double step, const bool final_continuous)
    {
        internal::performance::record_continuous_operator_call();
        const QmomNodeAdvectionResult continuous = advance_nodes(current_quadrature, start_time, step, velocity, options.node_options);
        record_continuous_diagnostics(diagnostics, final_continuous, continuous.status, continuous.diagnostics);
        if (!continuous.ok())
        {
            failure = qmom_failure(current_moments, current_quadrature, TransportStatus::continuous_failed,
                                   final_continuous ? TransportStage::second_continuous_half_step : TransportStage::first_continuous_half_step,
                                   diagnostics, start_time, continuous.message);
            return false;
        }

        try
        {
            internal::performance::record_moment_reconstruction_call();
            current_moments = reconstruct_raw_moments(continuous.quadrature, moments.size());
            if (final_continuous)
            {
                diagnostics.final_reconstruction_diagnostics = successful_reconstruction(current_moments.size());
            }
            else
            {
                diagnostics.first_reconstruction_diagnostics = successful_reconstruction(current_moments.size());
            }
            current_quadrature = qmom_from_quadrature(continuous.quadrature);
        }
        catch (const std::exception& error)
        {
            if (final_continuous)
            {
                diagnostics.final_reconstruction_diagnostics = {false, moments.size(), error.what()};
            }
            else
            {
                diagnostics.first_reconstruction_diagnostics = {false, moments.size(), error.what()};
            }
            failure = qmom_failure(current_moments, qmom_from_quadrature(continuous.quadrature), TransportStatus::moment_reconstruction_failed,
                                   final_continuous ? TransportStage::final_moment_reconstruction : TransportStage::first_moment_reconstruction,
                                   diagnostics, start_time + step, error.what());
            return false;
        }
        return true;
    };

    const auto apply_source = [&](const double start_time)
    {
        internal::performance::record_discontinuous_operator_call();
        // D-step: MomentSource and optional boundary contribution share the
        // existing SSP-RK3 source update path.
        const MomentSourceUpdateResult source_update =
            boundary == nullptr ? advance_moments_with_source(current_moments, start_time, time_step, source, options.source_options)
                                : advance_moments_with_source(current_moments, start_time, time_step, source, *boundary, options.source_options);
        record_source_update_diagnostics(diagnostics, source_update);
        if (!source_update.ok())
        {
            failure = qmom_failure(current_moments, current_quadrature, TransportStatus::source_failed, TransportStage::discontinuous_source_step,
                                   diagnostics, start_time, source_update.message);
            return false;
        }

        internal::performance::record_post_source_inversion();
        current_moments = source_update.moments;
        // Reconstruct quadrature after the moment update so realizability is
        // checked before any later continuous stage can use these moments.
        current_quadrature = compute_qmom(MomentSet(current_moments), node_count, options.qmom_options);
        diagnostics.post_source_qmom_status = current_quadrature.status;
        if (!current_quadrature.ok())
        {
            failure = qmom_failure(current_moments, current_quadrature, TransportStatus::inversion_failed, TransportStage::post_source_inversion,
                                   diagnostics, start_time, current_quadrature.message);
            return false;
        }
        return true;
    };

    if (!execute_transport_sequence(options.scheme, t, time_step, apply_continuous, apply_source))
    {
        return std::move(*failure);
    }

    finalize_diagnostics(diagnostics);
    return QmomTransportResult{
        std::move(current_moments), std::move(current_quadrature), TransportStatus::success, std::move(diagnostics), t + time_step, ""};
}

QmomTransportResult advance_transport_qmom(const std::vector<double>& moments, const std::size_t node_count, const double t, const double time_step,
                                           const PhaseSpaceVelocity& velocity, const MomentSource& source, const QmomTransportOptions& options)
{
    return advance_transport_qmom_impl(moments, node_count, t, time_step, velocity, source, nullptr, options);
}

QmomTransportResult advance_transport_qmom(const std::vector<double>& moments, const std::size_t node_count, const double t, const double time_step,
                                           const PhaseSpaceVelocity& velocity, const MomentSource& source, const MomentBoundaryContribution& boundary,
                                           const QmomTransportOptions& options)
{
    return advance_transport_qmom_impl(moments, node_count, t, time_step, velocity, source, &boundary, options);
}

// Public dense CQMOM overloads dispatch here; this helper applies the selected
// C/D numerical transport sequence while the caller owns full raw-tensor storage.
static CqmomDenseTransportResult advance_transport_cqmom_dense_impl(const RawMomentTensor& moments, const PhaseSpaceHierarchy& hierarchy,
                                                                    const double t, const double time_step, const PhaseSpaceVelocity& velocity,
                                                                    const MomentSource& source, const MomentBoundaryContribution* const boundary,
                                                                    const CqmomDenseTransportOptions& options)
{
    internal::performance::record_transport_macro_step();
    TransportDiagnostics diagnostics;
    diagnostics.scheme = options.scheme;
    try
    {
        validate_cqmom_input(moments, hierarchy, t, time_step, velocity, source, options.scheme);
        validate_boundary(boundary);
    }
    catch (const std::exception& error)
    {
        return cqmom_failure(moments, std::nullopt, TransportStatus::invalid_input, TransportStage::none, diagnostics, t, error.what());
    }

    CqmomResult initial_quadrature =
        compute_cqmom({hierarchy, moments, options.conditional_options, options.qmom_options, options.adaptive_branch_counts});
    diagnostics.initial_cqmom_status = initial_quadrature.status;
    diagnostics.initial_cqmom_diagnostics = initial_quadrature.diagnostics;
    if (!initial_quadrature.ok())
    {
        return cqmom_failure(moments, initial_quadrature.quadrature, TransportStatus::inversion_failed, TransportStage::initial_inversion,
                             diagnostics, t, initial_quadrature.message);
    }

    RawMomentTensor current_moments = moments;
    CqmomQuadrature current_quadrature = initial_quadrature.quadrature;
    std::optional<CqmomDenseTransportResult> failure;

    const auto apply_continuous = [&](const double start_time, const double step, const bool final_continuous)
    {
        internal::performance::record_continuous_operator_call();
        const CqmomNodeAdvectionResult continuous = advance_nodes(current_quadrature, start_time, step, velocity, options.node_options);
        record_continuous_diagnostics(diagnostics, final_continuous, continuous.status, continuous.diagnostics);
        if (!continuous.ok())
        {
            failure = cqmom_failure(std::move(current_moments), current_quadrature, TransportStatus::continuous_failed,
                                    final_continuous ? TransportStage::second_continuous_half_step : TransportStage::first_continuous_half_step,
                                    diagnostics, start_time, continuous.message);
            return false;
        }

        try
        {
            internal::performance::record_moment_reconstruction_call();
            current_moments = reconstruct_raw_moments(continuous.quadrature, moments.extents);
            if (final_continuous)
            {
                diagnostics.final_reconstruction_diagnostics = successful_reconstruction(current_moments.values.size());
            }
            else
            {
                diagnostics.first_reconstruction_diagnostics = successful_reconstruction(current_moments.values.size());
            }
            current_quadrature = continuous.quadrature;
        }
        catch (const std::exception& error)
        {
            if (final_continuous)
            {
                diagnostics.final_reconstruction_diagnostics = {false, moments.values.size(), error.what()};
            }
            else
            {
                diagnostics.first_reconstruction_diagnostics = {false, moments.values.size(), error.what()};
            }
            failure = cqmom_failure(std::move(current_moments), continuous.quadrature, TransportStatus::moment_reconstruction_failed,
                                    final_continuous ? TransportStage::final_moment_reconstruction : TransportStage::first_moment_reconstruction,
                                    diagnostics, start_time + step, error.what());
            return false;
        }
        return true;
    };

    const auto apply_source = [&](const double start_time)
    {
        internal::performance::record_discontinuous_operator_call();
        // D-step: MomentSource and optional boundary contribution share the
        // existing SSP-RK3 source update path.
        const MomentSourceUpdateResult source_update =
            boundary == nullptr
                ? advance_moments_with_source(current_moments.values, start_time, time_step, source, options.source_options)
                : advance_moments_with_source(current_moments.values, start_time, time_step, source, *boundary, options.source_options);
        record_source_update_diagnostics(diagnostics, source_update);
        if (!source_update.ok())
        {
            failure = cqmom_failure(std::move(current_moments), current_quadrature, TransportStatus::source_failed,
                                    TransportStage::discontinuous_source_step, diagnostics, start_time, source_update.message);
            return false;
        }

        current_moments = RawMomentTensor{moments.extents, source_update.moments};
        // Reconstruct quadrature after the moment update so CQMOM
        // realizability is checked before any later continuous stage.
        internal::performance::record_post_source_inversion();
        CqmomResult post_source_quadrature =
            compute_cqmom({hierarchy, current_moments, options.conditional_options, options.qmom_options, options.adaptive_branch_counts});
        diagnostics.post_source_cqmom_status = post_source_quadrature.status;
        diagnostics.post_source_cqmom_diagnostics = post_source_quadrature.diagnostics;
        if (!post_source_quadrature.ok())
        {
            failure = cqmom_failure(std::move(current_moments), post_source_quadrature.quadrature, TransportStatus::inversion_failed,
                                    TransportStage::post_source_inversion, diagnostics, start_time, post_source_quadrature.message);
            return false;
        }
        current_quadrature = post_source_quadrature.quadrature;
        return true;
    };

    if (!execute_transport_sequence(options.scheme, t, time_step, apply_continuous, apply_source))
    {
        return std::move(*failure);
    }

    finalize_diagnostics(diagnostics);
    return CqmomDenseTransportResult{
        std::move(current_moments), std::move(current_quadrature), TransportStatus::success, std::move(diagnostics), t + time_step, ""};
}

CqmomDenseTransportResult advance_transport_cqmom(const RawMomentTensor& moments, const PhaseSpaceHierarchy& hierarchy, const double t,
                                                  const double time_step, const PhaseSpaceVelocity& velocity, const MomentSource& source,
                                                  const CqmomDenseTransportOptions& options)
{
    return advance_transport_cqmom_dense_impl(moments, hierarchy, t, time_step, velocity, source, nullptr, options);
}

CqmomDenseTransportResult advance_transport_cqmom(const RawMomentTensor& moments, const PhaseSpaceHierarchy& hierarchy, const double t,
                                                  const double time_step, const PhaseSpaceVelocity& velocity, const MomentSource& source,
                                                  const MomentBoundaryContribution& boundary, const CqmomDenseTransportOptions& options)
{
    return advance_transport_cqmom_dense_impl(moments, hierarchy, t, time_step, velocity, source, &boundary, options);
}

namespace
{

// Lean and CqmomState transport share this path; it applies the selected C/D
// sequence to the fixed CQMOM moment layout and rebuilds quadrature at stage
// boundaries.
CqmomLeanTransportResult advance_transport_cqmom_lean_impl(const CqmomLeanTransportState& state, const double t, const double time_step,
                                                           const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                           const CqmomMomentBoundaryContribution* const boundary,
                                                           const CqmomLeanTransportOptions& options, const bool adaptive_branch_counts,
                                                           const CqmomDiagnosticLevel diagnostic_level, const VandermondeSolver vandermonde_solver,
                                                           const CqmomMetadataTrust metadata_trust)
{
    internal::performance::record_transport_macro_step();
    TransportDiagnostics diagnostics;
    diagnostics.scheme = options.scheme;
    try
    {
        if (metadata_trust == CqmomMetadataTrust::trusted_case_state)
        {
            validate_trusted_lean_cqmom_input(state, t, time_step, velocity, source, options.scheme);
        }
        else
        {
            validate_lean_cqmom_input(state, t, time_step, velocity, source, options.scheme);
        }
        validate_cqmom_boundary(boundary);
    }
    catch (const std::exception& error)
    {
        return lean_cqmom_failure(state, std::nullopt, TransportStatus::invalid_input, TransportStage::none, diagnostics, t, error.what());
    }

    const auto compute_stage_cqmom = [&](const CqmomMomentData& moments)
    {
        const CqmomLeanCase lean_case{moments,          options.conditional_options, options.qmom_options, adaptive_branch_counts,
                                      diagnostic_level, vandermonde_solver};
        if (metadata_trust == CqmomMetadataTrust::trusted_case_state)
        {
            return internal::compute_trusted_cqmom(lean_case);
        }
        return compute_cqmom(lean_case);
    };

    CqmomResult initial_quadrature = compute_stage_cqmom(state.moments);
    diagnostics.initial_cqmom_status = initial_quadrature.status;
    diagnostics.initial_cqmom_diagnostics = initial_quadrature.diagnostics;
    if (!initial_quadrature.ok())
    {
        return lean_cqmom_failure(state, initial_quadrature.quadrature, TransportStatus::inversion_failed, TransportStage::initial_inversion,
                                  diagnostics, t, initial_quadrature.message);
    }

    CqmomLeanTransportState current_state = state;
    CqmomQuadrature current_quadrature = initial_quadrature.quadrature;
    std::optional<CqmomLeanTransportResult> failure;

    const auto apply_continuous = [&](const double start_time, const double step, const bool final_continuous)
    {
        internal::performance::record_continuous_operator_call();
        const CqmomNodeAdvectionResult continuous = advance_nodes(current_quadrature, start_time, step, velocity, options.node_options);
        record_continuous_diagnostics(diagnostics, final_continuous, continuous.status, continuous.diagnostics);
        if (!continuous.ok())
        {
            failure = lean_cqmom_failure(std::move(current_state), current_quadrature, TransportStatus::continuous_failed,
                                         final_continuous ? TransportStage::second_continuous_half_step : TransportStage::first_continuous_half_step,
                                         diagnostics, start_time, continuous.message);
            return false;
        }

        try
        {
            internal::performance::record_moment_reconstruction_call();
            // Trusted state transport preserves the existing fixed layout and
            // reconstructs only updated moment values after node advection.
            current_state.moments = metadata_trust == CqmomMetadataTrust::trusted_case_state
                                        ? internal::reconstruct_trusted_cqmom_moments(continuous.quadrature, current_state.moments.layout)
                                        : reconstruct_cqmom_moments(continuous.quadrature, current_state.moments.layout);
            if (final_continuous)
            {
                diagnostics.final_reconstruction_diagnostics = successful_reconstruction(current_state.moments.values.size());
            }
            else
            {
                diagnostics.first_reconstruction_diagnostics = successful_reconstruction(current_state.moments.values.size());
            }
            current_quadrature = continuous.quadrature;
        }
        catch (const std::exception& error)
        {
            if (final_continuous)
            {
                diagnostics.final_reconstruction_diagnostics = {false, current_state.moments.values.size(), error.what()};
            }
            else
            {
                diagnostics.first_reconstruction_diagnostics = {false, current_state.moments.values.size(), error.what()};
            }
            failure = lean_cqmom_failure(std::move(current_state), continuous.quadrature, TransportStatus::moment_reconstruction_failed,
                                         final_continuous ? TransportStage::final_moment_reconstruction : TransportStage::first_moment_reconstruction,
                                         diagnostics, start_time + step, error.what());
            return false;
        }
        return true;
    };

    const auto apply_source = [&](const double start_time)
    {
        internal::performance::record_discontinuous_operator_call();
        // D-step: CqmomMomentSource and optional boundary contribution share
        // the existing SSP-RK3 source update path.
        const CqmomMomentSourceUpdateResult source_update =
            boundary == nullptr
                ? advance_cqmom_moments_with_source(current_state.moments, start_time, time_step, source, options.source_options)
                : advance_cqmom_moments_with_source(current_state.moments, start_time, time_step, source, *boundary, options.source_options);
        record_source_update_diagnostics(diagnostics, source_update);
        if (!source_update.ok())
        {
            failure = lean_cqmom_failure(std::move(current_state), current_quadrature, TransportStatus::source_failed,
                                         TransportStage::discontinuous_source_step, diagnostics, start_time, source_update.message);
            return false;
        }

        current_state.moments = source_update.moments;
        // Reconstruct quadrature after the moment update so CQMOM
        // realizability is checked before any later continuous stage.
        internal::performance::record_post_source_inversion();
        CqmomResult post_source_quadrature = compute_stage_cqmom(current_state.moments);
        diagnostics.post_source_cqmom_status = post_source_quadrature.status;
        diagnostics.post_source_cqmom_diagnostics = post_source_quadrature.diagnostics;
        if (!post_source_quadrature.ok())
        {
            failure = lean_cqmom_failure(std::move(current_state), post_source_quadrature.quadrature, TransportStatus::inversion_failed,
                                         TransportStage::post_source_inversion, diagnostics, start_time, post_source_quadrature.message);
            return false;
        }
        current_quadrature = post_source_quadrature.quadrature;
        return true;
    };

    if (!execute_transport_sequence(options.scheme, t, time_step, apply_continuous, apply_source))
    {
        return std::move(*failure);
    }

    finalize_diagnostics(diagnostics);
    return CqmomLeanTransportResult{
        std::move(current_state), std::move(current_quadrature), TransportStatus::success, std::move(diagnostics), t + time_step, ""};
}

} // namespace

static CqmomStateTransportResult advance_transport_cqmom_state_impl(const CqmomState& state, double t, double time_step,
                                                                    const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                                    const CqmomMomentBoundaryContribution* boundary,
                                                                    const CqmomStateTransportOptions& options);

CqmomLeanTransportResult advance_transport_cqmom(const CqmomLeanTransportState& state, const double t, const double time_step,
                                                 const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                 const CqmomLeanTransportOptions& options)
{
    return advance_transport_cqmom_lean_impl(state, t, time_step, velocity, source, nullptr, options, false, CqmomDiagnosticLevel::detailed,
                                             VandermondeSolver::automatic, CqmomMetadataTrust::public_input);
}

CqmomLeanTransportResult advance_transport_cqmom(const CqmomLeanTransportState& state, const double t, const double time_step,
                                                 const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                 const CqmomMomentBoundaryContribution& boundary, const CqmomLeanTransportOptions& options)
{
    return advance_transport_cqmom_lean_impl(state, t, time_step, velocity, source, &boundary, options, false, CqmomDiagnosticLevel::detailed,
                                             VandermondeSolver::automatic, CqmomMetadataTrust::public_input);
}

CqmomStateTransportResult advance_transport_cqmom(const CqmomState& state, const double t, const double time_step, const PhaseSpaceVelocity& velocity,
                                                  const CqmomMomentSource& source, const CqmomStateTransportOptions& options)
{
    return advance_transport_cqmom_state_impl(state, t, time_step, velocity, source, nullptr, options);
}

CqmomStateTransportResult advance_transport_cqmom(const CqmomState& state, const double t, const double time_step, const PhaseSpaceVelocity& velocity,
                                                  const CqmomMomentSource& source, const CqmomMomentBoundaryContribution& boundary,
                                                  const CqmomStateTransportOptions& options)
{
    return advance_transport_cqmom_state_impl(state, t, time_step, velocity, source, &boundary, options);
}

static CqmomStateTransportResult advance_transport_cqmom_state_impl(const CqmomState& state, const double t, const double time_step,
                                                                    const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                                    const CqmomMomentBoundaryContribution* const boundary,
                                                                    const CqmomStateTransportOptions& options)
{
    // State transport reuses the lean numerical sequence and then wraps the
    // updated fixed-layout moments back into the public representation-neutral state.
    const CqmomLeanCase& lean_case = internal::cqmom_state_case(state);
    const CqmomLeanTransportState transport_state{lean_case.moments};
    CqmomLeanTransportOptions lean_options;
    lean_options.conditional_options = lean_case.conditional_options;
    lean_options.qmom_options = lean_case.qmom_options;
    lean_options.node_options = options.node_options;
    lean_options.source_options = options.source_options;
    lean_options.scheme = options.scheme;

    CqmomLeanTransportResult lean_result =
        advance_transport_cqmom_lean_impl(transport_state, t, time_step, velocity, source, boundary, lean_options, lean_case.adaptive_branch_counts,
                                          lean_case.diagnostic_level, lean_case.vandermonde_solver, CqmomMetadataTrust::trusted_case_state);

    CqmomLeanCase updated_case{lean_result.state.moments,        lean_case.conditional_options, lean_case.qmom_options,
                               lean_case.adaptive_branch_counts, lean_case.diagnostic_level,    lean_case.vandermonde_solver};
    TransportDiagnostics diagnostics = neutral_cqmom_diagnostics(std::move(lean_result.diagnostics));
    std::string message = neutral_cqmom_message(std::move(lean_result.message));
    return CqmomStateTransportResult{internal::create_cqmom_state(std::move(updated_case)),
                                     std::move(lean_result.quadrature),
                                     lean_result.status,
                                     std::move(diagnostics),
                                     lean_result.final_time,
                                     std::move(message)};
}

CqmomLeanTransportResult advance_transport_cqmom(const CqmomLeanTransportState& state, const double t, const double time_step,
                                                 const PhaseSpaceVelocity& velocity, const CqmomLeanTransportOptions& options)
{
    const CqmomMomentSource zero_source = [](double, const CqmomSourceContext& moments)
    { return CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
    return advance_transport_cqmom(state, t, time_step, velocity, zero_source, options);
}

} // namespace qmomkit
