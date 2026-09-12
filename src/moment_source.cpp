#include "qmomkit/moment_source.hpp"

#include "internal/cqmom_case.hpp"
#include "internal/cqmom_validation.hpp"
#include "internal/moment_source.hpp"
#include "internal/numeric_utils.hpp"
#include "qmomkit/moment_reconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace qmomkit
{
namespace
{

[[nodiscard]] MomentSourceUpdateStatus status_from_ssp_rk3(const SspRk3Status status)
{
    switch (status)
    {
    case SspRk3Status::success:
        return MomentSourceUpdateStatus::success;
    case SspRk3Status::rejected_by_predicate:
        return MomentSourceUpdateStatus::rejected_by_predicate;
    case SspRk3Status::minimum_timestep_reached:
        return MomentSourceUpdateStatus::minimum_timestep_reached;
    }
    return MomentSourceUpdateStatus::source_failed;
}

[[nodiscard]] std::string message_from_status(const MomentSourceUpdateStatus status)
{
    switch (status)
    {
    case MomentSourceUpdateStatus::success:
        return "";
    case MomentSourceUpdateStatus::invalid_input:
        return "Invalid moment-source update input.";
    case MomentSourceUpdateStatus::source_failed:
        return "Moment-source callback failed.";
    case MomentSourceUpdateStatus::rejected_by_predicate:
        return "Moment-source update was rejected by the acceptability predicate.";
    case MomentSourceUpdateStatus::minimum_timestep_reached:
        return "Moment-source update reached the minimum timestep before an acceptable state was found.";
    }
    return "Moment-source update failed.";
}

[[nodiscard]] MomentSourceUpdateDiagnostics diagnostics_from_ssp_rk3(const SspRk3Result& result, const MomentSourceUpdateStatus status)
{
    return MomentSourceUpdateDiagnostics{result.status,
                                         result.status == SspRk3Status::success,
                                         result.retries > 0 || result.status != SspRk3Status::success,
                                         result.retries,
                                         result.time_step_used,
                                         result.error_estimate,
                                         message_from_status(status)};
}

[[nodiscard]] MomentSourceUpdateResult invalid_result(std::vector<double> moments, std::string message)
{
    MomentSourceUpdateDiagnostics diagnostics;
    diagnostics.message = message;
    return MomentSourceUpdateResult{std::move(moments), MomentSourceUpdateStatus::invalid_input, diagnostics, std::move(message)};
}

[[nodiscard]] CqmomMomentSourceUpdateResult invalid_moment_data_result(CqmomMomentData moments, std::string message)
{
    MomentSourceUpdateDiagnostics diagnostics;
    diagnostics.message = message;
    return CqmomMomentSourceUpdateResult{std::move(moments), MomentSourceUpdateStatus::invalid_input, diagnostics, std::move(message)};
}

[[nodiscard]] CqmomMomentSourceUpdateResult failed_moment_data_result(CqmomMomentData moments, std::string message)
{
    MomentSourceUpdateDiagnostics diagnostics;
    diagnostics.message = message;
    return CqmomMomentSourceUpdateResult{std::move(moments), MomentSourceUpdateStatus::source_failed, diagnostics, std::move(message)};
}

void validate_source_context(const PhaseSpaceHierarchy& hierarchy, const std::vector<MultiIndex>& indices, const std::vector<double>& values)
{
    if (indices.size() != values.size())
    {
        throw std::invalid_argument("CQMOM source context indices and values must have the same size.");
    }
    for (const MultiIndex& index : indices)
    {
        if (index.size() != hierarchy.number_of_dimensions())
        {
            throw std::invalid_argument("CQMOM source context moment indices must match the hierarchy dimension count.");
        }
    }
}

void validate_source_generalized_order(const double order, const char* const context)
{
    if (!std::isfinite(order))
    {
        throw std::invalid_argument(std::string(context) + " must be finite.");
    }
    if (order < 0.0)
    {
        throw std::invalid_argument(std::string(context) + " must be non-negative.");
    }
}

void validate_fixed_cqmom_moment_data(const CqmomMomentData& moments, const char* context)
{
    internal::validate_fixed_cqmom_moment_data(moments, {std::string(context) + " requires the fixed CQMOM moment layout.",
                                                         std::string(context) + " values must match the fixed layout size.",
                                                         std::string(context) + " values must be finite."});
}

void validate_rate(const CqmomMomentRate& rate, const std::size_t expected_size)
{
    internal::validate_cqmom_moment_rate(
        rate, expected_size,
        {"CQMOM moment-source rate values must match the source context size.", "CQMOM moment-source rate values must be finite."});
}

void validate_source_rate(const std::vector<double>& rate, const std::size_t expected_size)
{
    if (rate.size() != expected_size)
    {
        throw std::invalid_argument("Moment-source callback rate values must match the moment state size.");
    }
    if (!internal::all_finite(rate))
    {
        throw std::invalid_argument("Moment-source callback rate values must be finite.");
    }
}

void validate_boundary_rate(const MomentBoundaryRate& rate, const std::size_t expected_size)
{
    if (rate.values.size() != expected_size)
    {
        throw std::invalid_argument("Moment-boundary contribution rate values must match the moment state size.");
    }
    if (!internal::all_finite(rate.values))
    {
        throw std::invalid_argument("Moment-boundary contribution rate values must be finite.");
    }
}

void validate_cqmom_boundary_rate(const CqmomMomentRate& rate, const std::size_t expected_size)
{
    internal::validate_cqmom_moment_rate(
        rate, expected_size,
        {"CQMOM boundary contribution rate values must match the source context size.", "CQMOM boundary contribution rate values must be finite."});
}

void validate_cqmom_stage_state(const std::vector<double>& state, const CqmomMomentLayout& layout)
{
    if (state.size() != layout.required_indices.size())
    {
        throw std::invalid_argument("CQMOM moment-source stage values must match the fixed layout size.");
    }
    if (!internal::all_finite(state))
    {
        throw std::invalid_argument("CQMOM moment-source stage values must be finite.");
    }
}

void validate_coordinate_realizability_policies(const CqmomMomentLayout& layout, const std::vector<std::optional<RealizabilityPolicy>>& policies)
{
    const std::size_t dimensions = layout.hierarchy.number_of_dimensions();
    for (std::size_t coordinate = dimensions; coordinate < policies.size(); ++coordinate)
    {
        if (policies[coordinate].has_value())
        {
            throw std::invalid_argument("CQMOM coordinate realizability policy index is outside the phase-space dimension.");
        }
    }
}

[[nodiscard]] bool has_coordinate_realizability_policy(const std::vector<std::optional<RealizabilityPolicy>>& policies,
                                                       const RealizabilityCheckMode check_mode)
{
    return std::any_of(policies.begin(), policies.end(), [check_mode](const std::optional<RealizabilityPolicy>& policy)
                       { return policy.has_value() && policy->check_mode == check_mode; });
}

[[nodiscard]] std::vector<double> extract_coordinate_marginal(const std::vector<double>& state, const CqmomMomentLayout& layout,
                                                              const std::size_t coordinate)
{
    std::vector<std::pair<std::size_t, std::size_t>> orders_and_positions;
    for (std::size_t position = 0; position < layout.required_indices.size(); ++position)
    {
        const MultiIndex& index = layout.required_indices[position];
        bool pure_coordinate = true;
        for (std::size_t dimension = 0; dimension < index.size(); ++dimension)
        {
            if (dimension != coordinate && index[dimension] != 0)
            {
                pure_coordinate = false;
                break;
            }
        }
        if (pure_coordinate)
        {
            orders_and_positions.emplace_back(index[coordinate], position);
        }
    }

    std::sort(orders_and_positions.begin(), orders_and_positions.end());
    std::vector<double> marginal;
    marginal.reserve(orders_and_positions.size());
    for (const auto& order_and_position : orders_and_positions)
    {
        marginal.push_back(state[order_and_position.second]);
    }
    return marginal;
}

[[nodiscard]] bool accepts_coordinate_realizability_policies(const std::vector<double>& state, const CqmomMomentLayout& layout,
                                                             const std::vector<std::optional<RealizabilityPolicy>>& policies,
                                                             const RealizabilityCheckMode check_mode)
{
    for (std::size_t coordinate = 0; coordinate < policies.size(); ++coordinate)
    {
        const std::optional<RealizabilityPolicy>& policy = policies[coordinate];
        if (!policy.has_value() || policy->check_mode != check_mode)
        {
            continue;
        }

        const std::vector<double> marginal = extract_coordinate_marginal(state, layout, coordinate);
        if (!accepts_realizability_policy(*policy, marginal, [](const std::vector<double>& moments) { return MomentSet(moments); }))
        {
            return false;
        }
    }
    return true;
}

[[nodiscard]] std::invalid_argument boundary_callback_error(const std::exception& error)
{
    return std::invalid_argument(std::string("Moment-boundary contribution callback failed: ") + error.what());
}

[[nodiscard]] std::invalid_argument cqmom_boundary_callback_error(const std::exception& error)
{
    return std::invalid_argument(std::string("CQMOM boundary contribution callback failed: ") + error.what());
}

void validate_request(const std::vector<double>& moments, const double t, const double time_step, const MomentSource& source,
                      const MomentSourceUpdateOptions& options)
{
    if (moments.empty())
    {
        throw std::invalid_argument("Moment-source update requires at least one moment.");
    }
    if (!internal::all_finite(moments) || !std::isfinite(t) || !std::isfinite(time_step))
    {
        throw std::invalid_argument("Moment-source update input contains non-finite values.");
    }
    if (time_step <= 0.0)
    {
        throw std::invalid_argument("Moment-source update timestep must be positive.");
    }
    if (!source)
    {
        throw std::invalid_argument("Moment-source callback must be provided.");
    }
    if (!std::isfinite(options.ssp_rk3_options.min_time_step) || options.ssp_rk3_options.min_time_step < 0.0 ||
        (options.ssp_rk3_options.min_time_step > 0.0 && time_step < options.ssp_rk3_options.min_time_step))
    {
        throw std::invalid_argument("Moment-source update timestep controls are invalid.");
    }
    if (std::any_of(options.coordinate_realizability_policies.begin(), options.coordinate_realizability_policies.end(),
                    [](const std::optional<RealizabilityPolicy>& policy) { return policy.has_value(); }))
    {
        throw std::invalid_argument("QMOM moment-source update does not support coordinate_realizability_policies; this option is CQMOM-only.");
    }
}

void validate_request(const std::vector<double>& moments, const double t, const double time_step, const MomentSource& source,
                      const MomentBoundaryContribution& boundary, const MomentSourceUpdateOptions& options)
{
    validate_request(moments, t, time_step, source, options);
    if (!boundary)
    {
        throw std::invalid_argument("Moment-boundary contribution callback must be provided.");
    }
}

void validate_cqmom_source_update_controls(const double t, const double time_step, const CqmomMomentSource& source,
                                           const MomentSourceUpdateOptions& options)
{
    if (!std::isfinite(t) || !std::isfinite(time_step) || !std::isfinite(t + time_step))
    {
        throw std::invalid_argument("CQMOM moment-source update time input contains non-finite values.");
    }
    if (time_step <= 0.0)
    {
        throw std::invalid_argument("CQMOM moment-source update timestep must be positive.");
    }
    if (!source)
    {
        throw std::invalid_argument("CQMOM source callback must be provided.");
    }
    if (!std::isfinite(options.ssp_rk3_options.min_time_step) || options.ssp_rk3_options.min_time_step < 0.0 ||
        (options.ssp_rk3_options.min_time_step > 0.0 && time_step < options.ssp_rk3_options.min_time_step))
    {
        throw std::invalid_argument("CQMOM moment-source update timestep controls are invalid.");
    }
}

void validate_moment_data_request(const CqmomMomentData& moments, const double t, const double time_step, const CqmomMomentSource& source,
                                  const MomentSourceUpdateOptions& options)
{
    validate_fixed_cqmom_moment_data(moments, "CQMOM moment-source update");
    validate_cqmom_source_update_controls(t, time_step, source, options);
    validate_coordinate_realizability_policies(moments.layout, options.coordinate_realizability_policies);
}

void validate_trusted_moment_data_request(const CqmomMomentData& moments, const double t, const double time_step, const CqmomMomentSource& source,
                                          const MomentSourceUpdateOptions& options)
{
    internal::validate_cqmom_moment_values_match_layout(
        moments, {"", "CQMOM moment-source update values must match the fixed layout size.", "CQMOM moment-source update values must be finite."});
    validate_cqmom_source_update_controls(t, time_step, source, options);
    validate_coordinate_realizability_policies(moments.layout, options.coordinate_realizability_policies);
}

void validate_moment_data_request(const CqmomMomentData& moments, const double t, const double time_step, const CqmomMomentSource& source,
                                  const CqmomMomentBoundaryContribution& boundary, const MomentSourceUpdateOptions& options)
{
    validate_moment_data_request(moments, t, time_step, source, options);
    if (!boundary)
    {
        throw std::invalid_argument("CQMOM boundary contribution callback must be provided.");
    }
}

void validate_trusted_moment_data_request(const CqmomMomentData& moments, const double t, const double time_step, const CqmomMomentSource& source,
                                          const CqmomMomentBoundaryContribution& boundary, const MomentSourceUpdateOptions& options)
{
    validate_trusted_moment_data_request(moments, t, time_step, source, options);
    if (!boundary)
    {
        throw std::invalid_argument("CQMOM boundary contribution callback must be provided.");
    }
}

[[nodiscard]] MomentSource compose_source_and_boundary(const MomentSource& source, const MomentBoundaryContribution& boundary)
{
    // Boundary exchange is an additive rate inside the moment-source
    // operator, so SSP-RK3 still sees a single combined RHS.
    return [source, boundary](const double stage_time, const QmomSourceContext& context)
    {
        std::vector<double> source_rate = source(stage_time, context);
        validate_source_rate(source_rate, context.size());
        MomentBoundaryRate boundary_rate;
        try
        {
            boundary_rate = boundary(stage_time, context.values());
        }
        catch (const std::exception& error)
        {
            throw boundary_callback_error(error);
        }
        validate_boundary_rate(boundary_rate, context.size());
        for (std::size_t i = 0; i < source_rate.size(); ++i)
        {
            source_rate[i] += boundary_rate.values[i];
        }
        return source_rate;
    };
}

[[nodiscard]] MomentSourceUpdateResult advance_moments_with_source(const std::vector<double>& moments, const double t, const double time_step,
                                                                   const OdeRhs& source, const MomentSourceUpdateOptions& options)
{
    try
    {
        std::vector<double> state = moments;
        double current_time = t;
        double remaining_time = time_step;
        std::size_t retries = 0;
        bool rejected = false;
        MomentSourceUpdateDiagnostics diagnostics;

        while (true)
        {
            const SspRk3Result step = ssp_rk3_step(state, current_time, remaining_time, source, options.ssp_rk3_options);
            const MomentSourceUpdateStatus status = status_from_ssp_rk3(step.status);
            diagnostics = diagnostics_from_ssp_rk3(step, status);
            retries += step.retries;
            rejected = rejected || diagnostics.rejected;
            diagnostics.retries = retries;
            diagnostics.rejected = rejected;

            if (status != MomentSourceUpdateStatus::success)
            {
                return MomentSourceUpdateResult{moments, status, diagnostics, diagnostics.message};
            }

            // Positive steps smaller than the remaining interval's spacing do not complete the interval.
            const double next_remaining_time = remaining_time - step.time_step_used;
            if (!std::isfinite(step.time_step_used) || step.time_step_used <= 0.0 || !(next_remaining_time < remaining_time))
            {
                diagnostics.message = "Moment-source accepted timestep did not make representable interval progress.";
                return MomentSourceUpdateResult{moments, MomentSourceUpdateStatus::source_failed, diagnostics, diagnostics.message};
            }
            if (step.time_step_used == remaining_time)
            {
                return MomentSourceUpdateResult{step.state, status, diagnostics, diagnostics.message};
            }

            state = step.state;
            current_time += step.time_step_used;
            remaining_time = next_remaining_time;
        }
    }
    catch (const std::exception& error)
    {
        MomentSourceUpdateDiagnostics diagnostics;
        diagnostics.message = error.what();
        return MomentSourceUpdateResult{moments, MomentSourceUpdateStatus::source_failed, diagnostics, error.what()};
    }
}

[[nodiscard]] CqmomMomentSource compose_cqmom_source_and_boundary(const CqmomMomentSource& source, const CqmomMomentBoundaryContribution& boundary)
{
    // Boundary exchange is an additive rate inside the moment-source
    // operator, so SSP-RK3 still sees a single combined RHS.
    return [source, boundary](const double stage_time, const CqmomSourceContext& context)
    {
        CqmomMomentRate source_rate = source(stage_time, context);
        validate_rate(source_rate, context.size());
        CqmomMomentRate boundary_rate;
        try
        {
            boundary_rate = boundary(stage_time, context);
        }
        catch (const std::exception& error)
        {
            throw cqmom_boundary_callback_error(error);
        }
        validate_cqmom_boundary_rate(boundary_rate, context.size());
        for (std::size_t i = 0; i < source_rate.values.size(); ++i)
        {
            source_rate.values[i] += boundary_rate.values[i];
        }
        return source_rate;
    };
}

[[nodiscard]] CqmomMomentSourceUpdateResult
advance_validated_cqmom_moments_with_source(const CqmomMomentData& moments, const double t, const double time_step, const CqmomMomentSource& source,
                                            const MomentSourceUpdateOptions& options, const CqmomLeanCase* const reconstruction_case = nullptr,
                                            const std::vector<CqmomActiveBranchCount>* const active_branch_counts = nullptr,
                                            const internal::PreparedFixedCqmomMetadata* const prepared_fixed_metadata = nullptr)
{
    try
    {
        const CqmomMomentLayout& layout = moments.layout;
        // The public wrapper validated the fixed layout before this helper.
        // SSP-RK3 stages mutate only values, so stages reuse the layout and check
        // vector size/finite values instead of rebuilding CQMOM metadata.
        const auto rhs = [&layout, &source, reconstruction_case, active_branch_counts, prepared_fixed_metadata](const double stage_time,
                                                                                                                const std::vector<double>& state)
        {
            validate_cqmom_stage_state(state, layout);
            const CqmomSourceContext context =
                reconstruction_case == nullptr
                    ? CqmomSourceContext{stage_time, layout.hierarchy, layout.required_indices, state}
                    : internal::CqmomSourceContextAccess::create(stage_time, layout.hierarchy, layout.required_indices, state, *reconstruction_case,
                                                                 *active_branch_counts, prepared_fixed_metadata);
            CqmomMomentRate rate = source(stage_time, context);
            validate_rate(rate, context.size());
            return rate.values;
        };

        SspRk3Options ssp_rk3_options = options.ssp_rk3_options;
        const std::vector<std::optional<RealizabilityPolicy>>& coordinate_policies = options.coordinate_realizability_policies;
        if (has_coordinate_realizability_policy(coordinate_policies, RealizabilityCheckMode::FinalOnly))
        {
            const StateAcceptancePredicate existing_accept = ssp_rk3_options.accept;
            ssp_rk3_options.accept = [&layout, &coordinate_policies, existing_accept](const std::vector<double>& state)
            {
                return (!existing_accept || existing_accept(state)) &&
                       accepts_coordinate_realizability_policies(state, layout, coordinate_policies, RealizabilityCheckMode::FinalOnly);
            };
        }
        if (has_coordinate_realizability_policy(coordinate_policies, RealizabilityCheckMode::EveryStage))
        {
            const StageAcceptancePredicate existing_stage_accept = ssp_rk3_options.stage_accept;
            ssp_rk3_options.stage_accept =
                [&layout, &coordinate_policies, existing_stage_accept](const SspRk3Stage stage, const std::vector<double>& state)
            {
                return (!existing_stage_accept || existing_stage_accept(stage, state)) &&
                       accepts_coordinate_realizability_policies(state, layout, coordinate_policies, RealizabilityCheckMode::EveryStage);
            };
        }

        MomentSourceUpdateOptions composed_options = options;
        composed_options.ssp_rk3_options = std::move(ssp_rk3_options);
        MomentSourceUpdateResult update = advance_moments_with_source(moments.values, t, time_step, rhs, composed_options);
        return CqmomMomentSourceUpdateResult{CqmomMomentData{layout, std::move(update.moments)}, update.status, std::move(update.diagnostics),
                                             std::move(update.message)};
    }
    catch (const std::exception& error)
    {
        return failed_moment_data_result(moments, error.what());
    }
}

} // namespace

namespace internal
{

QmomSourceContext QmomSourceContextAccess::create(const double time, const std::vector<double>& values, const std::size_t reconstruction_node_count,
                                                  const QmomOptions& reconstruction_options, std::optional<AbscissaBounds> reconstruction_bounds)
{
    return QmomSourceContext{time, values, reconstruction_node_count, reconstruction_options, std::move(reconstruction_bounds)};
}

CqmomSourceContext CqmomSourceContextAccess::create(const double time, const PhaseSpaceHierarchy& hierarchy, const std::vector<MultiIndex>& indices,
                                                    const std::vector<double>& values, const CqmomLeanCase& reconstruction_case,
                                                    const std::vector<CqmomActiveBranchCount>& active_branch_counts,
                                                    const PreparedFixedCqmomMetadata* const prepared_fixed_metadata)
{
    return CqmomSourceContext{time, hierarchy, indices, values, reconstruction_case, active_branch_counts, prepared_fixed_metadata};
}

} // namespace internal

QmomSourceContext::QmomSourceContext(const double time, const std::vector<double>& values) : time_(time), values_(&values)
{
    if (!std::isfinite(time_))
    {
        throw std::invalid_argument("QMOM source context time must be finite.");
    }
}

QmomSourceContext::QmomSourceContext(const double time, const std::vector<double>& values, const std::size_t reconstruction_node_count,
                                     const QmomOptions& reconstruction_options, std::optional<AbscissaBounds> reconstruction_bounds)
    : QmomSourceContext(time, values)
{
    if (reconstruction_node_count == 0)
    {
        throw std::invalid_argument("QMOM source-stage quadrature reconstruction requires at least one node.");
    }
    reconstruction_node_count_ = reconstruction_node_count;
    reconstruction_options_ = &reconstruction_options;
    reconstruction_bounds_ = std::move(reconstruction_bounds);
}

double QmomSourceContext::time() const noexcept
{
    return time_;
}

std::size_t QmomSourceContext::moment_count() const noexcept
{
    return values_->size();
}

std::size_t QmomSourceContext::size() const noexcept
{
    return moment_count();
}

const std::vector<double>& QmomSourceContext::values() const noexcept
{
    return *values_;
}

double QmomSourceContext::value(const std::size_t order) const
{
    return values_->at(order);
}

double QmomSourceContext::moment(const std::size_t order) const
{
    if (order < values_->size())
    {
        return (*values_)[order];
    }
    if (reconstruction_options_ == nullptr)
    {
        throw std::out_of_range("Requested QMOM source moment is not tracked by this context.");
    }

    const QmomResult& quadrature = reconstructed_stage_quadrature();
    const std::vector<double> moments =
        reconstruct_raw_moments(Quadrature1D{quadrature.abscissas, quadrature.weights, QuadratureStatus::success, ""}, order + 1);
    return moments.at(order);
}

double QmomSourceContext::generalized_moment(const double order) const
{
    validate_source_generalized_order(order, "Generalized QMOM source moment order");
    return evaluate_generalized_moment(reconstructed_stage_quadrature(), order);
}

const QmomResult& QmomSourceContext::reconstructed_stage_quadrature() const
{
    if (reconstruction_options_ == nullptr)
    {
        throw std::out_of_range("Requested QMOM source generalized moment is not available without reconstruction metadata.");
    }
    if (!reconstructed_stage_quadrature_.has_value())
    {
        QmomResult quadrature = compute_qmom(MomentSet(*values_), reconstruction_node_count_, *reconstruction_options_, reconstruction_bounds_);
        if (!quadrature.ok())
        {
            throw std::runtime_error("QMOM source-stage quadrature reconstruction failed: " + quadrature.message);
        }
        reconstructed_stage_quadrature_ = std::move(quadrature);
    }
    return *reconstructed_stage_quadrature_;
}

CqmomSourceContext::CqmomSourceContext(const double time, const PhaseSpaceHierarchy& source_hierarchy, const std::vector<MultiIndex>& indices,
                                       const std::vector<double>& values)
    : time_(time), hierarchy_(&source_hierarchy), indices_(&indices), values_(&values)
{
    if (!std::isfinite(time_))
    {
        throw std::invalid_argument("CQMOM source context time must be finite.");
    }
    validate_source_context(source_hierarchy, indices, values);
}

CqmomSourceContext::CqmomSourceContext(const double time, const PhaseSpaceHierarchy& source_hierarchy, const std::vector<MultiIndex>& indices,
                                       const std::vector<double>& values, const CqmomLeanCase& reconstruction_case,
                                       const std::vector<CqmomActiveBranchCount>& active_branch_counts,
                                       const internal::PreparedFixedCqmomMetadata* const prepared_fixed_metadata)
    : time_(time), hierarchy_(&source_hierarchy), indices_(&indices), values_(&values), reconstruction_case_(&reconstruction_case),
      active_branch_counts_(&active_branch_counts), prepared_fixed_metadata_(prepared_fixed_metadata)
{
    if (!std::isfinite(time_))
    {
        throw std::invalid_argument("CQMOM source context time must be finite.");
    }
    validate_source_context(source_hierarchy, indices, values);
}

double CqmomSourceContext::time() const noexcept
{
    return time_;
}

std::size_t CqmomSourceContext::moment_count() const noexcept
{
    return values_->size();
}

std::size_t CqmomSourceContext::dimension_count() const noexcept
{
    return hierarchy_->number_of_dimensions();
}

std::size_t CqmomSourceContext::size() const noexcept
{
    return moment_count();
}

const std::vector<MultiIndex>& CqmomSourceContext::indices() const noexcept
{
    return *indices_;
}

const std::vector<double>& CqmomSourceContext::values() const noexcept
{
    return *values_;
}

const MultiIndex& CqmomSourceContext::index(const std::size_t i) const
{
    return indices_->at(i);
}

double CqmomSourceContext::value(const std::size_t i) const
{
    return values_->at(i);
}

std::optional<std::size_t> CqmomSourceContext::position(const MultiIndex& moment_index) const
{
    const auto found = std::find(indices_->begin(), indices_->end(), moment_index);
    if (found == indices_->end())
    {
        return std::nullopt;
    }
    return static_cast<std::size_t>(std::distance(indices_->begin(), found));
}

std::optional<double> CqmomSourceContext::value(const MultiIndex& moment_index) const
{
    const std::optional<std::size_t> entry = position(moment_index);
    if (!entry.has_value())
    {
        return std::nullopt;
    }
    return value(*entry);
}

double CqmomSourceContext::moment(const MultiIndex& moment_index) const
{
    if (const std::optional<double> stored_value = value(moment_index); stored_value.has_value())
    {
        return *stored_value;
    }
    if (reconstruction_case_ == nullptr || active_branch_counts_ == nullptr)
    {
        throw std::out_of_range("Requested CQMOM source moment is not tracked by this context.");
    }

    const SelectedRawMomentResult result = reconstruct_raw_moments(reconstructed_stage_quadrature(), RawMomentRequest{{moment_index}});
    return result.values.front();
}

double CqmomSourceContext::generalized_moment(const std::vector<double>& orders) const
{
    if (orders.size() != dimension_count())
    {
        throw std::invalid_argument("Generalized CQMOM source moment order count must match the phase-space dimension count.");
    }
    for (const double order : orders)
    {
        validate_source_generalized_order(order, "Generalized CQMOM source moment orders");
    }
    return evaluate_generalized_moment(reconstructed_stage_quadrature(), orders);
}

void CqmomSourceContext::for_each_terminal_node(const std::function<void(double, const std::vector<double>&)>& callback) const
{
    if (!callback)
    {
        throw std::invalid_argument("CQMOM terminal-node visitor callback must be provided.");
    }

    const CqmomQuadrature& quadrature = reconstructed_stage_quadrature();
    const std::vector<double>& weights = quadrature.weights();
    const std::vector<std::vector<double>>& abscissas = quadrature.abscissas();
    std::vector<double> coordinates(abscissas.size(), 0.0);
    for (std::size_t node = 0; node < weights.size(); ++node)
    {
        for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
        {
            coordinates[dimension] = abscissas[dimension][node];
        }
        callback(weights[node], coordinates);
    }
}

const CqmomQuadrature& CqmomSourceContext::reconstructed_stage_quadrature() const
{
    if (reconstruction_case_ == nullptr || active_branch_counts_ == nullptr)
    {
        throw std::out_of_range("CQMOM source-stage quadrature is unavailable for this context.");
    }

    if (!reconstructed_stage_quadrature_.has_value())
    {
        const CqmomLeanCase stage_case{CqmomMomentData{reconstruction_case_->moments.layout, *values_},
                                       reconstruction_case_->conditional_options,
                                       reconstruction_case_->qmom_options,
                                       reconstruction_case_->adaptive_branch_counts,
                                       reconstruction_case_->diagnostic_level,
                                       reconstruction_case_->vandermonde_solver,
                                       reconstruction_case_->coordinate_support};
        CqmomResult quadrature = internal::compute_cqmom_internal(stage_case, *active_branch_counts_, stage_case.adaptive_branch_counts,
                                                                  prepared_fixed_metadata_, !stage_case.adaptive_branch_counts);
        if (!quadrature.ok())
        {
            throw std::runtime_error("CQMOM source-stage quadrature reconstruction failed: " + quadrature.message);
        }
        reconstructed_stage_quadrature_ = std::move(quadrature.quadrature);
    }
    return *reconstructed_stage_quadrature_;
}

const PhaseSpaceHierarchy& CqmomSourceContext::hierarchy() const noexcept
{
    return *hierarchy_;
}

MomentSourceUpdateResult advance_moments_with_source(const std::vector<double>& moments, const double t, const double time_step,
                                                     const MomentSource& source, const MomentSourceUpdateOptions& options)
{
    try
    {
        validate_request(moments, t, time_step, source, options);
    }
    catch (const std::exception& error)
    {
        return invalid_result(moments, error.what());
    }

    const OdeRhs rhs = [&source](const double stage_time, const std::vector<double>& state)
    {
        const QmomSourceContext context{stage_time, state};
        return source(stage_time, context);
    };
    return advance_moments_with_source(moments, t, time_step, rhs, options);
}

MomentSourceUpdateResult advance_moments_with_source(const std::vector<double>& moments, const double t, const double time_step,
                                                     const MomentSource& source, const MomentBoundaryContribution& boundary,
                                                     const MomentSourceUpdateOptions& options)
{
    try
    {
        validate_request(moments, t, time_step, source, boundary, options);
    }
    catch (const std::exception& error)
    {
        return invalid_result(moments, error.what());
    }

    return advance_moments_with_source(moments, t, time_step, compose_source_and_boundary(source, boundary), options);
}

CqmomMomentSourceUpdateResult advance_cqmom_moments_with_source(const CqmomMomentData& moments, const double t, const double time_step,
                                                                const CqmomMomentSource& source, const MomentSourceUpdateOptions& options)
{
    try
    {
        validate_moment_data_request(moments, t, time_step, source, options);
    }
    catch (const std::exception& error)
    {
        return invalid_moment_data_result(moments, error.what());
    }

    return advance_validated_cqmom_moments_with_source(moments, t, time_step, source, options);
}

CqmomMomentSourceUpdateResult advance_cqmom_moments_with_source(const CqmomMomentData& moments, const double t, const double time_step,
                                                                const CqmomMomentSource& source, const CqmomMomentBoundaryContribution& boundary,
                                                                const MomentSourceUpdateOptions& options)
{
    try
    {
        validate_moment_data_request(moments, t, time_step, source, boundary, options);
    }
    catch (const std::exception& error)
    {
        return invalid_moment_data_result(moments, error.what());
    }

    return advance_validated_cqmom_moments_with_source(moments, t, time_step, compose_cqmom_source_and_boundary(source, boundary), options);
}

namespace internal
{

MomentSourceUpdateResult advance_trusted_qmom_moments_with_source(const std::vector<double>& moments, const double t, const double time_step,
                                                                  const MomentSource& source, const MomentSourceUpdateOptions& options,
                                                                  const std::size_t reconstruction_node_count,
                                                                  const QmomOptions& reconstruction_options,
                                                                  const std::optional<AbscissaBounds> reconstruction_bounds)
{
    try
    {
        validate_request(moments, t, time_step, source, options);
    }
    catch (const std::exception& error)
    {
        return invalid_result(moments, error.what());
    }

    const OdeRhs rhs = [&source, reconstruction_node_count, &reconstruction_options, reconstruction_bounds](const double stage_time,
                                                                                                            const std::vector<double>& state)
    {
        const QmomSourceContext context =
            QmomSourceContextAccess::create(stage_time, state, reconstruction_node_count, reconstruction_options, reconstruction_bounds);
        return source(stage_time, context);
    };
    return advance_moments_with_source(moments, t, time_step, rhs, options);
}

MomentSourceUpdateResult advance_trusted_qmom_moments_with_source(const std::vector<double>& moments, const double t, const double time_step,
                                                                  const MomentSource& source, const MomentBoundaryContribution& boundary,
                                                                  const MomentSourceUpdateOptions& options,
                                                                  const std::size_t reconstruction_node_count,
                                                                  const QmomOptions& reconstruction_options,
                                                                  const std::optional<AbscissaBounds> reconstruction_bounds)
{
    try
    {
        validate_request(moments, t, time_step, source, boundary, options);
    }
    catch (const std::exception& error)
    {
        return invalid_result(moments, error.what());
    }

    const MomentSource composed_source = compose_source_and_boundary(source, boundary);
    const OdeRhs rhs = [&composed_source, reconstruction_node_count, &reconstruction_options, reconstruction_bounds](const double stage_time,
                                                                                                                     const std::vector<double>& state)
    {
        const QmomSourceContext context =
            QmomSourceContextAccess::create(stage_time, state, reconstruction_node_count, reconstruction_options, reconstruction_bounds);
        return composed_source(stage_time, context);
    };
    return advance_moments_with_source(moments, t, time_step, rhs, options);
}

CqmomMomentSourceUpdateResult advance_trusted_cqmom_moments_with_source(const CqmomMomentData& moments, const double t, const double time_step,
                                                                        const CqmomMomentSource& source, const MomentSourceUpdateOptions& options)
{
    try
    {
        validate_trusted_moment_data_request(moments, t, time_step, source, options);
    }
    catch (const std::exception& error)
    {
        return invalid_moment_data_result(moments, error.what());
    }

    return advance_validated_cqmom_moments_with_source(moments, t, time_step, source, options);
}

CqmomMomentSourceUpdateResult advance_trusted_cqmom_moments_with_source(const CqmomMomentData& moments, const double t, const double time_step,
                                                                        const CqmomMomentSource& source, const MomentSourceUpdateOptions& options,
                                                                        const CqmomLeanCase& reconstruction_case,
                                                                        const std::vector<CqmomActiveBranchCount>& active_branch_counts,
                                                                        const PreparedFixedCqmomMetadata* const prepared_fixed_metadata)
{
    try
    {
        validate_trusted_moment_data_request(moments, t, time_step, source, options);
    }
    catch (const std::exception& error)
    {
        return invalid_moment_data_result(moments, error.what());
    }

    return advance_validated_cqmom_moments_with_source(moments, t, time_step, source, options, &reconstruction_case, &active_branch_counts,
                                                       prepared_fixed_metadata);
}

CqmomMomentSourceUpdateResult advance_trusted_cqmom_moments_with_source(const CqmomMomentData& moments, const double t, const double time_step,
                                                                        const CqmomMomentSource& source,
                                                                        const CqmomMomentBoundaryContribution& boundary,
                                                                        const MomentSourceUpdateOptions& options)
{
    try
    {
        validate_trusted_moment_data_request(moments, t, time_step, source, boundary, options);
    }
    catch (const std::exception& error)
    {
        return invalid_moment_data_result(moments, error.what());
    }

    return advance_validated_cqmom_moments_with_source(moments, t, time_step, compose_cqmom_source_and_boundary(source, boundary), options);
}

CqmomMomentSourceUpdateResult
advance_trusted_cqmom_moments_with_source(const CqmomMomentData& moments, const double t, const double time_step, const CqmomMomentSource& source,
                                          const CqmomMomentBoundaryContribution& boundary, const MomentSourceUpdateOptions& options,
                                          const CqmomLeanCase& reconstruction_case, const std::vector<CqmomActiveBranchCount>& active_branch_counts,
                                          const PreparedFixedCqmomMetadata* const prepared_fixed_metadata)
{
    try
    {
        validate_trusted_moment_data_request(moments, t, time_step, source, boundary, options);
    }
    catch (const std::exception& error)
    {
        return invalid_moment_data_result(moments, error.what());
    }

    return advance_validated_cqmom_moments_with_source(moments, t, time_step, compose_cqmom_source_and_boundary(source, boundary), options,
                                                       &reconstruction_case, &active_branch_counts, prepared_fixed_metadata);
}

} // namespace internal

} // namespace qmomkit
