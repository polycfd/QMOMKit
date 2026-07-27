#include "qmomkit/moment_source.hpp"

#include "internal/cqmom_validation.hpp"
#include "internal/numeric_utils.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
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

void validate_fixed_cqmom_moment_data(const CqmomMomentData& moments, const char* context)
{
    internal::validate_fixed_cqmom_moment_data(moments, {std::string(context) + " requires the fixed CQMOM moment layout.",
                                                         std::string(context) + " values must match the fixed layout size.",
                                                         std::string(context) + " values must be finite."});
}

void validate_rate(const CqmomMomentRate& rate, const std::size_t expected_size)
{
    internal::validate_cqmom_moment_rate(
        rate, expected_size, {"CQMOM source rate values must match the source context size.", "CQMOM source rate values must be finite."});
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
        throw std::invalid_argument("CQMOM source stage values must match the fixed layout size.");
    }
    if (!internal::all_finite(state))
    {
        throw std::invalid_argument("CQMOM source stage values must be finite.");
    }
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
    if (options.ssp_rk3_options.min_time_step < 0.0 ||
        (options.ssp_rk3_options.min_time_step > 0.0 && time_step < options.ssp_rk3_options.min_time_step))
    {
        throw std::invalid_argument("Moment-source update timestep controls are invalid.");
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

void validate_moment_data_request(const CqmomMomentData& moments, const double t, const double time_step, const CqmomMomentSource& source,
                                  const MomentSourceUpdateOptions& options)
{
    validate_fixed_cqmom_moment_data(moments, "CQMOM source update");
    if (!std::isfinite(t) || !std::isfinite(time_step) || !std::isfinite(t + time_step))
    {
        throw std::invalid_argument("CQMOM source update time input contains non-finite values.");
    }
    if (time_step <= 0.0)
    {
        throw std::invalid_argument("CQMOM source update timestep must be positive.");
    }
    if (!source)
    {
        throw std::invalid_argument("CQMOM source callback must be provided.");
    }
    if (options.ssp_rk3_options.min_time_step < 0.0 ||
        (options.ssp_rk3_options.min_time_step > 0.0 && time_step < options.ssp_rk3_options.min_time_step))
    {
        throw std::invalid_argument("CQMOM source update timestep controls are invalid.");
    }
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

[[nodiscard]] MomentSource compose_source_and_boundary(const MomentSource& source, const MomentBoundaryContribution& boundary)
{
    // Boundary exchange is an additive rate inside the discontinuous source
    // operator, so SSP-RK3 still sees a single combined RHS.
    return [source, boundary](const double stage_time, const std::vector<double>& state)
    {
        std::vector<double> source_rate = source(stage_time, state);
        validate_source_rate(source_rate, state.size());
        MomentBoundaryRate boundary_rate;
        try
        {
            boundary_rate = boundary(stage_time, state);
        }
        catch (const std::exception& error)
        {
            throw boundary_callback_error(error);
        }
        validate_boundary_rate(boundary_rate, state.size());
        for (std::size_t i = 0; i < source_rate.size(); ++i)
        {
            source_rate[i] += boundary_rate.values[i];
        }
        return source_rate;
    };
}

[[nodiscard]] CqmomMomentSource compose_cqmom_source_and_boundary(const CqmomMomentSource& source, const CqmomMomentBoundaryContribution& boundary)
{
    // Boundary exchange is an additive rate inside the discontinuous source
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

[[nodiscard]] CqmomMomentSourceUpdateResult advance_validated_cqmom_moments_with_source(const CqmomMomentData& moments, const double t,
                                                                                        const double time_step, const CqmomMomentSource& source,
                                                                                        const MomentSourceUpdateOptions& options)
{
    try
    {
        const CqmomMomentLayout& layout = moments.layout;
        // The public wrapper validated the fixed layout before this helper.
        // SSP-RK3 stages mutate only values, so stages reuse the layout and check
        // vector size/finite values instead of rebuilding CQMOM metadata.
        const auto rhs = [&layout, &source](const double stage_time, const std::vector<double>& state)
        {
            validate_cqmom_stage_state(state, layout);
            const CqmomSourceContext context{layout.hierarchy, layout.required_indices, state};
            CqmomMomentRate rate = source(stage_time, context);
            validate_rate(rate, context.size());
            return rate.values;
        };

        SspRk3Result step = ssp_rk3_step(moments.values, t, time_step, rhs, options.ssp_rk3_options);
        const MomentSourceUpdateStatus status = status_from_ssp_rk3(step.status);
        const MomentSourceUpdateDiagnostics diagnostics = diagnostics_from_ssp_rk3(step, status);
        const std::string message = diagnostics.message;
        return CqmomMomentSourceUpdateResult{CqmomMomentData{layout, std::move(step.state)}, status, diagnostics, message};
    }
    catch (const std::exception& error)
    {
        return failed_moment_data_result(moments, error.what());
    }
}

} // namespace

CqmomSourceContext::CqmomSourceContext(const PhaseSpaceHierarchy& source_hierarchy, const std::vector<MultiIndex>& indices,
                                       const std::vector<double>& values)
    : hierarchy_(&source_hierarchy), indices_(&indices), values_(&values)
{
}

std::size_t CqmomSourceContext::size() const noexcept
{
    return values_->size();
}

const MultiIndex& CqmomSourceContext::index(const std::size_t i) const
{
    return indices_->at(i);
}

double CqmomSourceContext::value(const std::size_t i) const
{
    return values_->at(i);
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

    try
    {
        const SspRk3Result step = ssp_rk3_step(moments, t, time_step, source, options.ssp_rk3_options);
        const MomentSourceUpdateStatus status = status_from_ssp_rk3(step.status);
        const MomentSourceUpdateDiagnostics diagnostics = diagnostics_from_ssp_rk3(step, status);
        const std::string message = diagnostics.message;
        return MomentSourceUpdateResult{step.state, status, diagnostics, message};
    }
    catch (const std::exception& error)
    {
        MomentSourceUpdateDiagnostics diagnostics;
        diagnostics.message = error.what();
        return MomentSourceUpdateResult{moments, MomentSourceUpdateStatus::source_failed, diagnostics, error.what()};
    }
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

} // namespace qmomkit
