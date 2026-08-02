#include "internal/discontinuous_moment_rates.hpp"

#include "internal/cqmom_validation.hpp"
#include "internal/numeric_utils.hpp"
#include "qmomkit/moment_rates.hpp"

#include <cmath>
#include <exception>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace qmomkit::internal
{
namespace
{

[[nodiscard]] QmomDiscontinuousMomentRateResult qmom_invalid_input(std::string message)
{
    return QmomDiscontinuousMomentRateResult{{}, DiscontinuousMomentRateStatus::invalid_input, std::move(message)};
}

[[nodiscard]] QmomDiscontinuousMomentRateResult qmom_source_failed(std::string message)
{
    return QmomDiscontinuousMomentRateResult{{}, DiscontinuousMomentRateStatus::source_failed, std::move(message)};
}

[[nodiscard]] QmomDiscontinuousMomentRateResult qmom_invalid_rate(std::string message)
{
    return QmomDiscontinuousMomentRateResult{{}, DiscontinuousMomentRateStatus::invalid_rate, std::move(message)};
}

[[nodiscard]] CqmomDiscontinuousMomentRateResult cqmom_invalid_input(std::string message)
{
    return CqmomDiscontinuousMomentRateResult{std::nullopt, DiscontinuousMomentRateStatus::invalid_input, std::move(message)};
}

[[nodiscard]] CqmomDiscontinuousMomentRateResult cqmom_source_failed(std::string message)
{
    return CqmomDiscontinuousMomentRateResult{std::nullopt, DiscontinuousMomentRateStatus::source_failed, std::move(message)};
}

[[nodiscard]] CqmomDiscontinuousMomentRateResult cqmom_invalid_rate(std::string message)
{
    return CqmomDiscontinuousMomentRateResult{std::nullopt, DiscontinuousMomentRateStatus::invalid_rate, std::move(message)};
}

[[nodiscard]] std::string qmom_invalid_rate_message(const std::vector<double>& rate, const std::size_t expected_size)
{
    if (rate.size() != expected_size)
    {
        return "Moment-source callback rate values must match the moment state size.";
    }
    if (!all_finite(rate))
    {
        return "Moment-source callback rate values must be finite.";
    }
    return "";
}

[[nodiscard]] std::string qmom_boundary_invalid_rate_message(const MomentBoundaryRate& rate, const std::size_t expected_size)
{
    if (rate.values.size() != expected_size)
    {
        return "Moment-boundary contribution rate values must match the moment state size.";
    }
    if (!all_finite(rate.values))
    {
        return "Moment-boundary contribution rate values must be finite.";
    }
    return "";
}

void validate_cqmom_stage_state(const CqmomMomentData& moments)
{
    validate_cqmom_moment_values_match_layout(
        moments, {"", "CQMOM source stage values must match the fixed layout size.", "CQMOM source stage values must be finite."});
}

} // namespace

QmomDiscontinuousMomentRateResult evaluate_qmom_discontinuous_moment_rates(const double time, const std::vector<double>& moments,
                                                                           const MomentSource& source)
{
    if (moments.empty())
    {
        return qmom_invalid_input("Moment-source rate evaluation requires at least one moment.");
    }
    if (!std::isfinite(time))
    {
        return qmom_invalid_input("Moment-source rate evaluation time must be finite.");
    }
    if (!source)
    {
        return qmom_invalid_input("Moment-source callback must be provided.");
    }

    std::vector<double> rate;
    try
    {
        rate = source(time, moments);
    }
    catch (const std::exception& error)
    {
        return qmom_source_failed(error.what());
    }

    const std::string invalid_message = qmom_invalid_rate_message(rate, moments.size());
    if (!invalid_message.empty())
    {
        return qmom_invalid_rate(invalid_message);
    }
    return QmomDiscontinuousMomentRateResult{std::move(rate), DiscontinuousMomentRateStatus::success, ""};
}

QmomDiscontinuousMomentRateResult evaluate_qmom_discontinuous_moment_rates(const double time, const std::vector<double>& moments,
                                                                           const MomentSource& source, const MomentBoundaryContribution& boundary)
{
    // Preserve the source-only evaluator's validation precedence while
    // rejecting a missing boundary before either callback is invoked.
    if (!moments.empty() && std::isfinite(time) && source && !boundary)
    {
        return qmom_invalid_input("Moment-boundary contribution callback must be provided.");
    }

    QmomDiscontinuousMomentRateResult source_result = evaluate_qmom_discontinuous_moment_rates(time, moments, source);
    if (!source_result.ok())
    {
        return source_result;
    }

    MomentBoundaryRate boundary_rate;
    try
    {
        boundary_rate = boundary(time, moments);
    }
    catch (const std::exception& error)
    {
        return qmom_source_failed(std::string("Moment-boundary contribution callback failed: ") + error.what());
    }

    const std::string invalid_message = qmom_boundary_invalid_rate_message(boundary_rate, moments.size());
    if (!invalid_message.empty())
    {
        return qmom_invalid_rate(invalid_message);
    }

    for (std::size_t i = 0; i < source_result.rates.size(); ++i)
    {
        source_result.rates[i] += boundary_rate.values[i];
    }
    if (!all_finite(source_result.rates))
    {
        return qmom_invalid_rate("Combined discontinuous moment rate values must be finite.");
    }
    return source_result;
}

CqmomDiscontinuousMomentRateResult evaluate_cqmom_discontinuous_moment_rates(const double time, const CqmomMomentData& moments,
                                                                             const CqmomMomentSource& source)
{
    if (!std::isfinite(time))
    {
        return cqmom_invalid_input("CQMOM source rate evaluation time must be finite.");
    }
    if (!source)
    {
        return cqmom_invalid_input("CQMOM source callback must be provided.");
    }
    try
    {
        validate_cqmom_stage_state(moments);
    }
    catch (const std::exception& error)
    {
        return cqmom_invalid_input(error.what());
    }

    const CqmomSourceContext context{moments.layout.hierarchy, moments.layout.required_indices, moments.values};
    CqmomMomentRate rate;
    try
    {
        rate = source(time, context);
    }
    catch (const std::exception& error)
    {
        return cqmom_source_failed(error.what());
    }

    try
    {
        validate_cqmom_moment_rate(rate, context.size(),
                                   {"CQMOM source rate values must match the source context size.", "CQMOM source rate values must be finite."});
    }
    catch (const std::exception& error)
    {
        return cqmom_invalid_rate(error.what());
    }

    return CqmomDiscontinuousMomentRateResult{CqmomMomentData{moments.layout, std::move(rate.values)}, DiscontinuousMomentRateStatus::success, ""};
}

} // namespace qmomkit::internal

namespace qmomkit
{
namespace
{

[[nodiscard]] DiscontinuousMomentRateStatus to_public_status(const internal::DiscontinuousMomentRateStatus status)
{
    switch (status)
    {
    case internal::DiscontinuousMomentRateStatus::success:
        return DiscontinuousMomentRateStatus::success;
    case internal::DiscontinuousMomentRateStatus::invalid_input:
        return DiscontinuousMomentRateStatus::invalid_input;
    case internal::DiscontinuousMomentRateStatus::source_failed:
        return DiscontinuousMomentRateStatus::source_failed;
    case internal::DiscontinuousMomentRateStatus::invalid_rate:
        return DiscontinuousMomentRateStatus::invalid_rate;
    }
    return DiscontinuousMomentRateStatus::invalid_input;
}

} // namespace

QmomDiscontinuousMomentRateResult evaluate_qmom_discontinuous_moment_rates(const double time, const std::vector<double>& moments,
                                                                           const MomentSource& source)
{
    internal::QmomDiscontinuousMomentRateResult result = internal::evaluate_qmom_discontinuous_moment_rates(time, moments, source);
    return QmomDiscontinuousMomentRateResult{std::move(result.rates), to_public_status(result.status), std::move(result.message)};
}

QmomDiscontinuousMomentRateResult evaluate_qmom_discontinuous_moment_rates(const double time, const std::vector<double>& moments,
                                                                           const MomentSource& source, const MomentBoundaryContribution& boundary)
{
    internal::QmomDiscontinuousMomentRateResult result = internal::evaluate_qmom_discontinuous_moment_rates(time, moments, source, boundary);
    return QmomDiscontinuousMomentRateResult{std::move(result.rates), to_public_status(result.status), std::move(result.message)};
}

CqmomDiscontinuousMomentRateResult evaluate_cqmom_discontinuous_moment_rates(const double time, const CqmomMomentData& moments,
                                                                             const CqmomMomentSource& source)
{
    internal::CqmomDiscontinuousMomentRateResult result = internal::evaluate_cqmom_discontinuous_moment_rates(time, moments, source);
    return CqmomDiscontinuousMomentRateResult{std::move(result.rates), to_public_status(result.status), std::move(result.message)};
}

} // namespace qmomkit
