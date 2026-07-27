#include "internal/discontinuous_moment_rates.hpp"

#include "internal/cqmom_validation.hpp"
#include "internal/numeric_utils.hpp"

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
