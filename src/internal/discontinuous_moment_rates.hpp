#pragma once

#include "qmomkit/moment_source.hpp"

#include <optional>
#include <string>
#include <vector>

namespace qmomkit::internal
{

enum class DiscontinuousMomentRateStatus
{
    success,
    invalid_input,
    source_failed,
    invalid_rate
};

struct QmomDiscontinuousMomentRateResult
{
    std::vector<double> rates;
    DiscontinuousMomentRateStatus status = DiscontinuousMomentRateStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == DiscontinuousMomentRateStatus::success;
    }
};

struct CqmomDiscontinuousMomentRateResult
{
    std::optional<CqmomMomentData> rates;
    DiscontinuousMomentRateStatus status = DiscontinuousMomentRateStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == DiscontinuousMomentRateStatus::success;
    }
};

[[nodiscard]] QmomDiscontinuousMomentRateResult evaluate_qmom_discontinuous_moment_rates(double time, const std::vector<double>& moments,
                                                                                         const MomentSource& source);

[[nodiscard]] QmomDiscontinuousMomentRateResult evaluate_qmom_discontinuous_moment_rates(double time, const std::vector<double>& moments,
                                                                                         const MomentSource& source,
                                                                                         const MomentBoundaryContribution& boundary);

[[nodiscard]] CqmomDiscontinuousMomentRateResult evaluate_cqmom_discontinuous_moment_rates(double time, const CqmomMomentData& moments,
                                                                                           const CqmomMomentSource& source);

} // namespace qmomkit::internal
