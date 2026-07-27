#pragma once

#include "qmomkit/cqmom.hpp"
#include "qmomkit/node_advection.hpp"
#include "qmomkit/qmom.hpp"

#include <optional>
#include <string>
#include <vector>

namespace qmomkit::internal
{

enum class ContinuousMomentRateStatus
{
    success,
    invalid_input,
    inversion_failed,
    velocity_failed
};

struct QmomContinuousMomentRateResult
{
    std::vector<double> rates;
    ContinuousMomentRateStatus status = ContinuousMomentRateStatus::success;
    QmomStatus inversion_status = QmomStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == ContinuousMomentRateStatus::success;
    }
};

struct CqmomContinuousMomentRateResult
{
    std::optional<CqmomMomentData> rates;
    ContinuousMomentRateStatus status = ContinuousMomentRateStatus::success;
    CqmomStatus inversion_status = CqmomStatus::success;
    std::string message;
    std::vector<CqmomDiagnostics> diagnostics;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == ContinuousMomentRateStatus::success;
    }
};

[[nodiscard]] std::vector<std::vector<double>> evaluate_terminal_velocities(const PhaseSpaceVelocity& velocity, double time,
                                                                            const std::vector<std::vector<double>>& terminal_abscissas);

[[nodiscard]] QmomContinuousMomentRateResult evaluate_qmom_continuous_moment_rates(const std::vector<double>& moments, std::size_t node_count,
                                                                                   double time, const PhaseSpaceVelocity& velocity,
                                                                                   const QmomOptions& qmom_options = {});

[[nodiscard]] CqmomContinuousMomentRateResult evaluate_trusted_cqmom_continuous_moment_rates(const CqmomLeanCase& lean_case, double time,
                                                                                             const PhaseSpaceVelocity& velocity);

} // namespace qmomkit::internal
