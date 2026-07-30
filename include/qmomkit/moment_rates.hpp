#pragma once

#include "qmomkit/cqmom.hpp"
#include "qmomkit/moment_source.hpp"
#include "qmomkit/node_advection.hpp"
#include "qmomkit/qmom.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace qmomkit
{

enum class ContinuousMomentRateStatus
{
    success,
    invalid_input,
    inversion_failed,
    velocity_failed
};

enum class DiscontinuousMomentRateStatus
{
    success,
    invalid_input,
    source_failed,
    invalid_rate
};

/**
 * Instantaneous one-dimensional raw-moment rate from continuous node motion.
 */
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

/**
 * Instantaneous fixed-layout CQMOM raw-moment rate from continuous node motion.
 */
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

/**
 * Instantaneous one-dimensional raw-moment rate from a discontinuous source.
 */
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

/**
 * Instantaneous fixed-layout CQMOM raw-moment rate from a discontinuous source.
 */
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

/**
 * Evaluate dm/dt for one-dimensional raw moments induced by continuous
 * phase-space node motion at one time.
 *
 * This function performs the QMOM inversion, evaluates the velocity callback
 * at terminal nodes, and reconstructs rates for the input moment orders. It
 * does not advance time, retry a step, or expose quadrature nodes.
 */
[[nodiscard]] QmomContinuousMomentRateResult evaluate_qmom_continuous_moment_rates(const std::vector<double>& moments, std::size_t node_count,
                                                                                   double time, const PhaseSpaceVelocity& velocity,
                                                                                   const QmomOptions& options = {});

/**
 * Evaluate dM/dt for fixed-layout CQMOM raw moments induced by continuous
 * phase-space node motion at one time.
 *
 * This function performs the CQMOM inversion, evaluates the velocity callback
 * at terminal nodes, and reconstructs rates in `case_.moments.layout` order.
 * It does not advance time, retry a step, or expose quadrature nodes.
 */
[[nodiscard]] CqmomContinuousMomentRateResult evaluate_cqmom_continuous_moment_rates(const CqmomLeanCase& case_, double time,
                                                                                     const PhaseSpaceVelocity& velocity);

/**
 * Evaluate dm/dt for one-dimensional raw moments by invoking one source
 * callback at one time.
 */
[[nodiscard]] QmomDiscontinuousMomentRateResult evaluate_qmom_discontinuous_moment_rates(double time, const std::vector<double>& moments,
                                                                                         const MomentSource& source);

/**
 * Evaluate dM/dt for fixed-layout CQMOM raw moments by invoking one source
 * callback at one time. Returned values preserve the input moment layout.
 */
[[nodiscard]] CqmomDiscontinuousMomentRateResult evaluate_cqmom_discontinuous_moment_rates(double time, const CqmomMomentData& moments,
                                                                                           const CqmomMomentSource& source);

} // namespace qmomkit
