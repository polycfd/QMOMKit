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

enum class AdvectionMomentRateStatus
{
    success,
    invalid_input,
    inversion_failed,
    velocity_failed
};

enum class MomentSourceRateStatus
{
    success,
    invalid_input,
    source_failed,
    invalid_rate
};

/**
 * Instantaneous one-dimensional advection-induced moment rate.
 */
struct QmomAdvectionMomentRateResult
{
    std::vector<double> rates;
    AdvectionMomentRateStatus status = AdvectionMomentRateStatus::success;
    QmomStatus inversion_status = QmomStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == AdvectionMomentRateStatus::success;
    }
};

/**
 * Instantaneous fixed-layout CQMOM advection-induced moment rate.
 */
struct CqmomAdvectionMomentRateResult
{
    std::optional<CqmomMomentData> rates;
    AdvectionMomentRateStatus status = AdvectionMomentRateStatus::success;
    CqmomStatus inversion_status = CqmomStatus::success;
    std::string message;
    std::vector<CqmomDiagnostics> diagnostics;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == AdvectionMomentRateStatus::success;
    }
};

/**
 * Instantaneous one-dimensional moment-source rate.
 */
struct QmomMomentSourceRateResult
{
    std::vector<double> rates;
    MomentSourceRateStatus status = MomentSourceRateStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == MomentSourceRateStatus::success;
    }
};

/**
 * Instantaneous fixed-layout CQMOM moment-source rate.
 */
struct CqmomMomentSourceRateResult
{
    std::optional<CqmomMomentData> rates;
    MomentSourceRateStatus status = MomentSourceRateStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == MomentSourceRateStatus::success;
    }
};

/**
 * Evaluate an advection-induced moment rate for one-dimensional raw moments at
 * one time.
 *
 * This function performs the QMOM inversion, evaluates the velocity callback
 * at terminal nodes, and reconstructs rates for the input moment orders. It
 * does not advance time, retry a step, or expose quadrature nodes.
 */
[[nodiscard]] QmomAdvectionMomentRateResult evaluate_qmom_advection_moment_rates(const std::vector<double>& moments, std::size_t node_count,
                                                                                 double time, const PhaseSpaceVelocity& velocity,
                                                                                 const QmomOptions& options = {});

/**
 * Evaluate an advection-induced moment rate for fixed-layout CQMOM raw
 * moments at one time.
 *
 * This function performs the CQMOM inversion, evaluates the velocity callback
 * at terminal nodes, and reconstructs rates in `case_.moments.layout` order.
 * It does not advance time, retry a step, or expose quadrature nodes.
 */
[[nodiscard]] CqmomAdvectionMomentRateResult evaluate_cqmom_advection_moment_rates(const CqmomLeanCase& case_, double time,
                                                                                   const PhaseSpaceVelocity& velocity);

/**
 * Evaluate a moment-source rate for one-dimensional raw moments by invoking a
 * moment-source callback at one time.
 */
[[nodiscard]] QmomMomentSourceRateResult evaluate_qmom_moment_source_rates(double time, const std::vector<double>& moments,
                                                                           const MomentSource& source);

/**
 * Evaluate dm/dt for one-dimensional raw moments from a source and an additive
 * phase-space boundary contribution at one time. The source callback is
 * invoked before the boundary callback; both receive the supplied time and
 * moment vector. Returned rates preserve the input moment ordering. This
 * operation performs no normalization, reconstruction, or time integration,
 * and does not modify `moments`. Both callbacks must return one finite rate per
 * input moment.
 */
[[nodiscard]] QmomMomentSourceRateResult evaluate_qmom_moment_source_rates(double time, const std::vector<double>& moments,
                                                                           const MomentSource& source, const MomentBoundaryContribution& boundary);

/**
 * Evaluate a moment-source rate for fixed-layout CQMOM raw moments by invoking
 * a moment-source callback at one time. Returned values preserve the input
 * moment layout.
 */
[[nodiscard]] CqmomMomentSourceRateResult evaluate_cqmom_moment_source_rates(double time, const CqmomMomentData& moments,
                                                                             const CqmomMomentSource& source);

} // namespace qmomkit
