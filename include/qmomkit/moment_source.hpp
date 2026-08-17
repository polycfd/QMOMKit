#pragma once

#include "qmomkit/cqmom.hpp"
#include "qmomkit/realizability_policy.hpp"
#include "qmomkit/runge_kutta.hpp"

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace qmomkit
{
namespace internal
{
class CqmomSourceContextAccess;
struct PreparedFixedCqmomMetadata;
} // namespace internal

/**
 * Generic moment-source callback.
 *
 * The state vector is an arbitrary raw-moment vector owned by the caller. The
 * returned vector must have the same size and represents dm/dt at the supplied
 * time. The library assigns no physical meaning to any component.
 */
using MomentSource = std::function<std::vector<double>(double, const std::vector<double>&)>;

/**
 * One-dimensional raw-moment rate from phase-space boundary exchange.
 *
 * Boundary contributions are additive rates for raw moments caused by
 * nucleation, injection, removal, evaporation/condensation exchange, or similar
 * phase-space boundary physics. They are not a generic PDE boundary-condition
 * framework.
 */
struct MomentBoundaryRate
{
    std::vector<double> values;
};

/**
 * User callback for phase-space boundary exchange contributions to raw-moment
 * rates.
 *
 * The callback receives the current time and the current raw-moment vector. It
 * returns one additive rate per moment entry, in the same order as the input
 * moment vector.
 */
using MomentBoundaryContribution = std::function<MomentBoundaryRate(double, const std::vector<double>&)>;

/**
 * Non-owning mathematical view of the current CQMOM source state.
 *
 * Entries are ordered consistently with the derivative values expected from a
 * `CqmomMomentSource`. The context exposes the CQMOM moment layout selected for
 * phase-space evolution, along with phase-space hierarchy, stage time, moment multi-indices,
 * and moment values. The time is explicit and matches the stage or evaluation
 * time passed as the first argument to the callback. That layout is independent
 * of the underlying dense or lean CQMOM moment storage backend.
 * Phase-space source contexts may cache a quadrature reconstructed from their
 * current stage values while evaluating an untracked requested moment.
 */
class CqmomSourceContext
{
public:
    CqmomSourceContext(double time, const PhaseSpaceHierarchy& hierarchy, const std::vector<MultiIndex>& indices, const std::vector<double>& values);

    [[nodiscard]] double time() const noexcept;
    [[nodiscard]] std::size_t moment_count() const noexcept;
    [[nodiscard]] std::size_t dimension_count() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] const std::vector<MultiIndex>& indices() const noexcept;
    [[nodiscard]] const std::vector<double>& values() const noexcept;
    [[nodiscard]] const MultiIndex& index(std::size_t i) const;
    [[nodiscard]] double value(std::size_t i) const;
    [[nodiscard]] std::optional<std::size_t> position(const MultiIndex& moment_index) const;
    [[nodiscard]] std::optional<double> value(const MultiIndex& moment_index) const;
    /**
     * Return the current raw moment for `moment_index`.
     *
     * A tracked entry returns the current source-stage value directly. During
     * CQMOM phase-space source evolution, an untracked entry is evaluated from
     * a quadrature reconstructed from the current SSP-RK stage. Other context
     * uses reject untracked requests.
     */
    [[nodiscard]] double moment(const MultiIndex& moment_index) const;
    [[nodiscard]] const PhaseSpaceHierarchy& hierarchy() const noexcept;

private:
    friend class internal::CqmomSourceContextAccess;

    CqmomSourceContext(double time, const PhaseSpaceHierarchy& hierarchy, const std::vector<MultiIndex>& indices, const std::vector<double>& values,
                       const CqmomLeanCase& reconstruction_case, const std::vector<CqmomActiveBranchCount>& active_branch_counts,
                       const internal::PreparedFixedCqmomMetadata* prepared_fixed_metadata);

    double time_ = 0.0;
    const PhaseSpaceHierarchy* hierarchy_;
    const std::vector<MultiIndex>* indices_;
    const std::vector<double>* values_;
    const CqmomLeanCase* reconstruction_case_ = nullptr;
    const std::vector<CqmomActiveBranchCount>* active_branch_counts_ = nullptr;
    const internal::PreparedFixedCqmomMetadata* prepared_fixed_metadata_ = nullptr;
    mutable std::optional<CqmomQuadrature> reconstructed_stage_quadrature_;
};

/**
 * CQMOM required-moment derivative.
 *
 * `values[i]` stores dM_alpha/dt for the moment entry exposed by
 * `CqmomSourceContext::index(i)` and `CqmomSourceContext::value(i)`.
 */
struct CqmomMomentRate
{
    std::vector<double> values;
};

/**
 * CQMOM-native source callback over the phase-space evolution CQMOM moment layout.
 *
 * The callback receives the current CQMOM source context and returns one
 * derivative value per context entry, in context order.
 */
using CqmomMomentSource = std::function<CqmomMomentRate(double, const CqmomSourceContext&)>;

/**
 * CQMOM-native phase-space boundary exchange contribution.
 *
 * The callback receives the same representation-neutral CQMOM context used by
 * `CqmomMomentSource` and returns one additive boundary rate per context entry,
 * in context order.
 */
using CqmomMomentBoundaryContribution = std::function<CqmomMomentRate(double, const CqmomSourceContext&)>;

enum class MomentSourceUpdateStatus
{
    success,
    invalid_input,
    source_failed,
    rejected_by_predicate,
    minimum_timestep_reached
};

struct MomentSourceUpdateOptions
{
    SspRk3Options ssp_rk3_options;
    std::vector<std::optional<RealizabilityPolicy>> coordinate_realizability_policies;
};

struct MomentSourceUpdateDiagnostics
{
    SspRk3Status ssp_rk3_status = SspRk3Status::success;
    bool accepted = false;
    bool rejected = false;
    std::size_t retries = 0;
    double time_step_used = 0.0;
    std::optional<double> error_estimate;
    std::string message;
};

struct MomentSourceUpdateResult
{
    std::vector<double> moments;
    MomentSourceUpdateStatus status = MomentSourceUpdateStatus::success;
    MomentSourceUpdateDiagnostics diagnostics;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == MomentSourceUpdateStatus::success;
    }
};

struct CqmomMomentSourceUpdateResult
{
    CqmomMomentData moments;
    MomentSourceUpdateStatus status = MomentSourceUpdateStatus::success;
    MomentSourceUpdateDiagnostics diagnostics;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == MomentSourceUpdateStatus::success;
    }
};

/**
 * Advanced composable operation: advance a raw-moment vector through a
 * moment-source operator using the existing SSP-RK3 integrator.
 * Optional acceptability checks are supplied through
 * `options.ssp_rk3_options.accept`; rejected candidates use the SSP-RK3 retry
 * policy with timestep halving.
 */
[[nodiscard]] MomentSourceUpdateResult advance_moments_with_source(const std::vector<double>& moments, double t, double time_step,
                                                                   const MomentSource& source, const MomentSourceUpdateOptions& options = {});

/**
 * Advanced composable operation: advance raw moments through a generic moment source
 * plus an additive phase-space boundary exchange contribution using the same
 * SSP-RK3 moment-source update path.
 */
[[nodiscard]] MomentSourceUpdateResult advance_moments_with_source(const std::vector<double>& moments, double t, double time_step,
                                                                   const MomentSource& source, const MomentBoundaryContribution& boundary,
                                                                   const MomentSourceUpdateOptions& options);

/**
 * Advanced composable operation: advance fixed-layout CQMOM moments through a
 * moment-source operator using the existing vector SSP-RK3 integrator. The RK stages
 * are exposed to the CQMOM source callback through `CqmomSourceContext`.
 */
[[nodiscard]] CqmomMomentSourceUpdateResult advance_cqmom_moments_with_source(const CqmomMomentData& moments, double t, double time_step,
                                                                              const CqmomMomentSource& source,
                                                                              const MomentSourceUpdateOptions& options = {});

/**
 * Advanced composable operation: advance fixed-layout CQMOM moments through a
 * generic moment source plus an additive phase-space boundary exchange contribution
 * using the same SSP-RK3 source-update path.
 */
[[nodiscard]] CqmomMomentSourceUpdateResult advance_cqmom_moments_with_source(const CqmomMomentData& moments, double t, double time_step,
                                                                              const CqmomMomentSource& source,
                                                                              const CqmomMomentBoundaryContribution& boundary,
                                                                              const MomentSourceUpdateOptions& options);

} // namespace qmomkit
