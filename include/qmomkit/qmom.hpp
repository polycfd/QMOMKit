#pragma once

#include "qmomkit/abscissa_bounds.hpp"
#include "qmomkit/moment_set.hpp"
#include "qmomkit/realizability.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace qmomkit
{

enum class QmomStatus
{
    success,
    invalid_input,
    not_realizable
};

/**
 * Numerical controls for one-dimensional QMOM quadrature construction.
 *
 * Adaptive node reduction is on by default; set `adaptive = false` to require
 * the requested node count.
 */
struct QmomOptions
{
    bool adaptive = true;
    double rmin = 1.0e-10;
    double eabs = 1.0e-10;
    double cutoff = 1.0e-30;
};

struct QmomResult
{
    std::vector<double> abscissas;
    std::vector<double> weights;
    QmomStatus status = QmomStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == QmomStatus::success;
    }
};

class QmomState;
struct RealizabilityPolicy;

namespace internal
{
class QmomStateAccess;
}

/**
 * Accepted one-dimensional QMOM phase-space evolution state.
 *
 * The stored moments are reconstructed from the accepted quadrature at
 * construction. The state retains the reconstruction policy inherited from
 * its case. Users can inspect the state but cannot mutate moments and
 * quadrature independently.
 */
class QmomState
{
public:
    QmomState(const QmomState&) = default;
    QmomState& operator=(const QmomState&) = default;
    QmomState(QmomState&&) noexcept = default;
    QmomState& operator=(QmomState&&) noexcept = default;
    ~QmomState() = default;

    [[nodiscard]] const std::vector<double>& moment_values() const noexcept
    {
        return moments_;
    }

    [[nodiscard]] const QmomResult& quadrature() const noexcept
    {
        return quadrature_;
    }

    [[nodiscard]] std::size_t active_node_count() const noexcept
    {
        return quadrature_.weights.size();
    }

    [[nodiscard]] std::size_t maximum_node_count() const noexcept
    {
        return maximum_node_count_;
    }

    [[nodiscard]] const QmomOptions& qmom_options() const noexcept
    {
        return qmom_options_;
    }

    [[nodiscard]] const std::optional<AbscissaBounds>& abscissa_bounds() const noexcept
    {
        return abscissa_bounds_;
    }

private:
    QmomState(std::vector<double> moments, QmomResult quadrature, std::size_t maximum_node_count, QmomOptions qmom_options,
              std::optional<AbscissaBounds> abscissa_bounds);
    friend class internal::QmomStateAccess;

    std::vector<double> moments_;
    QmomResult quadrature_;
    std::size_t maximum_node_count_ = 0;
    QmomOptions qmom_options_;
    std::optional<AbscissaBounds> abscissa_bounds_;
};

/**
 * Outcome of creating an evolution state from moments.
 *
 * `invalid_input`: the moments cannot form a state with the settings of the
 * case, for example because they are too few or not finite.
 * `realizability_rejected`: the requested realizability check did not accept
 * the moments. `reconstruction_failed`: no quadrature could be reconstructed
 * from the moments.
 */
enum class StateCandidateStatus
{
    success,
    invalid_input,
    realizability_rejected,
    reconstruction_failed
};

/**
 * A QMOM evolution state created from moments, or the reason none was.
 *
 * `state` holds the state only on success. `realizability` holds the verdict
 * of the requested realizability check, whether it accepted the moments or
 * not, and `reconstruction_status` the status of the quadrature
 * reconstruction, when it ran.
 */
struct QmomStateCandidate
{
    StateCandidateStatus status = StateCandidateStatus::invalid_input;
    std::optional<QmomState> state;
    std::optional<RealizabilityResult> realizability;
    std::optional<QmomStatus> reconstruction_status;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == StateCandidateStatus::success;
    }
};

/**
 * User-facing case for one-dimensional QMOM quadrature construction.
 *
 * `moments` stores unnormalized raw moments in ascending order:
 * `moments[k] == m_k`. A case with `node_count == n` requires at least
 * `2 * n` moments.
 */
struct QmomCase
{
    MomentSet moments;
    std::size_t node_count = 0;
    QmomOptions qmom_options;
    std::optional<AbscissaBounds> abscissa_bounds;

    /**
     * Create an evolution state from other moments, with the node count,
     * options and abscissa bounds of this case; the moments of the case are
     * not used. This is how a host rebuilds the state of moments it changed
     * itself, for example by transport in physical space. The quadrature is
     * reconstructed once, the stored moments are those of the accepted
     * quadrature, and a failure is reported in the result rather than thrown.
     */
    [[nodiscard]] QmomStateCandidate state_candidate(const MomentSet& candidate_moments) const;

    /**
     * Create an evolution state as above, after checking the moments against
     * `policy` unless its `check_mode` is `None`. The check accepts strictly
     * realizable moments, and boundary moments when `policy.accept_boundary`
     * is set.
     */
    [[nodiscard]] QmomStateCandidate state_candidate(const MomentSet& candidate_moments, const RealizabilityPolicy& policy) const;

    /**
     * Construct the accepted initial phase-space evolution state from the
     * moments of this case: `state_candidate(moments)`, which throws
     * `std::runtime_error` when no state can be created.
     */
    [[nodiscard]] QmomState initial_state() const;
};

/**
 * Construct a one-dimensional QMOM quadrature from raw moments.
 *
 * The input moments are raw moments in ascending order. Returned abscissas are
 * sorted in ascending order, and weights use the same ordering.
 */
[[nodiscard]] QmomResult compute_qmom(const MomentSet& moments, std::size_t node_count, const QmomOptions& options = {});

/**
 * Construct a one-dimensional QMOM quadrature and require returned abscissas
 * to lie in the optional closed interval.
 */
[[nodiscard]] QmomResult compute_qmom(const MomentSet& moments, std::size_t node_count, const QmomOptions& options,
                                      std::optional<AbscissaBounds> abscissa_bounds);

/**
 * Construct a one-dimensional QMOM quadrature from a case object.
 */
[[nodiscard]] QmomResult compute_qmom(const QmomCase& qmom_case);

} // namespace qmomkit
