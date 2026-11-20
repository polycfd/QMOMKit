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
 * The quadrature family of a QMOM reconstruction with node count N.
 * `classical` is QMOM: up to N nodes from m_0, ..., m_{2N-1}. `hyperbolic` is
 * HyQMOM (Fox & Laurent 2022): up to N + 1 nodes from m_0, ..., m_{2N}, on
 * the real line only. `generalized` is GQMOM (Fox, Laurent & Passalacqua
 * 2023): up to N_G nodes from m_0, ..., m_{2N}, with N_G > N set in
 * `QmomCase::generalized_options`.
 */
enum class QmomFamily
{
    classical,
    hyperbolic,
    generalized
};

/**
 * How GQMOM continues the recurrence coefficients of the moments to N_G nodes.
 * `normal` is Definition 1 of Fox, Laurent & Passalacqua (2023), whose
 * parameter nu controls the tails, `gamma` is Definition 2, `lognormal` is
 * Definition 3, and `beta` is Definition 4. Abscissa bounds [l, u] set the
 * interval of the last three: the number density of `gamma` and `lognormal`
 * starts at l instead of 0, that of `beta` lies in [l, u] instead of [0, 1],
 * and the nodes must lie within the bounds. `normal` takes no bounds.
 */
enum class GeneralizedContinuation
{
    normal,
    gamma,
    lognormal,
    beta
};

/**
 * The settings of the `generalized` family, which no other family reads.
 * `node_count` is N_G, the number of nodes of the continued quadrature. It has
 * no default and must exceed the node count N of the reconstruction. Only the
 * `normal` continuation reads `nu`.
 */
struct GeneralizedOptions
{
    std::size_t node_count = 0;
    GeneralizedContinuation continuation = GeneralizedContinuation::normal;
    double nu = 1.0;
};

/**
 * Numerical controls for one-dimensional QMOM quadrature construction.
 *
 * Adaptive node reduction is on by default. Set `adaptive = false` to require
 * the full quadrature of the family (`QmomCase::family`), N nodes for
 * `classical`, N + 1 for `hyperbolic` and N_G for `generalized`, or, for
 * `generalized` with moments on the boundary, the classical N nodes. `rmin`
 * and `eabs` apply only to classical quadratures, not to the N + 1 nodes of
 * `hyperbolic` or the continued nodes of `generalized`.
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
    /**
     * The number of nodes of the classical quadrature on which the result
     * rests, which never exceeds the node count N of the reconstruction. For
     * `classical`, it equals the number of nodes of the result. For
     * `hyperbolic`, it is N when the result holds the N + 1 nodes of J_{N+1},
     * because these are built from the recurrence of the classical N-node
     * quadrature of m_0, ..., m_{2N-1}. For `generalized`, it is the number
     * k of nodes of the classical quadrature that is continued to N_G nodes.
     * When the moments lie on or below the boundary and a classical
     * quadrature is returned, it is the number of nodes of that quadrature. A
     * failed result has a base node count of 0.
     */
    std::size_t base_node_count = 0;

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

    /**
     * The base node count of the state's quadrature, as defined for
     * `QmomResult::base_node_count`. For `classical`, it equals the active
     * node count. For `hyperbolic`, it is N while the state holds the N + 1
     * nodes of J_{N+1}, and the active node count otherwise. For
     * `generalized`, it is k while the state holds N_G nodes continued from k,
     * and the active node count otherwise.
     */
    [[nodiscard]] std::size_t base_node_count() const noexcept
    {
        return quadrature_.base_node_count;
    }

    /**
     * The maximum node count of the state's quadrature: N for `classical`,
     * N + 1 for `hyperbolic`, with N the node count of its case, and N_G for
     * `generalized`.
     */
    [[nodiscard]] std::size_t maximum_node_count() const noexcept
    {
        switch (family_)
        {
        case QmomFamily::classical:
            return node_count_;
        case QmomFamily::hyperbolic:
            return node_count_ + 1;
        case QmomFamily::generalized:
            return generalized_options_.node_count;
        }
        return node_count_;
    }

    /** The family of the state's case. */
    [[nodiscard]] QmomFamily family() const noexcept
    {
        return family_;
    }

    /** The settings of the generalized family of the state's case. */
    [[nodiscard]] const GeneralizedOptions& generalized_options() const noexcept
    {
        return generalized_options_;
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
    QmomState(std::vector<double> moments, QmomResult quadrature, std::size_t node_count, QmomOptions qmom_options,
              std::optional<AbscissaBounds> abscissa_bounds, QmomFamily family, GeneralizedOptions generalized_options);
    friend class internal::QmomStateAccess;

    std::vector<double> moments_;
    QmomResult quadrature_;
    std::size_t node_count_ = 0;
    QmomOptions qmom_options_;
    std::optional<AbscissaBounds> abscissa_bounds_;
    QmomFamily family_ = QmomFamily::classical;
    GeneralizedOptions generalized_options_;
};

/**
 * Characteristic speeds of the moment system of free transport closed by
 * HyQMOM (Fox & Laurent 2022), for a hyperbolic state with node count n:
 * ascending, 2n + 1 values when its moments give n + 1 nodes, 2k + 1 when they
 * give a classical quadrature of k <= n nodes, none on failure. `weights[i]` is
 * the weight of `speeds[i]`, and `closure_moment` is m_{2n+1}, the flux of
 * m_{2n}. On failure, `weights` is empty and `closure_moment` is zero.
 */
struct HyqmomCharacteristicSpeeds
{
    std::vector<double> speeds;
    std::vector<double> weights;
    double closure_moment = 0.0;
    QmomStatus status = QmomStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == QmomStatus::success;
    }
};

/**
 * The characteristic speeds of the stored moments of `state`: those of the
 * quadrature that `compute_qmom(...)` reconstructs from them with the state's
 * node count n, family and options. For n + 1 nodes, they are the eigenvalues
 * of J_n and K_{n+1}, which are those of the Jacobian of the closed system
 * (shown for n <= 9 by Fox & Laurent 2022, checked by QMOMKit's tests for
 * sample moments up to n = 10); for a classical quadrature of k <= n nodes,
 * its k abscissas, each twice, and their mean. The speeds thus count only the
 * nodes that the moments resolve, which after advection can be fewer than the
 * state holds. The speeds with their weights have the moments m_0, ...,
 * m_{2n+1} of the reconstructed quadrature, and `closure_moment` is its
 * m_{2n+1}. For n + 1 nodes, the weights are (n + 1) / (2n + 1) times the
 * Gauss weights of J_n and n / (2n + 1) times those of K_{n+1}. For a
 * classical quadrature, the two copies of each abscissa share its weight in
 * the ratio k + 1 to k, and the mean has weight zero. `invalid_input` for a
 * family other than `hyperbolic`; otherwise the status and message of a
 * failed reconstruction, such as `not_realizable` in fixed mode when the
 * moments no longer support n + 1 nodes, and `not_realizable` when an
 * eigensolve does not converge or a speed or the closure moment is too large
 * for a double.
 */
[[nodiscard]] HyqmomCharacteristicSpeeds compute_hyqmom_characteristic_speeds(const QmomState& state);

/**
 * Outcome of creating an evolution state from moments.
 *
 * `empty`: the total number density is below the reconstruction cutoff
 * (`QmomOptions::cutoff`), as in an empty cell, so there is nothing to
 * reconstruct. `invalid_input`: the moments cannot form a state with the
 * settings of the case, for example because they are too few or not finite,
 * or the realizability policy is invalid, for example a bounded interval
 * without bounds. `realizability_rejected`: the requested realizability check
 * did not accept the moments. `reconstruction_failed`: no quadrature could be
 * reconstructed from the moments, or the moments of the accepted quadrature
 * exceed the range of a double.
 */
enum class StateCandidateStatus
{
    success,
    empty,
    invalid_input,
    realizability_rejected,
    reconstruction_failed
};

/**
 * A QMOM evolution state created from moments, or the reason none was.
 *
 * `state` holds the state only on success: take it with `*state` and use it
 * like a state from `QmomCase::initial_state()`. `realizability` holds the
 * verdict of the requested realizability check, whether it accepted the
 * moments or not, and `reconstruction_status` the status of the quadrature
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
 * `moments[k] == m_k`. A case with `node_count == N` requires at least
 * `2 * N` moments, or `2 * N + 1` for the `hyperbolic` and `generalized`
 * families. `family` selects the quadrature family, and `generalized_options`
 * holds the settings of the `generalized` family.
 */
struct QmomCase
{
    MomentSet moments;
    std::size_t node_count = 0;
    QmomOptions qmom_options{};
    std::optional<AbscissaBounds> abscissa_bounds{};
    QmomFamily family = QmomFamily::classical;
    GeneralizedOptions generalized_options{};

    /**
     * Create an evolution state from other moments, with the node count,
     * options, abscissa bounds, family and `generalized_options` of this case.
     * The moments of the case are not used. This is how a host rebuilds the
     * state of moments it changed itself, for example by transport in physical
     * space. All moments must be finite, including those the reconstruction
     * does not use, above order 2N - 1, or 2N for `hyperbolic` and
     * `generalized`. The quadrature is reconstructed once, the stored moments
     * are those of the accepted quadrature, and a failure is reported in the
     * result rather than thrown.
     */
    [[nodiscard]] QmomStateCandidate state_candidate(const MomentSet& candidate_moments) const;

    /**
     * Create an evolution state as above, after checking the moments against
     * `policy` unless its `check_mode` is `None`. The check accepts strictly
     * realizable moments, boundary moments when `policy.accept_boundary` is
     * set, and indeterminate moments when `policy.accept_indeterminate` is
     * set. Moments with m0 below the cutoff are accepted whatever the verdict,
     * so that moments without number density are `empty` whatever the
     * policy. An invalid policy, such as a bounded interval without an
     * interval, is invalid input before the moments are checked, and moments
     * too few for the node count are invalid input, checked or not.
     */
    [[nodiscard]] QmomStateCandidate state_candidate(const MomentSet& candidate_moments, const RealizabilityPolicy& policy) const;

    /**
     * Construct the accepted initial phase-space evolution state from the
     * moments of this case, throwing `std::runtime_error` when they cannot be
     * reconstructed. As for `state_candidate`, all moments must be finite,
     * including those the reconstruction does not use.
     */
    [[nodiscard]] QmomState initial_state() const;
};

/**
 * Construct a one-dimensional QMOM quadrature from raw moments.
 *
 * The input moments are raw moments in ascending order. Returned abscissas are
 * sorted in ascending order, and weights use the same ordering. This and the
 * next overload reconstruct the `classical` family. The other families are
 * reconstructed from a `QmomCase`.
 */
[[nodiscard]] QmomResult compute_qmom(const MomentSet& moments, std::size_t node_count, const QmomOptions& options = {});

/**
 * Construct a one-dimensional QMOM quadrature and require returned abscissas
 * to lie in the optional closed interval.
 */
[[nodiscard]] QmomResult compute_qmom(const MomentSet& moments, std::size_t node_count, const QmomOptions& options,
                                      std::optional<AbscissaBounds> abscissa_bounds);

/**
 * Construct a one-dimensional QMOM quadrature from a case object, with its
 * family. The hyperbolic family and the normal continuation of the generalized
 * family take no abscissa bounds: with them, the result is `invalid_input`.
 */
[[nodiscard]] QmomResult compute_qmom(const QmomCase& qmom_case);

} // namespace qmomkit
