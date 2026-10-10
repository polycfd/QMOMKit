#pragma once

#include "qmomkit/abscissa_bounds.hpp"
#include "qmomkit/conditional_moments.hpp"
#include "qmomkit/conditioning_path.hpp"
#include "qmomkit/phase_space_hierarchy.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/raw_moment_tensor.hpp"
#include "qmomkit/vandermonde_solver.hpp"

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace qmomkit
{
class CqmomPreparedCase;
class CqmomQuadrature;
class CqmomState;
struct CqmomResult;
struct CqmomCase;
struct SelectedRawMomentResult;

namespace internal
{
class CqmomConstructionAccess;
} // namespace internal

/**
 * Multidimensional quadrature produced by the public CQMOM API.
 *
 * Fixed results use phase-space-order rectangular storage. Adaptive
 * branch-specific results expose terminal weights, phase-space-order terminal
 * abscissas, and terminal conditioning paths in deterministic traversal order.
 */
class CqmomQuadrature
{
public:
    // Defined here rather than in the library so that callers can inline them.
    // Reconstruction reads these once per node per moment, where an out-of-line
    // call for a stored value is a measurable share of the work.
    [[nodiscard]] const PhaseSpaceHierarchy& hierarchy() const noexcept
    {
        return hierarchy_;
    }
    [[nodiscard]] std::size_t number_of_dimensions() const noexcept
    {
        return hierarchy_.number_of_dimensions();
    }
    [[nodiscard]] std::size_t total_nodes() const noexcept
    {
        return weights_.size();
    }
    [[nodiscard]] const std::vector<double>& weights() const noexcept
    {
        return weights_;
    }
    [[nodiscard]] const std::vector<std::vector<double>>& abscissas() const noexcept
    {
        return abscissas_;
    }
    [[nodiscard]] bool adaptive_branch_counts_enabled() const noexcept
    {
        return adaptive_branch_counts_enabled_;
    }
    [[nodiscard]] bool has_rectangular_storage() const noexcept
    {
        return rectangular_storage_;
    }
    [[nodiscard]] const ConditioningPathTable& terminal_paths() const noexcept
    {
        return terminal_paths_;
    }

    [[nodiscard]] double weight(MultiIndexView node_index) const;
    [[nodiscard]] double abscissa(std::size_t dimension, MultiIndexView node_index) const;
    [[nodiscard]] std::vector<double> node(MultiIndexView node_index) const;

private:
    friend CqmomResult compute_cqmom(const CqmomCase& cqmom_case);
    friend class internal::CqmomConstructionAccess;

    CqmomQuadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights, std::vector<std::vector<double>> abscissas,
                    bool adaptive_branch_counts_enabled, bool rectangular_storage, ConditioningPathTable terminal_paths);

    PhaseSpaceHierarchy hierarchy_;
    std::vector<double> weights_;
    std::vector<std::vector<double>> abscissas_;
    ConditioningPathTable terminal_paths_;
    bool adaptive_branch_counts_enabled_ = false;
    bool rectangular_storage_ = true;
};

enum class CqmomStatus
{
    success,
    invalid_input,
    qmom_failed,
    conditional_moment_failed
};

enum class CqmomDiagnosticLevel
{
    none,
    basic,
    detailed
};

using CqmomAbscissaBounds = std::vector<std::optional<AbscissaBounds>>;

[[nodiscard]] std::string format_cqmom_settings(bool adaptive, VandermondeSolver vandermonde_solver, CqmomDiagnosticLevel diagnostic_level);

struct CqmomDiagnostics
{
    std::size_t conditioning_level = 0;
    ConditioningPath parent_path;
    std::size_t phase_space_dimension = 0;
    std::size_t requested_node_count = 0;
    std::size_t actual_node_count = 0;
    bool adaptive_branch_counts_enabled = false;
    bool adaptive_node_count_reduced = false;
    QmomStatus qmom_status = QmomStatus::success;
    ConditionalMomentStatus conditional_moment_status = ConditionalMomentStatus::success;
    ConditionalMomentDiagnostics conditional_moment_diagnostics;
    VandermondeSolver conditional_vandermonde_requested_solver = VandermondeSolver::automatic;
    VandermondeSolver conditional_vandermonde_selected_solver = VandermondeSolver::dense;
    std::string message;
    bool bounded_support_checked = false;
    bool bounded_support_interval_present = false;
    AbscissaBounds bounded_support_interval{0.0, 0.0};
    std::optional<std::size_t> bounded_previous_active_node_count;
    std::size_t bounded_candidate_node_count = 0;
    std::string bounded_candidate_attempt_type;
    std::size_t bounded_accepted_fallback_node_count = 0;
    std::vector<double> bounded_returned_nodes;
    std::vector<double> bounded_returned_weights;
};

/**
 * Current active node count for one recursive CQMOM branch.
 *
 * `configured_node_count` is the fixed maximum from the phase-space hierarchy.
 * `active_node_count` is the current branch-local count observed in an
 * adaptive ragged quadrature. `parent_path` is ordered by conditioning level
 * and has length `conditioning_level`.
 */
struct CqmomActiveBranchCount
{
    std::size_t conditioning_level = 0;
    ConditioningPath parent_path;
    std::size_t phase_space_dimension = 0;
    std::size_t configured_node_count = 0;
    std::size_t active_node_count = 0;
};

/**
 * Layout of the raw moments that fixed-node recursive CQMOM reconstruction
 * requires.
 *
 * The hierarchy owns phase-space dimensions, node counts, and conditioning
 * order. `required_indices` stores phase-space-order raw moment multi-indices
 * in deterministic recursive-requirement order. The layout stores no moment
 * values and is not an arbitrary sparse tensor selection.
 */
struct CqmomMomentLayout
{
    PhaseSpaceHierarchy hierarchy;
    MultiIndexTable required_indices;
};

/**
 * Values of the raw moments identified by a `CqmomMomentLayout`.
 *
 * `values[i]` stores the raw moment for `layout.required_indices[i]`. Only the
 * moments required by fixed-node recursive CQMOM reconstruction are stored.
 */
struct CqmomMomentData
{
    CqmomMomentData(CqmomMomentLayout layout, std::vector<double> values);

    CqmomMomentLayout layout;
    std::vector<double> values;

    [[nodiscard]] bool contains_required_index(MultiIndexView index) const;
    [[nodiscard]] double moment(MultiIndexView index) const;
};

/**
 * CQMOM case definition: the initial moments and every reconstruction setting,
 * in an open struct that mirrors `QmomCase`.
 *
 * Build or edit a case here, then construct a `CqmomPreparedCase` from it,
 * which is the normal CQMOM workflow. A case can also be passed to
 * `compute_cqmom(...)` or create a state with `initial_state()` directly. Each
 * such call builds the recursive branch structure again, which costs about as
 * much as preparing the case once, but the structure is then not shared with
 * any later call or state.
 *
 * `moments.layout` supplies the hierarchy, conditioning order, and required
 * phase-space-order raw moment indices. `moments.values` stores one value per
 * layout index, so no full rectangular raw moment tensor is needed.
 *
 * `qmom_options.adaptive` decides whether the reconstruction is adaptive: every
 * marginal and conditional inversion and every branch node count adapt
 * together, and the quadrature uses ragged storage. With it false, every
 * branch keeps its configured node count and the storage is rectangular.
 */
struct CqmomCase
{
    CqmomMomentData moments;
    ConditionalMomentOptions conditional_options{};
    QmomOptions qmom_options{};
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::none;
    VandermondeSolver vandermonde_solver = VandermondeSolver::automatic;
    CqmomAbscissaBounds abscissa_bounds{};

    /**
     * Construct the accepted initial phase-space evolution state and
     * canonicalize stored moments from the accepted quadrature.
     *
     * Each call builds the recursive branch structure for the hierarchy and
     * hands it to the new state, which reuses it in every later step. A state
     * created by `CqmomPreparedCase::initial_state()` shares the structure the
     * prepared case already holds instead, so prefer that route whenever more
     * than one state is created on a hierarchy. For the same moments both
     * routes give the same state. Throws
     * `std::runtime_error` when the initial moments cannot be reconstructed and
     * `std::invalid_argument` when `moments.layout` is not the required layout
     * of its hierarchy.
     */
    [[nodiscard]] CqmomState initial_state() const;
};

/**
 * Immutable CQMOM state that phase-space evolution advances.
 *
 * A state holds the CQMOM moments, the reconstruction settings of its case and
 * the accepted terminal quadrature, available through `quadrature()`. It is
 * created by `CqmomPreparedCase::initial_state()`, `CqmomCase::initial_state()`
 * or `state_candidate(...)`, and each evolution step returns a new state.
 */
class CqmomState
{
public:
    CqmomState(const CqmomState& other);
    CqmomState& operator=(const CqmomState& other);
    CqmomState(CqmomState&&) noexcept;
    CqmomState& operator=(CqmomState&&) noexcept;
    ~CqmomState();

    [[nodiscard]] const PhaseSpaceHierarchy& hierarchy() const noexcept;
    [[nodiscard]] const CqmomQuadrature& quadrature() const noexcept;
    [[nodiscard]] const QmomOptions& qmom_options() const noexcept;
    [[nodiscard]] const ConditionalMomentOptions& conditional_options() const noexcept;
    [[nodiscard]] CqmomDiagnosticLevel diagnostic_level() const noexcept;
    [[nodiscard]] VandermondeSolver vandermonde_solver() const noexcept;
    [[nodiscard]] const CqmomAbscissaBounds& abscissa_bounds() const noexcept;
    [[nodiscard]] const std::vector<CqmomActiveBranchCount>& active_branch_counts() const noexcept;
    [[nodiscard]] std::vector<MultiIndex> moment_indices() const;
    [[nodiscard]] const std::vector<double>& moment_values() const noexcept;
    [[nodiscard]] std::optional<double> moment(MultiIndexView index) const;

private:
    struct Data;

    friend class internal::CqmomConstructionAccess;

    explicit CqmomState(std::unique_ptr<Data> data);

    std::unique_ptr<Data> data_;
};

/**
 * A CQMOM evolution state created from moments, or the reason none was.
 *
 * `state` holds the state only on success: take it with `*state` and use it
 * like a state from `CqmomPreparedCase::initial_state()`. When a coordinate
 * realizability check rejected the moments, `rejected_coordinate` names the
 * coordinate and `realizability` holds the verdict on its marginal moments.
 * `reconstruction_status` holds the status of the quadrature reconstruction,
 * when it ran.
 */
struct CqmomStateCandidate
{
    StateCandidateStatus status = StateCandidateStatus::invalid_input;
    std::optional<CqmomState> state;
    std::optional<RealizabilityResult> realizability;
    std::optional<std::size_t> rejected_coordinate;
    std::optional<CqmomStatus> reconstruction_status;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == StateCandidateStatus::success;
    }
};

/**
 * The normal CQMOM case: initial moments and reconstruction settings together
 * with the prepared recursive branch structure.
 *
 * The branch structure follows from the phase-space hierarchy and
 * `qmom_options.adaptive` alone. A prepared case builds it once, and every
 * reconstruction, every state created with `initial_state()`, and every case
 * derived with `config::with_initial_moments(...)` shares it instead of
 * building it again. That reuse is the only difference from `CqmomCase`, and
 * both yield the same quadrature. Load a prepared case with
 * `config::load_cqmom_prepared(...)` or construct it from a `CqmomCase`,
 * inspect the hierarchy and numerical policies through its accessors, then pass
 * it to `compute_cqmom(...)` or create an evolution state with
 * `initial_state()`.
 */
class CqmomPreparedCase
{
public:
    explicit CqmomPreparedCase(CqmomCase cqmom_case);

    CqmomPreparedCase(const CqmomPreparedCase& other);
    CqmomPreparedCase& operator=(const CqmomPreparedCase& other);
    CqmomPreparedCase(CqmomPreparedCase&&) noexcept;
    CqmomPreparedCase& operator=(CqmomPreparedCase&&) noexcept;
    ~CqmomPreparedCase();

    [[nodiscard]] const PhaseSpaceHierarchy& hierarchy() const noexcept;
    [[nodiscard]] const QmomOptions& qmom_options() const noexcept;
    [[nodiscard]] const ConditionalMomentOptions& conditional_options() const noexcept;
    [[nodiscard]] CqmomDiagnosticLevel diagnostic_level() const noexcept;
    [[nodiscard]] VandermondeSolver vandermonde_solver() const noexcept;
    [[nodiscard]] const CqmomAbscissaBounds& abscissa_bounds() const noexcept;

    /**
     * Create an evolution state from other moments, with the hierarchy,
     * settings and branch structure of this case; the moments of the case are
     * not used. This is how a host rebuilds the state of moments it changed
     * itself, for example by transport in physical space. `candidate_moments`
     * must have the required moment layout of the hierarchy: create it once
     * with `create_cqmom_moment_layout(hierarchy())` and replace only the
     * values for each call. The quadrature is reconstructed once, the stored
     * moments are those of the accepted quadrature, and a failure is reported
     * in the result rather than thrown.
     */
    [[nodiscard]] CqmomStateCandidate state_candidate(const CqmomMomentData& candidate_moments) const;

    /**
     * Create an evolution state as above, after checking the marginal moments
     * of every coordinate that has a policy in `coordinate_policies`, unless
     * its `check_mode` is `None`. A check accepts strictly realizable moments,
     * boundary moments when its `accept_boundary` is set, and indeterminate
     * moments when its `accept_indeterminate` is set. Moment data with
     * M[0,...,0] below the cutoff is accepted whatever the verdict, so that
     * moment data without number density is `empty` whatever the policies.
     * Mixed and conditional moments are not checked; the reconstruction tests
     * them. A policy beyond the dimensions of the hierarchy, whatever its
     * mode, and an invalid policy, such as a bounded interval without an
     * interval, are invalid input before any coordinate is checked, as in the
     * moment-source update.
     */
    [[nodiscard]] CqmomStateCandidate state_candidate(const CqmomMomentData& candidate_moments,
                                                      const std::vector<std::optional<RealizabilityPolicy>>& coordinate_policies) const;

    /**
     * Construct the accepted initial phase-space evolution state from the
     * moments of this case, throwing `std::runtime_error` when they cannot be
     * reconstructed.
     */
    [[nodiscard]] CqmomState initial_state() const;

private:
    friend CqmomResult compute_cqmom(const CqmomPreparedCase& cqmom_prepared_case);
    struct Data;

    friend class internal::CqmomConstructionAccess;

    explicit CqmomPreparedCase(std::unique_ptr<Data> data);

    std::unique_ptr<Data> data_;
};

struct CqmomResult
{
    CqmomQuadrature quadrature;
    CqmomStatus status = CqmomStatus::success;
    std::string message;
    std::vector<CqmomDiagnostics> diagnostics;
    // The phase-space dimension whose conditioning level failed. It identifies
    // a failure during the reconstruction and has no meaning on success or when
    // the input was rejected before the reconstruction started.
    std::size_t failed_phase_space_dimension = 0;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == CqmomStatus::success;
    }
};

/**
 * Construct a multidimensional CQMOM quadrature directly from a CQMOM case
 * definition.
 *
 * Each call checks the moment layout and builds the recursive branch structure
 * again; reconstructing from a `CqmomPreparedCase` does neither. The
 * required-moment layout must match the hierarchy stored in
 * `cqmom_case.moments.layout`; a layout that does not is reported as a failed
 * result. `cqmom_case.qmom_options.adaptive`, on by default, decides whether
 * branch node counts and inversions adapt; set it to false for the fixed
 * hierarchy and its rectangular layout.
 */
[[nodiscard]] CqmomResult compute_cqmom(const CqmomCase& cqmom_case);

/**
 * Construct a multidimensional CQMOM quadrature from a prepared CQMOM case, the
 * normal CQMOM workflow.
 *
 * The prepared branch structure is reused, so repeated reconstructions on one
 * hierarchy do not rebuild it.
 */
[[nodiscard]] CqmomResult compute_cqmom(const CqmomPreparedCase& cqmom_prepared_case);

/**
 * Enumerate the minimal fixed-node raw moment multi-indices required by
 * recursive CQMOM reconstruction for the supplied hierarchy.
 */
[[nodiscard]] CqmomMomentLayout create_cqmom_moment_layout(const PhaseSpaceHierarchy& hierarchy);

/**
 * Extract the moments required by `layout` from a phase-space-order raw moment
 * tensor.
 */
[[nodiscard]] CqmomMomentData extract_cqmom_moments(const RawMomentTensor& moment_tensor, const CqmomMomentLayout& layout);

/**
 * Validate and adapt selected raw moments to required CQMOM moment data.
 *
 * The selected result must be successful, must contain exactly the layout's
 * required indices in layout order, and must provide one finite value per
 * required index.
 */
[[nodiscard]] CqmomMomentData create_cqmom_moment_data(const CqmomMomentLayout& layout, const SelectedRawMomentResult& result);

} // namespace qmomkit
