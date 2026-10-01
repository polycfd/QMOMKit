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
struct CqmomLeanCase;
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
    friend CqmomResult compute_cqmom(const CqmomLeanCase& lean_case);
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

using CqmomCoordinateSupport = std::vector<std::optional<AbscissaBounds>>;

[[nodiscard]] std::string format_cqmom_settings(bool adaptive_branch_counts, VandermondeSolver vandermonde_solver,
                                                CqmomDiagnosticLevel diagnostic_level);

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
 * Advanced explicit-representation layout for the raw moments required by
 * fixed-node recursive CQMOM reconstruction.
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
 * Advanced explicit lean representation values for the raw moments identified
 * by a `CqmomMomentLayout`.
 *
 * `values[i]` stores the raw moment for `layout.required_indices[i]`. This
 * representation stores only moments required by fixed-node recursive CQMOM
 * reconstruction.
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
 * Advanced explicit lean representation for static CQMOM inversion from
 * required CQMOM moments.
 *
 * `moments.layout` supplies the hierarchy, conditioning order, and required
 * phase-space-order raw moment indices. `moments.values` stores one value per
 * layout index. This path avoids requiring a full rectangular raw moment
 * tensor.
 */
struct CqmomLeanCase
{
    CqmomMomentData moments;
    ConditionalMomentOptions conditional_options;
    QmomOptions qmom_options;
    bool adaptive_branch_counts = false;
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::none;
    VandermondeSolver vandermonde_solver = VandermondeSolver::automatic;
    CqmomCoordinateSupport coordinate_support;
};

/**
 * User-facing mutable CQMOM runtime state.
 *
 * A state owns the current CQMOM moments and reconstruction policies. Normal
 * CQMOM phase-space evolution creates this state from `CqmomPreparedCase`; its accepted
 * terminal quadrature is available through `quadrature()`. The state can be
 * passed to `advance_cqmom_phase_space(...)`.
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
    [[nodiscard]] const CqmomQuadrature& quadrature() const;
    [[nodiscard]] const QmomOptions& qmom_options() const noexcept;
    [[nodiscard]] const ConditionalMomentOptions& conditional_options() const noexcept;
    [[nodiscard]] bool adaptive_branch_counts() const noexcept;
    [[nodiscard]] CqmomDiagnosticLevel diagnostic_level() const noexcept;
    [[nodiscard]] VandermondeSolver vandermonde_solver() const noexcept;
    [[nodiscard]] const CqmomCoordinateSupport& coordinate_support() const noexcept;
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
 * CQMOM case that carries prepared reconstruction structure.
 *
 * A prepared case owns the initial moments and reconstruction options, and in
 * addition the recursive branch structure implied by the phase-space
 * hierarchy. Building that structure once lets every later reconstruction of
 * the same hierarchy skip it; that reuse is the only difference from
 * `CqmomLeanCase`, and both yield the same quadrature. Construct it from a
 * `CqmomLeanCase`, inspect the hierarchy and numerical policies through its
 * accessors, then pass it to `compute_cqmom(...)` or create an evolution state
 * with `initial_state()`.
 */
class CqmomPreparedCase
{
public:
    explicit CqmomPreparedCase(CqmomLeanCase lean_case);

    CqmomPreparedCase(const CqmomPreparedCase& other);
    CqmomPreparedCase& operator=(const CqmomPreparedCase& other);
    CqmomPreparedCase(CqmomPreparedCase&&) noexcept;
    CqmomPreparedCase& operator=(CqmomPreparedCase&&) noexcept;
    ~CqmomPreparedCase();

    [[nodiscard]] const PhaseSpaceHierarchy& hierarchy() const noexcept;
    [[nodiscard]] const QmomOptions& qmom_options() const noexcept;
    [[nodiscard]] const ConditionalMomentOptions& conditional_options() const noexcept;
    [[nodiscard]] bool adaptive_branch_counts() const noexcept;
    [[nodiscard]] CqmomDiagnosticLevel diagnostic_level() const noexcept;
    [[nodiscard]] VandermondeSolver vandermonde_solver() const noexcept;
    [[nodiscard]] const CqmomCoordinateSupport& coordinate_support() const noexcept;
    [[nodiscard]] CqmomState initial_state() const;

private:
    friend CqmomResult compute_cqmom(const CqmomPreparedCase& case_);
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
    std::size_t failed_phase_space_dimension = 0;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == CqmomStatus::success;
    }
};

/**
 * Advanced explicit lean-representation static CQMOM inversion from required
 * raw moment data.
 *
 * The required-moment layout must match the hierarchy stored in
 * `lean_case.moments.layout`. Adaptive branch-specific node counts are enabled
 * explicitly with `CqmomLeanCase::adaptive_branch_counts` and automatically
 * enable adaptive Wheeler for all marginal and conditional inversions.
 * Standalone QMOM remains independently configurable.
 */
[[nodiscard]] CqmomResult compute_cqmom(const CqmomLeanCase& lean_case);

/**
 * Construct a multidimensional CQMOM quadrature from the normal user-facing
 * CQMOM case.
 */
[[nodiscard]] CqmomResult compute_cqmom(const CqmomPreparedCase& case_);

/**
 * Enumerate the minimal fixed-node raw moment multi-indices required by
 * recursive CQMOM reconstruction for the supplied hierarchy.
 */
[[nodiscard]] CqmomMomentLayout create_cqmom_moment_layout(const PhaseSpaceHierarchy& hierarchy);

/**
 * Extract the moments required by `layout` from a dense phase-space-order raw
 * moment tensor.
 */
[[nodiscard]] CqmomMomentData extract_cqmom_moments(const RawMomentTensor& moments, const CqmomMomentLayout& layout);

/**
 * Validate and adapt selected raw moments to required CQMOM moment data.
 *
 * The selected result must be successful, must contain exactly the layout's
 * required indices in layout order, and must provide one finite value per
 * required index.
 */
[[nodiscard]] CqmomMomentData create_cqmom_moment_data(const CqmomMomentLayout& layout, const SelectedRawMomentResult& result);

} // namespace qmomkit
