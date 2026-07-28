#pragma once

#include "qmomkit/conditional_moments.hpp"
#include "qmomkit/conditioning_path.hpp"
#include "qmomkit/phase_space_hierarchy.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/vandermonde_solver.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace qmomkit
{
class CqmomCase;
class CqmomQuadrature;
class CqmomState;
struct CqmomDenseCase;
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
    [[nodiscard]] const PhaseSpaceHierarchy& hierarchy() const noexcept;
    [[nodiscard]] std::size_t number_of_dimensions() const noexcept;
    [[nodiscard]] std::size_t total_nodes() const noexcept;

    [[nodiscard]] const std::vector<double>& weights() const noexcept;
    [[nodiscard]] const std::vector<std::vector<double>>& abscissas() const noexcept;
    [[nodiscard]] bool adaptive_branch_counts_enabled() const noexcept;
    [[nodiscard]] bool has_rectangular_storage() const noexcept;
    [[nodiscard]] const std::vector<ConditioningPath>& terminal_paths() const noexcept;

    [[nodiscard]] double weight(const MultiIndex& node_index) const;
    [[nodiscard]] double abscissa(std::size_t dimension, const MultiIndex& node_index) const;
    [[nodiscard]] std::vector<double> node(const MultiIndex& node_index) const;

private:
    friend CqmomResult compute_cqmom(const CqmomDenseCase& dense_case);
    friend CqmomResult compute_cqmom(const CqmomLeanCase& lean_case);
    friend class internal::CqmomConstructionAccess;

    CqmomQuadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights, std::vector<std::vector<double>> abscissas,
                    bool adaptive_branch_counts_enabled, bool rectangular_storage, std::vector<ConditioningPath> terminal_paths);

    PhaseSpaceHierarchy hierarchy_;
    std::vector<double> weights_;
    std::vector<std::vector<double>> abscissas_;
    std::vector<ConditioningPath> terminal_paths_;
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
    bool dense_vandermonde_fallback_occurred = false;
    std::string message;
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
 * Structural description of the raw moments required by fixed-node recursive
 * CQMOM reconstruction.
 *
 * The hierarchy owns phase-space dimensions, node counts, and conditioning
 * order. `required_indices` stores phase-space-order raw moment multi-indices
 * in deterministic recursive-requirement order. The layout stores no moment
 * values and is not an arbitrary sparse tensor selection.
 */
struct CqmomMomentLayout
{
    PhaseSpaceHierarchy hierarchy;
    std::vector<MultiIndex> required_indices;
};

/**
 * Values for the raw moments identified by a `CqmomMomentLayout`.
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

    [[nodiscard]] bool contains_required_index(const MultiIndex& index) const;
    [[nodiscard]] double moment(const MultiIndex& index) const;
};

/**
 * Dense CQMOM case for multidimensional quadrature construction.
 *
 * `moments` is a raw multidimensional moment tensor flattened in phase-space
 * dimension order. `hierarchy` supplies the phase-space node counts and the
 * conditioning order used by the CQMOM construction. Diagnostic level
 * `detailed` preserves the default residual and condition-estimate reporting;
 * lower levels may leave optional expensive diagnostic values unavailable.
 */
struct CqmomDenseCase
{
    PhaseSpaceHierarchy hierarchy;
    RawMomentTensor moments;
    ConditionalMomentOptions conditional_options;
    QmomOptions qmom_options;
    bool adaptive_branch_counts = false;
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::detailed;
    VandermondeSolver vandermonde_solver = VandermondeSolver::automatic;
};

/**
 * Lean CQMOM case for quadrature from required CQMOM moments.
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
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::detailed;
    VandermondeSolver vandermonde_solver = VandermondeSolver::automatic;
};

/**
 * User-facing mutable CQMOM runtime state.
 *
 * A state owns the current CQMOM moments and reconstruction policies while
 * hiding the storage representation used internally. It is the neutral state
 * object passed through CQMOM transport workflows.
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
    [[nodiscard]] const QmomOptions& qmom_options() const noexcept;
    [[nodiscard]] const ConditionalMomentOptions& conditional_options() const noexcept;
    [[nodiscard]] bool adaptive_branch_counts() const noexcept;
    [[nodiscard]] CqmomDiagnosticLevel diagnostic_level() const noexcept;
    [[nodiscard]] VandermondeSolver vandermonde_solver() const noexcept;
    [[nodiscard]] const std::vector<CqmomActiveBranchCount>& active_branch_counts() const noexcept;

private:
    struct Data;

    friend class internal::CqmomConstructionAccess;

    explicit CqmomState(std::unique_ptr<Data> data);

    std::unique_ptr<Data> data_;
};

/**
 * User-facing CQMOM case for reconstruction workflows.
 *
 * A case owns the initialized CQMOM state and reconstruction options while
 * hiding the storage representation used internally. Users can inspect the
 * phase-space hierarchy and numerical policies, then pass the case directly to
 * `compute_cqmom(...)`.
 */
class CqmomCase
{
public:
    CqmomCase(const CqmomCase& other);
    CqmomCase& operator=(const CqmomCase& other);
    CqmomCase(CqmomCase&&) noexcept;
    CqmomCase& operator=(CqmomCase&&) noexcept;
    ~CqmomCase();

    [[nodiscard]] const PhaseSpaceHierarchy& hierarchy() const noexcept;
    [[nodiscard]] const QmomOptions& qmom_options() const noexcept;
    [[nodiscard]] const ConditionalMomentOptions& conditional_options() const noexcept;
    [[nodiscard]] bool adaptive_branch_counts() const noexcept;
    [[nodiscard]] CqmomDiagnosticLevel diagnostic_level() const noexcept;
    [[nodiscard]] VandermondeSolver vandermonde_solver() const noexcept;
    [[nodiscard]] CqmomState initial_state() const;

private:
    friend CqmomResult compute_cqmom(const CqmomCase& case_);
    struct Data;

    friend class internal::CqmomConstructionAccess;

    explicit CqmomCase(std::unique_ptr<Data> data);

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
 * Construct a multidimensional CQMOM quadrature for a fixed phase-space
 * hierarchy and conditioning order.
 *
 * Fixed rectangular hierarchies remain the default. Adaptive branch-specific
 * node counts are enabled explicitly with `CqmomDenseCase::adaptive_branch_counts`
 * and adaptive QMOM options. This does not implement transport or operator
 * splitting.
 */
[[nodiscard]] CqmomResult compute_cqmom(const CqmomDenseCase& dense_case);

/**
 * Construct a multidimensional CQMOM quadrature from required raw moment data.
 *
 * The required-moment layout must match the hierarchy stored in
 * `lean_case.moments.layout`. Adaptive branch-specific node counts are enabled
 * explicitly with `CqmomLeanCase::adaptive_branch_counts` and
 * adaptive QMOM options.
 */
[[nodiscard]] CqmomResult compute_cqmom(const CqmomLeanCase& lean_case);

/**
 * Construct a multidimensional CQMOM quadrature from a user-facing CQMOM case.
 */
[[nodiscard]] CqmomResult compute_cqmom(const CqmomCase& case_);

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
