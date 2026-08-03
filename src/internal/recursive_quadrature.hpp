#pragma once

#include "qmomkit/conditioning_path.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/phase_space_hierarchy.hpp"
#include "qmomkit/quadrature.hpp"
#include "qmomkit/wheeler.hpp"
#include "vandermonde.hpp"

#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace qmomkit::internal
{

/**
 * Ordered indexed storage for CQMOM moment data.
 *
 * Indices and values retain deterministic construction order. Lookup uses a
 * sorted flat-key offset table over the supplied key extents; the object is
 * still a CQMOM requirement container, not a general sparse tensor API.
 */
class IndexedMomentSet
{
public:
    IndexedMomentSet() = default;
    IndexedMomentSet(std::vector<std::size_t> key_extents, std::vector<MultiIndex> indices, std::vector<double> values);

    [[nodiscard]] const std::vector<std::size_t>& key_extents() const noexcept;
    [[nodiscard]] const std::vector<MultiIndex>& indices() const noexcept;
    [[nodiscard]] const std::vector<double>& values() const noexcept;
    [[nodiscard]] std::vector<double>& mutable_values() noexcept;
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] bool contains(const MultiIndex& index) const;
    [[nodiscard]] std::size_t offset(const MultiIndex& index) const;
    [[nodiscard]] double moment(const MultiIndex& index) const;

private:
    std::vector<std::size_t> key_extents_;
    std::vector<MultiIndex> indices_;
    std::vector<double> values_;
    std::vector<std::size_t> lookup_offsets_;
};

/**
 * Non-owning view over CQMOM raw moment storage.
 *
 * The caller is responsible for requesting moments that exist in the adapted
 * storage. This keeps reconstruction code independent of the owning dense or
 * moment_data representation. Referenced dense or indexed storage must outlive
 * the view.
 */
class CqmomMomentView
{
public:
    explicit CqmomMomentView(const RawMomentTensor& moments);
    explicit CqmomMomentView(const IndexedMomentSet& moments);

    [[nodiscard]] double moment(const MultiIndex& index) const;

private:
    enum class SourceKind
    {
        dense,
        moment_data
    };

    const RawMomentTensor* dense_moments_ = nullptr;
    const IndexedMomentSet* moment_data_ = nullptr;
    SourceKind source_kind_ = SourceKind::dense;
};

[[nodiscard]] IndexedMomentSet create_indexed_cqmom_moment_set(const CqmomMomentData& moments);

struct BranchMomentRhsColumn
{
    // Child moment multi-index in branch-local conditioning order.
    MultiIndex child_index;
    // Parent-slice offsets for each parent moment order; these form one RHS
    // column of the conditional Vandermonde solve.
    std::vector<std::size_t> parent_offsets_by_order;
};

struct BranchMomentLayout;

struct BranchMomentChildLayout
{
    // Layout case used when the current branch has this many parent nodes.
    std::size_t parent_node_count = 0;
    // RHS columns assembled from the parent slice to recover child moments.
    std::vector<BranchMomentRhsColumn> rhs_columns;
    // Branch-local layout of the remaining conditioning axes after the solve.
    std::shared_ptr<const BranchMomentLayout> child_layout;
};

struct BranchMomentLayout
{
    // Branch-local moment extents and indices are ordered by conditioning level,
    // not by phase-space dimension.
    std::vector<std::size_t> extents;
    std::vector<MultiIndex> indices;
    // Offsets of the leading-axis moment sequence consumed by Wheeler.
    std::vector<std::size_t> wheeler_offsets;
    // Child solve layouts keyed by the parent node count realized by Wheeler.
    std::vector<BranchMomentChildLayout> child_layouts_by_parent_node_count;
};

/**
 * Immutable fixed-topology reconstruction metadata prepared by a normal
 * CQMOM case. The layout owns no numerical values or solver workspaces.
 */
struct PreparedFixedCqmomMetadata
{
    const CqmomMomentLayout moment_layout;
    const std::shared_ptr<const BranchMomentLayout> root_branch_layout;

private:
    PreparedFixedCqmomMetadata(CqmomMomentLayout moment_layout, std::shared_ptr<const BranchMomentLayout> root_branch_layout);

    friend std::shared_ptr<const PreparedFixedCqmomMetadata> create_prepared_fixed_cqmom_metadata(CqmomMomentLayout moment_layout);
};

struct BranchMomentSlice
{
    // The reconstruction entry point retains the root layout owner for the
    // full traversal. Child layouts are reachable from that root, so slices
    // only need a borrowed layout pointer and never alter reference counts.
    const BranchMomentLayout* layout = nullptr;
    std::vector<double> values;
};

[[nodiscard]] std::shared_ptr<const BranchMomentLayout> create_branch_moment_layout(std::vector<std::size_t> extents, std::vector<MultiIndex> indices,
                                                                                    bool adaptive_parent_node_counts = true);
[[nodiscard]] std::shared_ptr<const PreparedFixedCqmomMetadata> create_prepared_fixed_cqmom_metadata(CqmomMomentLayout moment_layout);
[[nodiscard]] bool prepared_fixed_cqmom_metadata_is_compatible(const PreparedFixedCqmomMetadata& metadata, const PhaseSpaceHierarchy& hierarchy,
                                                               const CqmomMomentLayout* expected_layout = nullptr);
[[nodiscard]] std::vector<double> extract_wheeler_moments(const BranchMomentSlice& slice);
[[nodiscard]] DenseMatrix make_branch_rhs_matrix(const BranchMomentSlice& slice, std::size_t parent_node_count);
[[nodiscard]] const BranchMomentChildLayout& child_layout_for_parent_node_count(const BranchMomentLayout& layout, std::size_t parent_node_count);

/**
 * Stores multidimensional quadrature data for recursive CQMOM.
 *
 * Fixed results keep rectangular phase-space-order storage. Adaptive
 * branch-specific results keep a ragged tree plus terminal quadrature arrays in
 * deterministic traversal order. This container performs no Wheeler inversion,
 * conditional-moment reconstruction, or transport update.
 */
class RecursiveQuadrature
{
public:
    explicit RecursiveQuadrature(PhaseSpaceHierarchy hierarchy);
    RecursiveQuadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights, std::vector<std::vector<double>> abscissas);

    [[nodiscard]] const PhaseSpaceHierarchy& hierarchy() const noexcept;
    [[nodiscard]] std::size_t number_of_dimensions() const noexcept;
    [[nodiscard]] std::size_t total_nodes() const noexcept;

    [[nodiscard]] const std::vector<double>& weights() const noexcept;
    [[nodiscard]] const std::vector<std::vector<double>>& abscissas() const noexcept;
    [[nodiscard]] bool adaptive_branch_counts_enabled() const noexcept;
    [[nodiscard]] bool has_rectangular_storage() const noexcept;

    struct TreeNode
    {
        std::size_t id = 0;
        std::size_t parent_id = std::numeric_limits<std::size_t>::max();
        ConditioningPath path;
        std::size_t conditioning_level = 0;
        std::size_t phase_space_dimension = 0;
        std::size_t local_node_index = 0;
        double local_weight = 1.0;
        double accumulated_weight = 1.0;
        std::vector<double> abscissas_by_dimension;
        std::vector<std::size_t> children;
        bool terminal = false;
    };

    [[nodiscard]] const std::vector<TreeNode>& tree_nodes() const noexcept;
    [[nodiscard]] const std::vector<std::size_t>& terminal_node_ids() const noexcept;
    [[nodiscard]] const std::vector<ConditioningPath>& terminal_paths() const noexcept;

    [[nodiscard]] double weight(const MultiIndex& node_index) const;
    [[nodiscard]] double abscissa(std::size_t dimension, const MultiIndex& node_index) const;
    [[nodiscard]] std::vector<double> node(const MultiIndex& node_index) const;

    void set_weight(const MultiIndex& node_index, double weight);
    void set_abscissa(std::size_t dimension, const MultiIndex& node_index, double value);
    void set_adaptive_branch_counts_enabled(bool enabled);
    std::size_t add_tree_node(TreeNode node);
    void set_terminal_quadrature(std::vector<double> weights, std::vector<std::vector<double>> abscissas, std::vector<ConditioningPath> paths,
                                 std::vector<std::size_t> terminal_node_ids);

private:
    PhaseSpaceHierarchy hierarchy_;
    std::vector<double> weights_;
    std::vector<std::vector<double>> abscissas_;
    std::vector<TreeNode> tree_nodes_;
    std::vector<std::size_t> terminal_node_ids_;
    std::vector<ConditioningPath> terminal_paths_;
    bool adaptive_branch_counts_enabled_ = false;
    bool rectangular_storage_ = true;
};

enum class RecursiveQuadratureStatus
{
    success,
    invalid_input,
    wheeler_failed,
    conditional_moment_failed
};

struct RecursiveQuadratureRequest
{
    PhaseSpaceHierarchy hierarchy;
    RawMomentTensor raw_moments;
    ConditionalMomentOptions conditional_options;
    WheelerOptions wheeler_options;
    bool adaptive_branch_counts = false;
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::detailed;
    VandermondeSolver conditional_vandermonde_solver = VandermondeSolver::automatic;
    CqmomCoordinateSupport coordinate_support;
};

struct RecursiveQuadratureMomentViewRequest
{
    PhaseSpaceHierarchy hierarchy;
    CqmomMomentView moments;
    ConditionalMomentOptions conditional_options;
    WheelerOptions wheeler_options;
    bool adaptive_branch_counts = false;
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::detailed;
    VandermondeSolver conditional_vandermonde_solver = VandermondeSolver::automatic;
    const CqmomMomentLayout* validated_layout = nullptr;
    std::vector<CqmomActiveBranchCount> active_branch_counts;
    bool allow_active_branch_count_growth = false;
    CqmomCoordinateSupport coordinate_support;
    const PreparedFixedCqmomMetadata* prepared_fixed_metadata = nullptr;
    // Normal fixed case/state reconstruction requires prepared metadata.
    // Explicit standalone reconstruction leaves this false and remains
    // dynamically prepared.
    bool require_prepared_fixed_metadata = false;
};

struct RecursiveQuadratureMomentDataRequest
{
    PhaseSpaceHierarchy hierarchy;
    CqmomMomentData moments;
    ConditionalMomentOptions conditional_options;
    WheelerOptions wheeler_options;
    bool adaptive_branch_counts = false;
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::detailed;
    VandermondeSolver conditional_vandermonde_solver = VandermondeSolver::automatic;
    std::vector<CqmomActiveBranchCount> active_branch_counts;
    bool allow_active_branch_count_growth = false;
    CqmomCoordinateSupport coordinate_support;
    const PreparedFixedCqmomMetadata* prepared_fixed_metadata = nullptr;
    bool require_prepared_fixed_metadata = false;
};

struct RecursiveQuadratureBranchDiagnostics
{
    std::size_t conditioning_level = 0;
    ConditioningPath path;
    std::size_t phase_space_dimension = 0;
    std::size_t requested_node_count = 0;
    std::size_t actual_node_count = 0;
    bool adaptive_branch_counts_enabled = false;
    bool adaptive_node_count_reduced = false;
    QuadratureStatus wheeler_status = QuadratureStatus::success;
    ConditionalMomentStatus conditional_moment_status = ConditionalMomentStatus::success;
    ConditionalMomentDiagnostics conditional_moment_diagnostics;
    VandermondeSolver conditional_vandermonde_requested_solver = VandermondeSolver::automatic;
    VandermondeSolver conditional_vandermonde_selected_solver = VandermondeSolver::dense;
    VandermondeAlgorithm conditional_vandermonde_algorithm_used = VandermondeAlgorithm::dense;
    bool dense_vandermonde_fallback_occurred = false;
    std::string message;
    bool bounded_support_checked = false;
    bool bounded_support_interval_present = false;
    BoundedInterval bounded_support_interval;
    std::optional<std::size_t> bounded_previous_active_node_count;
    std::size_t bounded_candidate_node_count = 0;
    std::string bounded_candidate_attempt_type;
    std::vector<double> bounded_candidate_moments;
    RealizabilityStatus bounded_realizability_status = RealizabilityStatus::invalid_input;
    double bounded_effective_tolerance = 0.0;
    std::string bounded_failed_matrix;
    double bounded_minimum_margin = 0.0;
    std::size_t bounded_accepted_fallback_node_count = 0;
    std::vector<double> bounded_returned_nodes;
    std::vector<double> bounded_returned_weights;
    std::vector<CqmomBoundedSupportDeterminantDiagnostics> bounded_support_determinants;
};

struct RecursiveQuadratureResult
{
    RecursiveQuadrature quadrature;
    RecursiveQuadratureStatus status = RecursiveQuadratureStatus::success;
    std::string message;
    std::vector<RecursiveQuadratureBranchDiagnostics> diagnostics;
    ConditioningPath failed_path;
    std::size_t failed_conditioning_level = 0;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == RecursiveQuadratureStatus::success;
    }
};

/**
 * Construct a full recursive CQMOM quadrature.
 *
 * Traversal follows `hierarchy.conditioning_order()`. Raw moments are supplied
 * in phase-space dimension order. Internally, branch-local moment slices are
 * ordered by conditioning level, so Wheeler consumes the leading axis of that
 * conditioning-order slice, not necessarily phase-space dimension 0.
 * Fixed-mode output uses rectangular phase-space-order storage. Adaptive
 * branch-count output uses terminal quadrature arrays and a ragged tree.
 */
[[nodiscard]] RecursiveQuadratureResult compute_recursive_quadrature(const RecursiveQuadratureRequest& request);

/**
 * Construct recursive CQMOM quadrature from an abstract CQMOM moment view.
 * The view supplies fixed-node moments from existing CQMOM storage.
 */
[[nodiscard]] RecursiveQuadratureResult compute_recursive_quadrature_from_moment_view(const RecursiveQuadratureMomentViewRequest& request);

/**
 * Construct recursive CQMOM quadrature from CQMOM moment data.
 *
 * This internal path owns the indexed CqmomMomentData storage for the duration of
 * reconstruction, then delegates to the same moment-view recursive algorithm
 * used by the dense path.
 */
[[nodiscard]] RecursiveQuadratureResult compute_recursive_quadrature_from_moment_data(const RecursiveQuadratureMomentDataRequest& request);

/**
 * Construct recursive CQMOM quadrature from QMOMKit-owned CQMOM moment data.
 *
 * The caller guarantees that immutable layout metadata has already been
 * validated through a QMOMKit-owned case/state construction path. This skips
 * repeated layout/hierarchy validation, but still validates changing moment
 * values and numerical options before using the existing recursive algorithm.
 */
[[nodiscard]] RecursiveQuadratureResult compute_recursive_quadrature_from_trusted_moment_data(const RecursiveQuadratureMomentDataRequest& request);

} // namespace qmomkit::internal
