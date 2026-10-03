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
#include <utility>
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
    IndexedMomentSet(std::vector<std::size_t> key_extents, MultiIndexTable indices, std::vector<double> values);

    [[nodiscard]] const MultiIndexTable& indices() const noexcept;
    [[nodiscard]] const std::vector<double>& values() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] bool contains(const MultiIndex& index) const;
    [[nodiscard]] std::size_t offset(MultiIndexView index) const;
    [[nodiscard]] double moment(MultiIndexView index) const;

private:
    std::vector<std::size_t> key_extents_;
    MultiIndexTable indices_;
    std::vector<double> values_;
    // Flat keys for `lookup_offsets_`, in the same order, so neither the sort
    // comparator nor a lookup has to re-derive a key from its multi-index.
    std::vector<std::size_t> lookup_offsets_;
    std::vector<std::size_t> lookup_keys_;
};

/**
 * Non-owning view over CQMOM moments held in an `IndexedMomentSet`.
 *
 * Recursive reconstruction reads its root moments through this view. The
 * caller requests only moments the set holds, and the set must outlive the
 * view.
 */
class CqmomMomentView
{
public:
    explicit CqmomMomentView(const IndexedMomentSet& moments);

    [[nodiscard]] double moment(MultiIndexView index) const;

private:
    const IndexedMomentSet* moments_ = nullptr;
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
    // Position of each entry of `indices`, keyed by its mixed-radix code over
    // `extents` and sorted by that code. Looking an index up by scanning
    // `indices` compares whole vectors, which is why this is materialized once
    // per layout instead: the child-case builder performs one lookup per RHS
    // column per parent node count, and an adaptive layout builds a case for
    // every node count Wheeler could realize.
    std::vector<std::pair<std::size_t, std::size_t>> offsets_by_code;
    // Child solve layouts keyed by the parent node count realized by Wheeler.
    std::vector<BranchMomentChildLayout> child_layouts_by_parent_node_count;
};

/**
 * Immutable fixed-topology reconstruction metadata prepared by a normal
 * CQMOM case. The layout owns no numerical values or solver workspaces.
 */
struct PreparedCqmomMetadata
{
    const CqmomMomentLayout moment_layout;
    const std::shared_ptr<const BranchMomentLayout> root_branch_layout;
    // The branch-layout tree differs between fixed and adaptive reconstruction:
    // an adaptive tree carries a child case for every parent node count Wheeler
    // could realize, a fixed one only the single configured count. Metadata is
    // therefore usable only by a request whose flag matches the one it was
    // built with.
    const bool adaptive_branch_counts = false;

private:
    PreparedCqmomMetadata(CqmomMomentLayout moment_layout, std::shared_ptr<const BranchMomentLayout> root_branch_layout, bool adaptive_branch_counts);

    friend std::shared_ptr<const PreparedCqmomMetadata> create_prepared_cqmom_metadata(CqmomMomentLayout moment_layout, bool adaptive_branch_counts);
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
[[nodiscard]] std::shared_ptr<const PreparedCqmomMetadata> create_prepared_cqmom_metadata(CqmomMomentLayout moment_layout,
                                                                                          bool adaptive_branch_counts);
[[nodiscard]] bool prepared_cqmom_metadata_is_compatible(const PreparedCqmomMetadata& metadata, const PhaseSpaceHierarchy& hierarchy,
                                                         bool adaptive_branch_counts, const CqmomMomentLayout* expected_layout = nullptr);
[[nodiscard]] std::vector<double> extract_wheeler_moments(const BranchMomentSlice& slice);
[[nodiscard]] DenseMatrix make_branch_rhs_matrix(const BranchMomentSlice& slice, std::size_t parent_node_count);
[[nodiscard]] const BranchMomentChildLayout& child_layout_for_parent_node_count(const BranchMomentLayout& layout, std::size_t parent_node_count);

/**
 * Stores multidimensional quadrature data for recursive CQMOM.
 *
 * Fixed results keep rectangular phase-space-order storage. Adaptive
 * branch-specific results keep a ragged tree plus terminal quadrature arrays in
 * deterministic traversal order. This container performs no Wheeler inversion,
 * conditional-moment reconstruction, or phase-space evolution update.
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
    [[nodiscard]] const ConditioningPathTable& terminal_paths() const noexcept;

    [[nodiscard]] double weight(const MultiIndex& node_index) const;
    [[nodiscard]] double abscissa(std::size_t dimension, const MultiIndex& node_index) const;
    [[nodiscard]] std::vector<double> node(const MultiIndex& node_index) const;

    void set_weight(const MultiIndex& node_index, double weight);
    void set_abscissa(std::size_t dimension, const MultiIndex& node_index, double value);
    void set_adaptive_branch_counts_enabled(bool enabled);
    std::size_t add_tree_node(TreeNode node);
    void set_terminal_quadrature(std::vector<double> weights, std::vector<std::vector<double>> abscissas, ConditioningPathTable paths,
                                 std::vector<std::size_t> terminal_node_ids);

private:
    PhaseSpaceHierarchy hierarchy_;
    std::vector<double> weights_;
    std::vector<std::vector<double>> abscissas_;
    std::vector<TreeNode> tree_nodes_;
    std::vector<std::size_t> terminal_node_ids_;
    ConditioningPathTable terminal_paths_;
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

struct RecursiveQuadratureMomentViewRequest
{
    PhaseSpaceHierarchy hierarchy;
    CqmomMomentView moments;
    ConditionalMomentOptions conditional_options;
    WheelerOptions wheeler_options;
    bool adaptive_branch_counts = false;
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::none;
    VandermondeSolver conditional_vandermonde_solver = VandermondeSolver::automatic;
    const CqmomMomentLayout* validated_layout = nullptr;
    std::vector<CqmomActiveBranchCount> active_branch_counts;
    bool allow_active_branch_count_growth = false;
    CqmomCoordinateSupport coordinate_support;
    const PreparedCqmomMetadata* prepared_metadata = nullptr;
    // Normal fixed case/state reconstruction requires prepared metadata.
    // Explicit standalone reconstruction leaves this false and remains
    // dynamically prepared.
    bool require_prepared_metadata = false;
};

struct RecursiveQuadratureMomentDataRequest
{
    PhaseSpaceHierarchy hierarchy;
    CqmomMomentData moments;
    ConditionalMomentOptions conditional_options;
    WheelerOptions wheeler_options;
    bool adaptive_branch_counts = false;
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::none;
    VandermondeSolver conditional_vandermonde_solver = VandermondeSolver::automatic;
    std::vector<CqmomActiveBranchCount> active_branch_counts;
    bool allow_active_branch_count_growth = false;
    CqmomCoordinateSupport coordinate_support;
    const PreparedCqmomMetadata* prepared_metadata = nullptr;
    bool require_prepared_metadata = false;
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
 * Construct recursive CQMOM quadrature from a CQMOM moment view.
 *
 * Traversal follows `hierarchy.conditioning_order()`. Raw moments are supplied
 * in phase-space dimension order. Internally, branch-local moment slices are
 * ordered by conditioning level, so Wheeler consumes the leading axis of that
 * conditioning-order slice, not necessarily phase-space dimension 0.
 * Fixed-mode output uses rectangular phase-space-order storage. Adaptive
 * branch-count output uses terminal quadrature arrays and a ragged tree.
 */
[[nodiscard]] RecursiveQuadratureResult compute_recursive_quadrature_from_moment_view(const RecursiveQuadratureMomentViewRequest& request);

/**
 * Construct recursive CQMOM quadrature from CQMOM moment data.
 *
 * This internal path owns the indexed CqmomMomentData storage for the duration of
 * reconstruction, then delegates to the moment-view recursive algorithm.
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
