#include "internal/recursive_quadrature.hpp"

#include "internal/cqmom_validation.hpp"
#include "internal/linear_algebra.hpp"
#include "internal/numeric_utils.hpp"
#include "internal/performance.hpp"
#include "internal/vandermonde.hpp"
#include "qmomkit/distributions.hpp"
#include "qmomkit/moment_set.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace qmomkit
{
namespace internal
{

IndexedMomentSet::IndexedMomentSet(std::vector<std::size_t> key_extents, std::vector<MultiIndex> indices, std::vector<double> values)
    : key_extents_(std::move(key_extents)), indices_(std::move(indices)), values_(std::move(values))
{
    if (indices_.size() != values_.size())
    {
        throw std::invalid_argument("Indexed CQMOM moment indices and values must have the same size.");
    }

    lookup_offsets_.resize(indices_.size());
    std::iota(lookup_offsets_.begin(), lookup_offsets_.end(), 0);
    std::sort(lookup_offsets_.begin(), lookup_offsets_.end(), [this](const std::size_t lhs, const std::size_t rhs)
              { return flatten_multi_index(indices_[lhs], key_extents_) < flatten_multi_index(indices_[rhs], key_extents_); });

    for (std::size_t offset = 1; offset < lookup_offsets_.size(); ++offset)
    {
        const std::size_t previous_key = flatten_multi_index(indices_[lookup_offsets_[offset - 1]], key_extents_);
        const std::size_t current_key = flatten_multi_index(indices_[lookup_offsets_[offset]], key_extents_);
        if (previous_key == current_key)
        {
            throw std::invalid_argument("Indexed CQMOM moment storage must not contain duplicate indices.");
        }
    }
}

const std::vector<std::size_t>& IndexedMomentSet::key_extents() const noexcept
{
    return key_extents_;
}

const std::vector<MultiIndex>& IndexedMomentSet::indices() const noexcept
{
    return indices_;
}

const std::vector<double>& IndexedMomentSet::values() const noexcept
{
    return values_;
}

std::vector<double>& IndexedMomentSet::mutable_values() noexcept
{
    return values_;
}

std::size_t IndexedMomentSet::size() const noexcept
{
    return values_.size();
}

bool IndexedMomentSet::contains(const MultiIndex& index) const
{
    try
    {
        (void)offset(index);
    }
    catch (const std::invalid_argument&)
    {
        return false;
    }
    catch (const std::out_of_range&)
    {
        return false;
    }
    return true;
}

std::size_t IndexedMomentSet::offset(const MultiIndex& index) const
{
    const std::size_t key = flatten_multi_index(index, key_extents_);
    const auto found =
        std::lower_bound(lookup_offsets_.begin(), lookup_offsets_.end(), key, [this](const std::size_t candidate, const std::size_t value)
                         { return flatten_multi_index(indices_[candidate], key_extents_) < value; });
    if (found == lookup_offsets_.end() || flatten_multi_index(indices_[*found], key_extents_) != key)
    {
        throw std::out_of_range("Indexed CQMOM moment storage does not contain the requested moment.");
    }
    return *found;
}

double IndexedMomentSet::moment(const MultiIndex& index) const
{
    return values_.at(offset(index));
}

CqmomMomentView::CqmomMomentView(const RawMomentTensor& moments) : dense_moments_(&moments), source_kind_(SourceKind::dense) {}

CqmomMomentView::CqmomMomentView(const IndexedMomentSet& moments) : moment_data_(&moments), source_kind_(SourceKind::moment_data) {}

double CqmomMomentView::moment(const MultiIndex& index) const
{
    if (source_kind_ == SourceKind::moment_data)
    {
        return moment_data_->moment(index);
    }
    return dense_moments_->values[flatten_multi_index(index, dense_moments_->extents)];
}

IndexedMomentSet create_indexed_cqmom_moment_set(const CqmomMomentData& moments)
{
    return IndexedMomentSet(moments.layout.hierarchy.moment_extents(), moments.layout.required_indices, moments.values);
}

namespace
{

void append_unique_index(std::vector<MultiIndex>& indices, MultiIndex index)
{
    if (std::find(indices.begin(), indices.end(), index) == indices.end())
    {
        indices.push_back(std::move(index));
    }
}

[[nodiscard]] std::size_t local_branch_offset(const std::vector<MultiIndex>& indices, const MultiIndex& index)
{
    const auto found = std::find(indices.begin(), indices.end(), index);
    if (found == indices.end())
    {
        throw std::out_of_range("Branch-local moment layout does not contain the requested moment.");
    }
    return static_cast<std::size_t>(std::distance(indices.begin(), found));
}

[[nodiscard]] BranchMomentChildLayout create_child_layout_case(const BranchMomentLayout& layout, const std::size_t parent_node_count,
                                                               const bool adaptive_parent_node_counts)
{
    BranchMomentChildLayout child_case;
    child_case.parent_node_count = parent_node_count;

    // Keep only mixed moments whose leading conditioning-axis order can be
    // represented by the current parent nodes; the remaining indices become
    // the child branch coordinates.
    std::vector<MultiIndex> child_indices;
    for (const MultiIndex& parent_index : layout.indices)
    {
        if (!parent_index.empty() && parent_index.front() < parent_node_count)
        {
            append_unique_index(child_indices, MultiIndex(parent_index.begin() + 1, parent_index.end()));
        }
    }

    child_case.rhs_columns.reserve(child_indices.size());
    // Each child index becomes one RHS column. Rows sweep parent moment order
    // so the Vandermonde solve recovers one child slice per parent node.
    for (const MultiIndex& child_index : child_indices)
    {
        BranchMomentRhsColumn column;
        column.child_index = child_index;
        column.parent_offsets_by_order.reserve(parent_node_count);
        MultiIndex parent_index(layout.extents.size(), 0);
        for (std::size_t axis = 0; axis < child_index.size(); ++axis)
        {
            parent_index[axis + 1] = child_index[axis];
        }
        for (std::size_t parent_order = 0; parent_order < parent_node_count; ++parent_order)
        {
            parent_index[0] = parent_order;
            column.parent_offsets_by_order.push_back(local_branch_offset(layout.indices, parent_index));
        }
        child_case.rhs_columns.push_back(std::move(column));
    }

    const std::vector<std::size_t> child_extents(layout.extents.begin() + 1, layout.extents.end());
    child_case.child_layout = create_branch_moment_layout(child_extents, std::move(child_indices), adaptive_parent_node_counts);
    return child_case;
}

} // namespace

std::shared_ptr<const BranchMomentLayout> create_branch_moment_layout(std::vector<std::size_t> extents, std::vector<MultiIndex> indices,
                                                                      const bool adaptive_parent_node_counts)
{
    internal::performance::record_moment_layout_build();
    if (extents.empty())
    {
        throw std::invalid_argument("Branch-local moment layout must have at least one extent.");
    }

    auto layout = std::make_shared<BranchMomentLayout>();
    layout->extents = std::move(extents);
    layout->indices = std::move(indices);

    // A branch-local layout is ordered by conditioning level. Its leading axis
    // supplies the one-dimensional moment sequence for the next Wheeler solve.
    layout->wheeler_offsets.reserve(layout->extents.front());
    MultiIndex wheeler_index(layout->extents.size(), 0);
    for (std::size_t order = 0; order < layout->extents.front(); ++order)
    {
        wheeler_index.front() = order;
        layout->wheeler_offsets.push_back(local_branch_offset(layout->indices, wheeler_index));
    }

    if (layout->extents.size() == 1)
    {
        return layout;
    }

    const std::size_t max_parent_node_count = layout->extents.front() / 2;
    layout->child_layouts_by_parent_node_count.resize(max_parent_node_count + 1);
    if (adaptive_parent_node_counts)
    {
        for (std::size_t parent_node_count = 1; parent_node_count <= max_parent_node_count; ++parent_node_count)
        {
            layout->child_layouts_by_parent_node_count[parent_node_count] =
                create_child_layout_case(*layout, parent_node_count, adaptive_parent_node_counts);
        }
    }
    else
    {
        layout->child_layouts_by_parent_node_count[max_parent_node_count] =
            create_child_layout_case(*layout, max_parent_node_count, adaptive_parent_node_counts);
    }

    return layout;
}

std::vector<double> extract_wheeler_moments(const BranchMomentSlice& slice)
{
    if (!slice.layout)
    {
        throw std::invalid_argument("Branch-local moment slice must have a layout.");
    }
    if (slice.values.size() != slice.layout->indices.size())
    {
        throw std::invalid_argument("Branch-local moment slice value count does not match its layout.");
    }

    std::vector<double> moments(slice.layout->wheeler_offsets.size(), 0.0);
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        moments[order] = slice.values.at(slice.layout->wheeler_offsets[order]);
    }
    return moments;
}

const BranchMomentChildLayout& child_layout_for_parent_node_count(const BranchMomentLayout& layout, const std::size_t parent_node_count)
{
    if (parent_node_count == 0 || parent_node_count >= layout.child_layouts_by_parent_node_count.size())
    {
        throw std::out_of_range("Branch-local moment layout does not support the requested parent node count.");
    }
    const BranchMomentChildLayout& child_case = layout.child_layouts_by_parent_node_count[parent_node_count];
    if (child_case.parent_node_count != parent_node_count)
    {
        throw std::out_of_range("Branch-local moment layout does not support the requested parent node count.");
    }
    return child_case;
}

DenseMatrix make_branch_rhs_matrix(const BranchMomentSlice& slice, const std::size_t parent_node_count)
{
    if (!slice.layout)
    {
        throw std::invalid_argument("Branch-local moment slice must have a layout.");
    }
    if (slice.values.size() != slice.layout->indices.size())
    {
        throw std::invalid_argument("Branch-local moment slice value count does not match its layout.");
    }

    const BranchMomentChildLayout& child_case = child_layout_for_parent_node_count(*slice.layout, parent_node_count);
    DenseMatrix rhs(parent_node_count, child_case.rhs_columns.size());
    for (std::size_t col = 0; col < child_case.rhs_columns.size(); ++col)
    {
        const BranchMomentRhsColumn& rhs_column = child_case.rhs_columns[col];
        for (std::size_t parent_order = 0; parent_order < parent_node_count; ++parent_order)
        {
            rhs(parent_order, col) = slice.values.at(rhs_column.parent_offsets_by_order[parent_order]);
        }
    }
    return rhs;
}

} // namespace internal

namespace
{

void validate_positive_extents(const std::vector<std::size_t>& extents, const char* name)
{
    if (extents.empty())
    {
        throw std::invalid_argument(std::string(name) + " must not be empty.");
    }
    if (std::any_of(extents.begin(), extents.end(), [](const std::size_t extent) { return extent == 0; }))
    {
        throw std::invalid_argument(std::string(name) + " must contain only positive extents.");
    }
}

[[nodiscard]] std::size_t checked_multiply(const std::size_t lhs, const std::size_t rhs)
{
    if (rhs != 0 && lhs > std::numeric_limits<std::size_t>::max() / rhs)
    {
        throw std::invalid_argument("Index extents overflow size_t.");
    }
    return lhs * rhs;
}

[[nodiscard]] std::vector<std::size_t> default_conditioning_order(const std::size_t dimensions)
{
    std::vector<std::size_t> order(dimensions, 0);
    std::iota(order.begin(), order.end(), 0);
    return order;
}

void validate_conditioning_order(const std::vector<std::size_t>& order, const std::size_t dimensions)
{
    if (order.size() != dimensions)
    {
        throw std::invalid_argument("Conditioning order must contain one entry per phase-space dimension.");
    }

    std::vector<bool> seen(dimensions, false);
    for (const std::size_t dimension : order)
    {
        if (dimension >= dimensions)
        {
            throw std::invalid_argument("Conditioning order contains an out-of-range dimension.");
        }
        if (seen[dimension])
        {
            throw std::invalid_argument("Conditioning order must not contain duplicate dimensions.");
        }
        seen[dimension] = true;
    }
}

[[nodiscard]] std::vector<std::size_t> ordered_extents(const std::vector<std::size_t>& extents, const std::vector<std::size_t>& order)
{
    std::vector<std::size_t> result;
    result.reserve(order.size());
    for (const std::size_t dimension : order)
    {
        result.push_back(extents[dimension]);
    }
    return result;
}

[[nodiscard]] MultiIndex from_ordered_index(const MultiIndex& ordered_index, const std::vector<std::size_t>& order)
{
    MultiIndex result(order.size(), 0);
    for (std::size_t i = 0; i < order.size(); ++i)
    {
        result[order[i]] = ordered_index[i];
    }
    return result;
}

void validate_abscissas(const PhaseSpaceHierarchy& hierarchy, const std::vector<std::vector<double>>& abscissas)
{
    if (abscissas.size() != hierarchy.number_of_dimensions())
    {
        throw std::invalid_argument("Recursive quadrature must store abscissas for every phase-space dimension.");
    }

    for (const std::vector<double>& dimension_abscissas : abscissas)
    {
        if (dimension_abscissas.size() != hierarchy.total_nodes())
        {
            throw std::invalid_argument("Every abscissa dimension must contain one value per quadrature node.");
        }
    }
}

using internal::BranchMomentChildLayout;
using internal::BranchMomentSlice;
using internal::DenseMatrix;
using internal::LinearSolveResult;

struct BranchState
{
    ConditioningPath path;
    std::size_t tree_node_id = 0;
    double accumulated_weight = 0.0;
    std::vector<double> abscissas_by_dimension;
    BranchMomentSlice moment_slice;
};

struct ConditionalSliceResult
{
    std::vector<BranchMomentSlice> child_slices;
    ConditionalMomentStatus status = ConditionalMomentStatus::success;
    std::string message;
    ConditionalMomentDiagnostics diagnostics;
    VandermondeSolver vandermonde_requested_solver = VandermondeSolver::automatic;
    VandermondeSolver vandermonde_selected_solver = VandermondeSolver::dense;
    internal::VandermondeAlgorithm vandermonde_algorithm_used = internal::VandermondeAlgorithm::dense;
    bool dense_vandermonde_fallback_occurred = false;
};

[[nodiscard]] bool valid_conditional_options(const ConditionalMomentOptions& options)
{
    return std::isfinite(options.singular_tolerance) && options.singular_tolerance > 0.0 && std::isfinite(options.conditioning_warning_threshold) &&
           options.conditioning_warning_threshold > 0.0;
}

[[nodiscard]] bool valid_diagnostic_level(const CqmomDiagnosticLevel level)
{
    switch (level)
    {
    case CqmomDiagnosticLevel::none:
    case CqmomDiagnosticLevel::basic:
    case CqmomDiagnosticLevel::detailed:
        return true;
    }
    return false;
}

[[nodiscard]] bool records_success_branch_diagnostics(const CqmomDiagnosticLevel level)
{
    return level != CqmomDiagnosticLevel::none;
}

[[nodiscard]] BranchMomentSlice create_ordered_root_slice(const CqmomMomentLayout& layout, const internal::CqmomMomentView& moments,
                                                          const bool adaptive_branch_counts)
{
    BranchMomentSlice slice;
    const PhaseSpaceHierarchy& hierarchy = layout.hierarchy;
    const std::vector<std::size_t> moment_extents = hierarchy.moment_extents();
    std::vector<std::size_t> ordered_extents;
    ordered_extents.reserve(hierarchy.number_of_dimensions());
    std::vector<MultiIndex> ordered_indices;
    std::vector<double> values;
    ordered_indices.reserve(layout.required_indices.size());
    values.reserve(layout.required_indices.size());
    for (const std::size_t dimension : hierarchy.conditioning_order())
    {
        ordered_extents.push_back(moment_extents[dimension]);
    }

    // Branch-local slices are ordered by conditioning level. Raw/global moment
    // storage remains in phase-space dimension order.
    for (const MultiIndex& phase_space_index : layout.required_indices)
    {
        MultiIndex ordered_index(hierarchy.number_of_dimensions(), 0);
        for (std::size_t axis = 0; axis < hierarchy.number_of_dimensions(); ++axis)
        {
            ordered_index[axis] = phase_space_index[hierarchy.conditioning_order()[axis]];
        }
        ordered_indices.push_back(std::move(ordered_index));
        values.push_back(moments.moment(phase_space_index));
    }
    slice.layout = internal::create_branch_moment_layout(std::move(ordered_extents), std::move(ordered_indices), adaptive_branch_counts);
    slice.values = std::move(values);
    return slice;
}

[[nodiscard]] std::vector<double> one_dimensional_moments(const BranchMomentSlice& slice)
{
    return internal::extract_wheeler_moments(slice);
}

[[nodiscard]] ConditionalSliceResult reconstruct_child_slices(const BranchMomentSlice& parent_slice, const std::vector<double>& parent_abscissas,
                                                              const std::vector<double>& parent_weights, const ConditionalMomentOptions& options,
                                                              const CqmomDiagnosticLevel diagnostic_level, const VandermondeSolver vandermonde_solver)
{
    ConditionalSliceResult result;
    if (!parent_slice.layout || parent_slice.layout->extents.size() < 2)
    {
        result.status = ConditionalMomentStatus::invalid_input;
        result.message = "Branch-local slice must contain a target axis.";
        return result;
    }
    if (parent_abscissas.empty() || parent_abscissas.size() != parent_weights.size())
    {
        result.status = ConditionalMomentStatus::invalid_input;
        result.message = "Parent abscissas and weights must be non-empty and same-sized.";
        return result;
    }
    if (parent_slice.layout->extents.front() < parent_abscissas.size())
    {
        result.status = ConditionalMomentStatus::invalid_input;
        result.message = "Branch-local slice does not contain enough parent moment orders.";
        return result;
    }
    if (!internal::all_finite(parent_abscissas) || !internal::all_finite(parent_weights) || !internal::all_finite(parent_slice.values))
    {
        result.status = ConditionalMomentStatus::invalid_input;
        result.message = "Branch-local conditional moment inputs must be finite.";
        return result;
    }

    const std::size_t parent_node_count = parent_abscissas.size();
    const BranchMomentChildLayout* child_case = nullptr;
    try
    {
        child_case = &internal::child_layout_for_parent_node_count(*parent_slice.layout, parent_node_count);
    }
    catch (const std::exception& ex)
    {
        result.status = ConditionalMomentStatus::invalid_input;
        result.message = ex.what();
        return result;
    }
    const std::size_t rhs_cols = child_case->rhs_columns.size();
    DenseMatrix rhs = internal::make_branch_rhs_matrix(parent_slice, parent_node_count);

    // Solve (V^T diag(w)) C = RHS for C: branch-local conditional moments of
    // the remaining axes at each parent node.
    const internal::VandermondeSolveOptions solve_options{vandermonde_solver, options.singular_tolerance, false, 1.0e-10};
    internal::VandermondeSolveResult solve = internal::solve_weighted_vandermonde(parent_abscissas, parent_weights, rhs, solve_options,
                                                                                  "Branch-local weighted Vandermonde system is singular.");
    result.vandermonde_requested_solver = solve.requested_solver;
    result.vandermonde_selected_solver = solve.selected_solver;
    result.vandermonde_algorithm_used = solve.algorithm_used;
    result.dense_vandermonde_fallback_occurred = solve.dense_fallback_occurred;
    if (!solve.success)
    {
        result.status = ConditionalMomentStatus::singular_system;
        result.message = solve.message;
        return result;
    }

    if (diagnostic_level == CqmomDiagnosticLevel::detailed)
    {
        DenseMatrix system = internal::make_weighted_vandermonde(parent_abscissas, parent_weights);
        const double residual = internal::residual_norm(system, solve.solution, rhs);
        const double condition = internal::condition_estimate(system, options.singular_tolerance);
        result.diagnostics.residual_norm = residual;
        result.diagnostics.condition_estimate = condition;
        result.diagnostics.conditioning_warning = condition > options.conditioning_warning_threshold;
    }

    const BranchMomentSlice empty_child_slice{child_case->child_layout, std::vector<double>(rhs_cols, 0.0)};
    result.child_slices.assign(parent_node_count, empty_child_slice);
    for (BranchMomentSlice& slice : result.child_slices)
    {
        slice.values.assign(rhs_cols, 0.0);
    }
    for (std::size_t node = 0; node < parent_node_count; ++node)
    {
        // Solution row 'node' is the complete moment slice for the child
        // branch conditioned on that parent node.
        for (std::size_t col = 0; col < rhs_cols; ++col)
        {
            result.child_slices[node].values[col] = solve.solution(node, col);
        }
    }

    result.status = ConditionalMomentStatus::success;
    result.message = "Branch-local conditional moments reconstructed.";
    return result;
}

[[nodiscard]] internal::RecursiveQuadratureResult create_recursive_result(const PhaseSpaceHierarchy& hierarchy)
{
    return internal::RecursiveQuadratureResult{internal::RecursiveQuadrature(hierarchy), internal::RecursiveQuadratureStatus::success, "", {}, {}, 0};
}

void record_failure(internal::RecursiveQuadratureResult& result, const internal::RecursiveQuadratureStatus status, const ConditioningPath& path,
                    const std::size_t level, const std::string& message)
{
    result.status = status;
    result.message = message;
    result.failed_path = path;
    result.failed_conditioning_level = level;
}

[[nodiscard]] bool validate_recursive_request_preflight(const ConditionalMomentOptions& conditional_options,
                                                        const CqmomDiagnosticLevel diagnostic_level, const WheelerOptions& wheeler_options,
                                                        const bool adaptive_branch_counts, internal::RecursiveQuadratureResult& result)
{
    if (!valid_conditional_options(conditional_options))
    {
        record_failure(result, internal::RecursiveQuadratureStatus::invalid_input, {}, 0, "Conditional moment options are invalid.");
        return false;
    }
    if (!valid_diagnostic_level(diagnostic_level))
    {
        record_failure(result, internal::RecursiveQuadratureStatus::invalid_input, {}, 0, "CQMOM diagnostic level is invalid.");
        return false;
    }
    if (wheeler_options.adaptive && !adaptive_branch_counts)
    {
        record_failure(result, internal::RecursiveQuadratureStatus::invalid_input, {}, 0,
                       "Adaptive branch-specific node counts are not supported by recursive "
                       "fixed-hierarchy quadrature unless adaptive_branch_counts is enabled.");
        return false;
    }
    return true;
}

[[nodiscard]] WheelerOptions resolve_cqmom_wheeler_options(WheelerOptions options, const bool adaptive_branch_counts)
{
    if (adaptive_branch_counts)
    {
        options.adaptive = true;
    }
    return options;
}

struct ActiveBranchKey
{
    std::size_t conditioning_level = 0;
    ConditioningPath parent_path;
    std::size_t phase_space_dimension = 0;
};

[[nodiscard]] bool operator<(const ActiveBranchKey& lhs, const ActiveBranchKey& rhs) noexcept
{
    if (lhs.conditioning_level != rhs.conditioning_level)
    {
        return lhs.conditioning_level < rhs.conditioning_level;
    }
    if (lhs.parent_path != rhs.parent_path)
    {
        return lhs.parent_path < rhs.parent_path;
    }
    return lhs.phase_space_dimension < rhs.phase_space_dimension;
}

struct ActiveBranchTopology
{
    bool constrained = false;
    std::map<ActiveBranchKey, CqmomActiveBranchCount> counts;
};

struct BranchNodeCountLimit
{
    std::size_t ceiling = 0;
    std::optional<std::size_t> active_node_count;
};

[[nodiscard]] bool path_has_prefix(const ConditioningPath& path, const ConditioningPath& prefix)
{
    if (prefix.size() > path.size())
    {
        return false;
    }
    for (std::size_t level = 0; level < prefix.size(); ++level)
    {
        if (path.node_index_at_level(level) != prefix.node_index_at_level(level))
        {
            return false;
        }
    }
    return true;
}

[[nodiscard]] bool branch_limit_invalidated(const ConditioningPath& parent_path, const std::vector<ConditioningPath>& invalidated_prefixes)
{
    for (const ConditioningPath& prefix : invalidated_prefixes)
    {
        if (path_has_prefix(parent_path, prefix))
        {
            return true;
        }
    }
    return false;
}

[[nodiscard]] ActiveBranchTopology validate_active_branch_topology(const PhaseSpaceHierarchy& hierarchy, const bool adaptive_branch_counts,
                                                                   const std::vector<CqmomActiveBranchCount>& active_branch_counts)
{
    ActiveBranchTopology topology;
    if (!adaptive_branch_counts || active_branch_counts.empty())
    {
        return topology;
    }
    topology.constrained = true;

    const std::size_t dimensions = hierarchy.number_of_dimensions();
    for (const CqmomActiveBranchCount& count : active_branch_counts)
    {
        if (count.conditioning_level >= dimensions)
        {
            throw std::invalid_argument("CQMOM active branch count metadata contains a conditioning level outside the hierarchy.");
        }
        if (count.parent_path.size() != count.conditioning_level)
        {
            throw std::invalid_argument("CQMOM active branch count metadata parent path length must match its conditioning level.");
        }
        const std::size_t expected_dimension = hierarchy.conditioning_order()[count.conditioning_level];
        if (count.phase_space_dimension != expected_dimension)
        {
            throw std::invalid_argument("CQMOM active branch count metadata phase-space dimension does not match the conditioning order.");
        }
        const std::size_t configured_count = hierarchy.node_count(expected_dimension);
        if (count.configured_node_count != configured_count)
        {
            throw std::invalid_argument("CQMOM active branch count metadata does not match the configured hierarchy node count.");
        }
        if (count.active_node_count == 0 || count.active_node_count > configured_count)
        {
            throw std::invalid_argument("CQMOM active branch count metadata is outside the configured node-count capacity.");
        }
        const ActiveBranchKey key{count.conditioning_level, count.parent_path, count.phase_space_dimension};
        if (!topology.counts.emplace(key, count).second)
        {
            throw std::invalid_argument("CQMOM active branch count metadata contains duplicate entries for a branch.");
        }
    }

    std::size_t reachable_entries = 0;
    std::vector<ConditioningPath> parents{ConditioningPath{}};
    for (std::size_t level = 0; level < dimensions; ++level)
    {
        const std::size_t phase_space_dimension = hierarchy.conditioning_order()[level];
        std::vector<ConditioningPath> children;
        for (const ConditioningPath& parent : parents)
        {
            const ActiveBranchKey key{level, parent, phase_space_dimension};
            const auto count = topology.counts.find(key);
            if (count == topology.counts.end())
            {
                throw std::invalid_argument("CQMOM active branch count metadata is missing a reachable branch.");
            }
            ++reachable_entries;
            for (std::size_t child = 0; child < count->second.active_node_count; ++child)
            {
                children.push_back(parent.with_child(child));
            }
        }
        parents = std::move(children);
    }
    if (reachable_entries != topology.counts.size())
    {
        throw std::invalid_argument("CQMOM active branch count metadata contains unreachable or stale branches.");
    }
    return topology;
}

[[nodiscard]] BranchNodeCountLimit recursive_branch_node_count_limit(const ActiveBranchTopology& topology, const std::size_t conditioning_level,
                                                                     const ConditioningPath& parent_path, const std::size_t phase_space_dimension,
                                                                     const std::size_t configured_node_count,
                                                                     const std::vector<ConditioningPath>& invalidated_prefixes)
{
    if (!topology.constrained || branch_limit_invalidated(parent_path, invalidated_prefixes))
    {
        return {configured_node_count, std::nullopt};
    }
    const auto active = topology.counts.find(ActiveBranchKey{conditioning_level, parent_path, phase_space_dimension});
    if (active == topology.counts.end())
    {
        throw std::invalid_argument("CQMOM active branch count metadata is missing a reachable branch.");
    }
    (void)configured_node_count;
    return {active->second.active_node_count, active->second.active_node_count};
}

[[nodiscard]] Quadrature1D compute_branch_quadrature(std::vector<double> moments, const std::size_t node_count, const WheelerOptions& wheeler_options)
{
    return compute_quadrature(MomentSet(std::move(moments)), node_count, wheeler_options);
}

[[nodiscard]] bool validate_coordinate_support(const PhaseSpaceHierarchy& hierarchy, const CqmomCoordinateSupport& support, std::string& message)
{
    if (support.empty())
    {
        return true;
    }
    if (support.size() != hierarchy.number_of_dimensions())
    {
        message = "CQMOM coordinate_support must be empty or contain one entry per phase-space dimension.";
        return false;
    }
    for (std::size_t dimension = 0; dimension < support.size(); ++dimension)
    {
        if (!support[dimension].has_value())
        {
            continue;
        }
        const BoundedInterval interval = *support[dimension];
        if (!std::isfinite(interval.lower) || !std::isfinite(interval.upper) || !(interval.lower < interval.upper))
        {
            message = "CQMOM coordinate_support bounded intervals must have finite lower < upper values.";
            return false;
        }
    }
    return true;
}

[[nodiscard]] std::optional<BoundedInterval> support_interval_for_dimension(const CqmomCoordinateSupport& support,
                                                                            const std::size_t phase_space_dimension)
{
    if (support.empty())
    {
        return std::nullopt;
    }
    return support.at(phase_space_dimension);
}

struct BoundedBranchCandidate
{
    std::size_t node_count = 0;
    std::string attempt_type;
    std::vector<double> moments;
    RealizabilityResult realizability;
    Quadrature1D quadrature;
    bool inversion_attempted = false;
    bool nodes_inside_support = false;
    bool accepted = false;
    std::string message;
};

[[nodiscard]] bool nodes_inside_interval(const std::vector<double>& nodes, const BoundedInterval interval, const double tolerance)
{
    return std::all_of(nodes.begin(), nodes.end(), [interval, tolerance](const double node)
                       { return std::isfinite(node) && node >= interval.lower - tolerance && node <= interval.upper + tolerance; });
}

[[nodiscard]] DenseMatrix bounded_hankel_hadamard_matrix(const std::vector<double>& moments, const std::size_t order)
{
    DenseMatrix matrix(order + 1, order + 1);
    for (std::size_t row = 0; row <= order; ++row)
    {
        for (std::size_t col = 0; col <= order; ++col)
        {
            matrix(row, col) = moments[row + col];
        }
    }
    return matrix;
}

[[nodiscard]] DenseMatrix bounded_linear_localizing_matrix(const std::vector<double>& moments, const BoundedInterval interval,
                                                           const std::size_t order, const bool lower_bound)
{
    DenseMatrix matrix(order + 1, order + 1);
    for (std::size_t row = 0; row <= order; ++row)
    {
        for (std::size_t col = 0; col <= order; ++col)
        {
            const std::size_t index = row + col;
            matrix(row, col) =
                lower_bound ? moments[index + 1] - interval.lower * moments[index] : interval.upper * moments[index] - moments[index + 1];
        }
    }
    return matrix;
}

[[nodiscard]] DenseMatrix bounded_quadratic_localizing_matrix(const std::vector<double>& moments, const BoundedInterval interval,
                                                              const std::size_t order)
{
    DenseMatrix matrix(order + 1, order + 1);
    for (std::size_t row = 0; row <= order; ++row)
    {
        for (std::size_t col = 0; col <= order; ++col)
        {
            const std::size_t index = row + col;
            matrix(row, col) =
                (interval.lower + interval.upper) * moments[index + 1] - interval.lower * interval.upper * moments[index] - moments[index + 2];
        }
    }
    return matrix;
}

[[nodiscard]] double determinant(DenseMatrix matrix)
{
    if (matrix.rows() != matrix.cols())
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double value = 1.0;
    for (std::size_t pivot_col = 0; pivot_col < matrix.rows(); ++pivot_col)
    {
        std::size_t pivot_row = pivot_col;
        double pivot_abs = std::abs(matrix(pivot_col, pivot_col));
        for (std::size_t row = pivot_col + 1; row < matrix.rows(); ++row)
        {
            const double candidate = std::abs(matrix(row, pivot_col));
            if (candidate > pivot_abs)
            {
                pivot_abs = candidate;
                pivot_row = row;
            }
        }
        if (pivot_abs == 0.0)
        {
            return 0.0;
        }
        if (pivot_row != pivot_col)
        {
            for (std::size_t col = 0; col < matrix.cols(); ++col)
            {
                std::swap(matrix(pivot_col, col), matrix(pivot_row, col));
            }
            value = -value;
        }
        const double pivot = matrix(pivot_col, pivot_col);
        value *= pivot;
        for (std::size_t row = pivot_col + 1; row < matrix.rows(); ++row)
        {
            const double factor = matrix(row, pivot_col) / pivot;
            for (std::size_t col = pivot_col + 1; col < matrix.cols(); ++col)
            {
                matrix(row, col) -= factor * matrix(pivot_col, col);
            }
        }
    }
    return value;
}

[[nodiscard]] std::vector<CqmomBoundedSupportDeterminantDiagnostics> bounded_support_determinants(const std::vector<double>& moments,
                                                                                                  const BoundedInterval interval)
{
    std::vector<CqmomBoundedSupportDeterminantDiagnostics> diagnostics;
    if (moments.empty())
    {
        return diagnostics;
    }

    const std::size_t checked_order = moments.size() - 1;
    const std::size_t order = checked_order / 2;
    if (checked_order == 0)
    {
        diagnostics.push_back({"Hankel-Hadamard H_0", 0, determinant(bounded_hankel_hadamard_matrix(moments, 0))});
    }
    else if (checked_order % 2 == 1)
    {
        diagnostics.push_back({"lower localizing L_lower_" + std::to_string(order), order,
                               determinant(bounded_linear_localizing_matrix(moments, interval, order, true))});
        diagnostics.push_back({"upper localizing L_upper_" + std::to_string(order), order,
                               determinant(bounded_linear_localizing_matrix(moments, interval, order, false))});
    }
    else
    {
        diagnostics.push_back({"Hankel-Hadamard H_" + std::to_string(order), order, determinant(bounded_hankel_hadamard_matrix(moments, order))});
        if (order > 0)
        {
            diagnostics.push_back({"quadratic interval localizing L_interval_" + std::to_string(order - 1), order - 1,
                                   determinant(bounded_quadratic_localizing_matrix(moments, interval, order - 1))});
        }
    }
    return diagnostics;
}

[[nodiscard]] BoundedBranchCandidate evaluate_bounded_branch_candidate(const std::vector<double>& moment_sequence, const std::size_t node_count,
                                                                       const BoundedInterval interval, const WheelerOptions& wheeler_options,
                                                                       std::string attempt_type)
{
    BoundedBranchCandidate candidate;
    candidate.node_count = node_count;
    candidate.attempt_type = std::move(attempt_type);
    const std::size_t required_moment_count = 2 * node_count;
    if (moment_sequence.size() < required_moment_count)
    {
        candidate.quadrature.status = QuadratureStatus::invalid_input;
        candidate.message = "Bounded-support branch candidate does not have enough conditional moments.";
        return candidate;
    }
    candidate.moments.assign(moment_sequence.begin(), moment_sequence.begin() + static_cast<std::ptrdiff_t>(required_moment_count));

    const RealizabilityOptions realizability_options{};
    candidate.realizability = check_realizability(MomentSet(candidate.moments), interval, realizability_options);
    if (!candidate.realizability.realizable())
    {
        candidate.quadrature.status = QuadratureStatus::not_realizable;
        candidate.message = "Bounded-support realizability check failed: " + candidate.realizability.message;
        return candidate;
    }

    WheelerOptions fixed_candidate_options = wheeler_options;
    fixed_candidate_options.adaptive = false;
    candidate.inversion_attempted = true;
    candidate.quadrature = compute_branch_quadrature(candidate.moments, node_count, fixed_candidate_options);
    if (!candidate.quadrature.ok())
    {
        candidate.message = "Bounded-support Wheeler inversion failed: " + candidate.quadrature.message;
        return candidate;
    }
    if (candidate.quadrature.weights.size() != node_count || candidate.quadrature.abscissas.size() != node_count)
    {
        candidate.quadrature.status = QuadratureStatus::not_realizable;
        candidate.message = "Bounded-support Wheeler inversion returned an unexpected node count.";
        return candidate;
    }

    candidate.nodes_inside_support = nodes_inside_interval(candidate.quadrature.abscissas, interval, realizability_options.tolerance);
    if (!candidate.nodes_inside_support)
    {
        candidate.quadrature.status = QuadratureStatus::not_realizable;
        candidate.message = "Bounded-support Wheeler inversion returned nodes outside the configured interval.";
        return candidate;
    }

    candidate.accepted = true;
    candidate.message = "Bounded-support branch candidate accepted.";
    return candidate;
}

void record_bounded_candidate_diagnostics(internal::RecursiveQuadratureBranchDiagnostics& diagnostics, const BoundedBranchCandidate& candidate,
                                          const BoundedInterval interval, const std::optional<std::size_t> previous_active_node_count,
                                          const std::size_t accepted_fallback_node_count, const bool include_determinants)
{
    diagnostics.bounded_support_checked = true;
    diagnostics.bounded_support_interval_present = true;
    diagnostics.bounded_support_interval = interval;
    diagnostics.bounded_previous_active_node_count = previous_active_node_count;
    diagnostics.bounded_candidate_node_count = candidate.node_count;
    diagnostics.bounded_candidate_attempt_type = candidate.attempt_type;
    diagnostics.bounded_candidate_moments = candidate.moments;
    diagnostics.bounded_realizability_status = candidate.realizability.status;
    diagnostics.bounded_effective_tolerance = candidate.realizability.tolerance_used;
    diagnostics.bounded_failed_matrix = candidate.realizability.failed_matrix;
    diagnostics.bounded_minimum_margin = candidate.realizability.minimum_pivot;
    diagnostics.bounded_accepted_fallback_node_count = accepted_fallback_node_count;
    if (candidate.inversion_attempted)
    {
        diagnostics.bounded_returned_nodes = candidate.quadrature.abscissas;
        diagnostics.bounded_returned_weights = candidate.quadrature.weights;
    }
    if (!candidate.accepted || include_determinants)
    {
        diagnostics.bounded_support_determinants = bounded_support_determinants(candidate.moments, interval);
    }
}

struct BoundedBranchSelection
{
    Quadrature1D quadrature;
    std::size_t selected_node_count_ceiling = 0;
    bool reduced = false;
    bool success = false;
    std::string message;
};

[[nodiscard]] bool branch_growth_trial_allowed(bool allow_active_branch_count_growth, bool adaptive_branch_counts,
                                               const WheelerOptions& wheeler_options, const BranchNodeCountLimit& limit,
                                               std::size_t configured_node_count, const Quadrature1D& current);

[[nodiscard]] BoundedBranchSelection select_bounded_branch_quadrature(const std::vector<double>& moment_sequence,
                                                                      const BranchNodeCountLimit& node_count_limit,
                                                                      const std::size_t configured_node_count, const BoundedInterval interval,
                                                                      const WheelerOptions& wheeler_options, const bool adaptive_branch_counts,
                                                                      const bool allow_active_branch_count_growth, const bool include_determinants,
                                                                      internal::RecursiveQuadratureBranchDiagnostics& diagnostics)
{
    BoundedBranchSelection selection;
    const std::size_t first_count = node_count_limit.active_node_count.value_or(node_count_limit.ceiling);
    std::optional<BoundedBranchCandidate> rejected_candidate;

    if (!adaptive_branch_counts)
    {
        BoundedBranchCandidate fixed = evaluate_bounded_branch_candidate(moment_sequence, configured_node_count, interval, wheeler_options, "fixed");
        if (!fixed.accepted)
        {
            record_bounded_candidate_diagnostics(diagnostics, fixed, interval, node_count_limit.active_node_count, 0, include_determinants);
            selection.quadrature = std::move(fixed.quadrature);
            selection.selected_node_count_ceiling = configured_node_count;
            selection.message = fixed.message;
            return selection;
        }
        record_bounded_candidate_diagnostics(diagnostics, fixed, interval, node_count_limit.active_node_count, configured_node_count,
                                             include_determinants);
        selection.quadrature = std::move(fixed.quadrature);
        selection.selected_node_count_ceiling = configured_node_count;
        selection.success = true;
        return selection;
    }

    for (std::size_t candidate_count = first_count; candidate_count >= 1; --candidate_count)
    {
        const char* attempt_type =
            candidate_count == first_count ? (node_count_limit.active_node_count.has_value() ? "current" : "configured") : "fallback";
        BoundedBranchCandidate candidate =
            evaluate_bounded_branch_candidate(moment_sequence, candidate_count, interval, wheeler_options, attempt_type);
        if (candidate.accepted)
        {
            if (candidate_count < first_count)
            {
                selection.reduced = true;
                record_bounded_candidate_diagnostics(diagnostics, rejected_candidate.value_or(candidate), interval,
                                                     node_count_limit.active_node_count, candidate_count, include_determinants);
            }
            else
            {
                record_bounded_candidate_diagnostics(diagnostics, candidate, interval, node_count_limit.active_node_count, candidate_count,
                                                     include_determinants);
            }
            selection.quadrature = std::move(candidate.quadrature);
            selection.selected_node_count_ceiling = candidate_count;
            selection.success = true;
            break;
        }
        if (!rejected_candidate.has_value())
        {
            rejected_candidate = candidate;
        }
        if (candidate_count == 1)
        {
            record_bounded_candidate_diagnostics(diagnostics, *rejected_candidate, interval, node_count_limit.active_node_count, 0,
                                                 include_determinants);
            selection.quadrature = std::move(candidate.quadrature);
            selection.selected_node_count_ceiling = first_count;
            selection.message = candidate.message;
            return selection;
        }
    }

    if (branch_growth_trial_allowed(allow_active_branch_count_growth, adaptive_branch_counts, wheeler_options, node_count_limit,
                                    configured_node_count, selection.quadrature) &&
        selection.quadrature.weights.size() == *node_count_limit.active_node_count)
    {
        const std::size_t growth_count = *node_count_limit.active_node_count + 1;
        BoundedBranchCandidate growth = evaluate_bounded_branch_candidate(moment_sequence, growth_count, interval, wheeler_options, "growth");
        if (growth.accepted)
        {
            record_bounded_candidate_diagnostics(diagnostics, growth, interval, node_count_limit.active_node_count, growth_count,
                                                 include_determinants);
            selection.quadrature = std::move(growth.quadrature);
            selection.selected_node_count_ceiling = growth_count;
            selection.reduced = false;
        }
        else
        {
            record_bounded_candidate_diagnostics(diagnostics, growth, interval, node_count_limit.active_node_count,
                                                 selection.quadrature.weights.size(), include_determinants);
        }
    }

    return selection;
}

[[nodiscard]] bool branch_growth_trial_allowed(const bool allow_active_branch_count_growth, const bool adaptive_branch_counts,
                                               const WheelerOptions& wheeler_options, const BranchNodeCountLimit& limit,
                                               const std::size_t configured_node_count, const Quadrature1D& current)
{
    return allow_active_branch_count_growth && adaptive_branch_counts && wheeler_options.adaptive && limit.active_node_count.has_value() &&
           current.ok() && current.weights.size() == *limit.active_node_count && *limit.active_node_count < configured_node_count;
}

} // namespace

std::size_t product_extents(const std::vector<std::size_t>& extents)
{
    validate_positive_extents(extents, "Extents");

    std::size_t product = 1;
    for (const std::size_t extent : extents)
    {
        product = checked_multiply(product, extent);
    }
    return product;
}

std::size_t flatten_multi_index(const MultiIndex& multi_index, const std::vector<std::size_t>& extents)
{
    (void)product_extents(extents);
    if (multi_index.size() != extents.size())
    {
        throw std::invalid_argument("Multi-index rank must match extents rank.");
    }

    std::size_t flat_index = 0;
    for (std::size_t axis = 0; axis < extents.size(); ++axis)
    {
        if (multi_index[axis] >= extents[axis])
        {
            throw std::out_of_range("Multi-index component is outside its extent.");
        }
        flat_index = checked_multiply(flat_index, extents[axis]) + multi_index[axis];
    }
    return flat_index;
}

MultiIndex unflatten_multi_index(const std::size_t flat_index, const std::vector<std::size_t>& extents)
{
    const std::size_t total = product_extents(extents);
    if (flat_index >= total)
    {
        throw std::out_of_range("Flat index is outside the indexed extent product.");
    }

    MultiIndex multi_index(extents.size(), 0);
    std::size_t remainder = flat_index;
    for (std::size_t reverse_axis = extents.size(); reverse_axis > 0; --reverse_axis)
    {
        const std::size_t axis = reverse_axis - 1;
        multi_index[axis] = remainder % extents[axis];
        remainder /= extents[axis];
    }
    return multi_index;
}

std::vector<MultiIndex> enumerate_multi_indices(const std::vector<std::size_t>& extents)
{
    const std::size_t total = product_extents(extents);
    std::vector<MultiIndex> indices;
    indices.reserve(total);
    for (std::size_t flat_index = 0; flat_index < total; ++flat_index)
    {
        indices.push_back(unflatten_multi_index(flat_index, extents));
    }
    return indices;
}

PhaseSpaceHierarchy::PhaseSpaceHierarchy(std::vector<std::size_t> node_counts)
    : node_counts_(std::move(node_counts)), conditioning_order_(default_conditioning_order(node_counts_.size())),
      total_nodes_(product_extents(node_counts_))
{
    validate_conditioning_order(conditioning_order_, node_counts_.size());
}

PhaseSpaceHierarchy::PhaseSpaceHierarchy(std::vector<std::size_t> node_counts, std::vector<std::size_t> conditioning_order)
    : node_counts_(std::move(node_counts)), conditioning_order_(std::move(conditioning_order)), total_nodes_(product_extents(node_counts_))
{
    validate_conditioning_order(conditioning_order_, node_counts_.size());
}

std::size_t PhaseSpaceHierarchy::number_of_dimensions() const noexcept
{
    return node_counts_.size();
}

const std::vector<std::size_t>& PhaseSpaceHierarchy::node_counts() const noexcept
{
    return node_counts_;
}

const std::vector<std::size_t>& PhaseSpaceHierarchy::conditioning_order() const noexcept
{
    return conditioning_order_;
}

std::size_t PhaseSpaceHierarchy::node_count(const std::size_t dimension) const
{
    if (dimension >= node_counts_.size())
    {
        throw std::out_of_range("Phase-space dimension is out of range.");
    }
    return node_counts_[dimension];
}

std::size_t PhaseSpaceHierarchy::total_nodes() const noexcept
{
    return total_nodes_;
}

std::vector<std::size_t> PhaseSpaceHierarchy::moment_extents() const
{
    std::vector<std::size_t> extents;
    extents.reserve(node_counts_.size());
    for (const std::size_t node_count : node_counts_)
    {
        extents.push_back(checked_multiply(2, node_count));
    }
    return extents;
}

std::size_t PhaseSpaceHierarchy::flatten_node_index(const MultiIndex& node_index) const
{
    if (node_index.size() != number_of_dimensions())
    {
        throw std::invalid_argument("Node multi-index rank must match the hierarchy dimension.");
    }
    return flatten_multi_index(node_index, node_counts_);
}

MultiIndex PhaseSpaceHierarchy::unflatten_node_index(const std::size_t flat_index) const
{
    return unflatten_multi_index(flat_index, node_counts_);
}

std::vector<MultiIndex> enumerate_moment_multi_indices(const PhaseSpaceHierarchy& hierarchy)
{
    return enumerate_multi_indices(hierarchy.moment_extents());
}

CqmomMomentLayout create_cqmom_moment_layout(const PhaseSpaceHierarchy& hierarchy)
{
    internal::performance::record_moment_layout_build();
    const std::size_t dimensions = hierarchy.number_of_dimensions();
    const std::vector<std::size_t>& order = hierarchy.conditioning_order();

    std::vector<MultiIndex> indices;
    for (std::size_t level = 0; level < dimensions; ++level)
    {
        // Required moments are enumerated in conditioning order for one
        // recursive level at a time: previous axes span their node counts,
        // the active axis spans the 2*N Wheeler moments, and later axes stay at
        // order zero until they become active.
        std::vector<std::size_t> ordered_ranges(dimensions, 1);
        for (std::size_t axis = 0; axis < dimensions; ++axis)
        {
            const std::size_t phase_space_dimension = order[axis];
            if (axis < level)
            {
                ordered_ranges[axis] = hierarchy.node_count(phase_space_dimension);
            }
            else if (axis == level)
            {
                ordered_ranges[axis] = checked_multiply(2, hierarchy.node_count(phase_space_dimension));
            }
        }

        const std::size_t level_index_count = product_extents(ordered_ranges);
        for (std::size_t flat_index = 0; flat_index < level_index_count; ++flat_index)
        {
            const MultiIndex ordered_index = unflatten_multi_index(flat_index, ordered_ranges);
            // Public layouts store phase-space-order indices even though the
            // recursive enumeration is performed in conditioning order.
            const MultiIndex phase_space_index = from_ordered_index(ordered_index, order);
            if (std::find(indices.begin(), indices.end(), phase_space_index) == indices.end())
            {
                indices.push_back(phase_space_index);
            }
        }
    }

    return CqmomMomentLayout{hierarchy, std::move(indices)};
}

CqmomMomentData::CqmomMomentData(CqmomMomentLayout moment_layout, std::vector<double> moment_values)
    : layout(std::move(moment_layout)), values(std::move(moment_values))
{
}

bool CqmomMomentData::contains_required_index(const MultiIndex& index) const
{
    return std::find(layout.required_indices.begin(), layout.required_indices.end(), index) != layout.required_indices.end();
}

double CqmomMomentData::moment(const MultiIndex& index) const
{
    if (values.size() != layout.required_indices.size())
    {
        throw std::invalid_argument("CQMOM moment data values must match the layout size.");
    }

    const auto found = std::find(layout.required_indices.begin(), layout.required_indices.end(), index);
    if (found == layout.required_indices.end())
    {
        throw std::out_of_range("Requested moment is not part of the CQMOM moment layout.");
    }

    return values[static_cast<std::size_t>(std::distance(layout.required_indices.begin(), found))];
}

CqmomMomentData extract_cqmom_moments(const RawMomentTensor& moments, const CqmomMomentLayout& layout)
{
    if (moments.extents.size() != layout.hierarchy.number_of_dimensions())
    {
        throw std::invalid_argument("Dense moment tensor extents must match the CQMOM moment layout dimension.");
    }

    const std::size_t expected_values = product_extents(moments.extents);
    if (moments.values.size() != expected_values)
    {
        throw std::invalid_argument("Dense moment tensor value count does not match its extents.");
    }

    std::vector<double> values;
    values.reserve(layout.required_indices.size());
    for (const MultiIndex& index : layout.required_indices)
    {
        if (!internal::multi_index_within_extents(index, moments.extents))
        {
            throw std::out_of_range("Dense moment tensor does not contain a required CQMOM moment.");
        }
        values.push_back(moments.values[flatten_multi_index(index, moments.extents)]);
    }

    return CqmomMomentData{layout, std::move(values)};
}

CqmomMomentData create_cqmom_moment_data(const CqmomMomentLayout& layout, const SelectedRawMomentResult& result)
{
    if (!result.ok())
    {
        throw std::invalid_argument("Selected raw moment result must be successful before conversion to CQMOM moment data.");
    }
    if (result.indices != layout.required_indices)
    {
        throw std::invalid_argument("Selected raw moment indices must exactly match the CQMOM moment layout.");
    }
    if (result.values.size() != layout.required_indices.size())
    {
        throw std::invalid_argument("Selected raw moment value count must match the CQMOM moment layout size.");
    }
    if (!internal::all_finite(result.values))
    {
        throw std::invalid_argument("Selected raw moment values must be finite before conversion to CQMOM moment data.");
    }

    return CqmomMomentData{layout, result.values};
}

std::vector<MultiIndex> enumerate_quadrature_node_multi_indices(const PhaseSpaceHierarchy& hierarchy)
{
    std::vector<MultiIndex> indices;
    indices.reserve(hierarchy.total_nodes());
    const std::vector<std::size_t> traversal_extents = ordered_extents(hierarchy.node_counts(), hierarchy.conditioning_order());
    for (std::size_t flat_index = 0; flat_index < hierarchy.total_nodes(); ++flat_index)
    {
        indices.push_back(from_ordered_index(unflatten_multi_index(flat_index, traversal_extents), hierarchy.conditioning_order()));
    }
    return indices;
}

namespace
{

[[nodiscard]] ConditioningPath phase_space_node_index_to_path(const PhaseSpaceHierarchy& hierarchy, const MultiIndex& node_index)
{
    std::vector<std::size_t> path(hierarchy.number_of_dimensions(), 0);
    const std::vector<std::size_t>& conditioning_order = hierarchy.conditioning_order();
    for (std::size_t level = 0; level < conditioning_order.size(); ++level)
    {
        path[level] = node_index.at(conditioning_order[level]);
    }
    return ConditioningPath(std::move(path));
}

[[nodiscard]] std::vector<ConditioningPath> rectangular_terminal_paths(const PhaseSpaceHierarchy& hierarchy)
{
    std::vector<ConditioningPath> paths;
    paths.reserve(hierarchy.total_nodes());
    for (const MultiIndex& node_index : enumerate_quadrature_node_multi_indices(hierarchy))
    {
        paths.push_back(phase_space_node_index_to_path(hierarchy, node_index));
    }
    return paths;
}

} // namespace

namespace internal
{

RecursiveQuadrature::RecursiveQuadrature(PhaseSpaceHierarchy hierarchy)
    : hierarchy_(std::move(hierarchy)), weights_(hierarchy_.total_nodes(), 0.0),
      abscissas_(hierarchy_.number_of_dimensions(), std::vector<double>(hierarchy_.total_nodes(), 0.0)),
      terminal_paths_(rectangular_terminal_paths(hierarchy_))
{
}

RecursiveQuadrature::RecursiveQuadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights, std::vector<std::vector<double>> abscissas)
    : hierarchy_(std::move(hierarchy)), weights_(std::move(weights)), abscissas_(std::move(abscissas)),
      terminal_paths_(rectangular_terminal_paths(hierarchy_))
{
    if (weights_.size() != hierarchy_.total_nodes())
    {
        throw std::invalid_argument("Recursive quadrature must contain one weight per quadrature node.");
    }
    validate_abscissas(hierarchy_, abscissas_);
}

const PhaseSpaceHierarchy& RecursiveQuadrature::hierarchy() const noexcept
{
    return hierarchy_;
}

std::size_t RecursiveQuadrature::number_of_dimensions() const noexcept
{
    return hierarchy_.number_of_dimensions();
}

std::size_t RecursiveQuadrature::total_nodes() const noexcept
{
    return weights_.size();
}

const std::vector<double>& RecursiveQuadrature::weights() const noexcept
{
    return weights_;
}

const std::vector<std::vector<double>>& RecursiveQuadrature::abscissas() const noexcept
{
    return abscissas_;
}

bool RecursiveQuadrature::adaptive_branch_counts_enabled() const noexcept
{
    return adaptive_branch_counts_enabled_;
}

bool RecursiveQuadrature::has_rectangular_storage() const noexcept
{
    return rectangular_storage_;
}

const std::vector<RecursiveQuadrature::TreeNode>& RecursiveQuadrature::tree_nodes() const noexcept
{
    return tree_nodes_;
}

const std::vector<std::size_t>& RecursiveQuadrature::terminal_node_ids() const noexcept
{
    return terminal_node_ids_;
}

const std::vector<ConditioningPath>& RecursiveQuadrature::terminal_paths() const noexcept
{
    return terminal_paths_;
}

double RecursiveQuadrature::weight(const MultiIndex& node_index) const
{
    if (!rectangular_storage_)
    {
        throw std::logic_error("Ragged recursive quadrature cannot be addressed by rectangular node multi-index.");
    }
    return weights_.at(hierarchy_.flatten_node_index(node_index));
}

double RecursiveQuadrature::abscissa(const std::size_t dimension, const MultiIndex& node_index) const
{
    if (!rectangular_storage_)
    {
        throw std::logic_error("Ragged recursive quadrature cannot be addressed by rectangular node multi-index.");
    }
    if (dimension >= number_of_dimensions())
    {
        throw std::out_of_range("Phase-space dimension is out of range.");
    }
    return abscissas_[dimension].at(hierarchy_.flatten_node_index(node_index));
}

std::vector<double> RecursiveQuadrature::node(const MultiIndex& node_index) const
{
    std::vector<double> values(number_of_dimensions(), 0.0);
    for (std::size_t dimension = 0; dimension < number_of_dimensions(); ++dimension)
    {
        values[dimension] = abscissa(dimension, node_index);
    }
    return values;
}

void RecursiveQuadrature::set_weight(const MultiIndex& node_index, const double weight)
{
    if (!rectangular_storage_)
    {
        throw std::logic_error("Ragged recursive quadrature cannot be addressed by rectangular node multi-index.");
    }
    weights_.at(hierarchy_.flatten_node_index(node_index)) = weight;
}

void RecursiveQuadrature::set_abscissa(const std::size_t dimension, const MultiIndex& node_index, const double value)
{
    if (!rectangular_storage_)
    {
        throw std::logic_error("Ragged recursive quadrature cannot be addressed by rectangular node multi-index.");
    }
    if (dimension >= number_of_dimensions())
    {
        throw std::out_of_range("Phase-space dimension is out of range.");
    }
    abscissas_[dimension].at(hierarchy_.flatten_node_index(node_index)) = value;
}

void RecursiveQuadrature::set_adaptive_branch_counts_enabled(const bool enabled)
{
    adaptive_branch_counts_enabled_ = enabled;
    if (enabled)
    {
        rectangular_storage_ = false;
        weights_.clear();
        abscissas_.assign(number_of_dimensions(), {});
        terminal_paths_.clear();
        terminal_node_ids_.clear();
    }
}

std::size_t RecursiveQuadrature::add_tree_node(TreeNode node)
{
    node.id = tree_nodes_.size();
    const std::size_t parent_id = node.parent_id;
    tree_nodes_.push_back(std::move(node));
    if (parent_id != std::numeric_limits<std::size_t>::max())
    {
        tree_nodes_.at(parent_id).children.push_back(tree_nodes_.back().id);
    }
    return tree_nodes_.back().id;
}

void RecursiveQuadrature::set_terminal_quadrature(std::vector<double> weights, std::vector<std::vector<double>> abscissas,
                                                  std::vector<ConditioningPath> paths, std::vector<std::size_t> terminal_node_ids)
{
    if (abscissas.size() != number_of_dimensions())
    {
        throw std::invalid_argument("Terminal quadrature must store abscissas for every phase-space dimension.");
    }
    for (const std::vector<double>& dimension_abscissas : abscissas)
    {
        if (dimension_abscissas.size() != weights.size())
        {
            throw std::invalid_argument("Every terminal abscissa dimension must contain one value per terminal node.");
        }
    }
    if (paths.size() != weights.size() || terminal_node_ids.size() != weights.size())
    {
        throw std::invalid_argument("Terminal paths and ids must contain one value per terminal node.");
    }

    weights_ = std::move(weights);
    abscissas_ = std::move(abscissas);
    terminal_paths_ = std::move(paths);
    terminal_node_ids_ = std::move(terminal_node_ids);
    rectangular_storage_ = !adaptive_branch_counts_enabled_;
}

RecursiveQuadratureResult compute_recursive_quadrature_from_moment_view(const RecursiveQuadratureMomentViewRequest& request)
{
    const WheelerOptions wheeler_options = resolve_cqmom_wheeler_options(request.wheeler_options, request.adaptive_branch_counts);
    internal::performance::record_cqmom_reconstruction();
    // Recursive CQMOM traversal:
    // 1. adapt root moments into conditioning-order branch storage;
    // 2. run Wheeler on the leading axis of each branch;
    // 3. reconstruct child conditional moments for the remaining axes;
    // 4. multiply local weights along each path to assemble terminal nodes.
    RecursiveQuadratureResult result = create_recursive_result(request.hierarchy);
    result.quadrature.set_adaptive_branch_counts_enabled(request.adaptive_branch_counts);
    const std::size_t dimensions = request.hierarchy.number_of_dimensions();

    if (!validate_recursive_request_preflight(request.conditional_options, request.diagnostic_level, wheeler_options, request.adaptive_branch_counts,
                                              result))
    {
        return result;
    }
    std::string coordinate_support_message;
    if (!validate_coordinate_support(request.hierarchy, request.coordinate_support, coordinate_support_message))
    {
        record_failure(result, RecursiveQuadratureStatus::invalid_input, {}, 0, coordinate_support_message);
        return result;
    }
    ActiveBranchTopology active_topology;
    try
    {
        active_topology = validate_active_branch_topology(request.hierarchy, request.adaptive_branch_counts, request.active_branch_counts);
    }
    catch (const std::exception& ex)
    {
        record_failure(result, RecursiveQuadratureStatus::invalid_input, {}, 0, ex.what());
        return result;
    }
    BranchMomentSlice root_slice;
    try
    {
        std::optional<CqmomMomentLayout> generated_root_layout;
        const CqmomMomentLayout* root_layout = request.validated_layout;
        if (root_layout == nullptr)
        {
            generated_root_layout.emplace(create_cqmom_moment_layout(request.hierarchy));
            root_layout = &*generated_root_layout;
        }
        root_slice = create_ordered_root_slice(*root_layout, request.moments, request.adaptive_branch_counts);
    }
    catch (const std::exception& ex)
    {
        record_failure(result, RecursiveQuadratureStatus::invalid_input, {}, 0, ex.what());
        return result;
    }

    std::vector<BranchState> current_branches;
    RecursiveQuadrature::TreeNode root_node;
    root_node.id = 0;
    root_node.parent_id = std::numeric_limits<std::size_t>::max();
    root_node.path = ConditioningPath{};
    root_node.abscissas_by_dimension.assign(dimensions, 0.0);
    const std::size_t root_node_id = result.quadrature.add_tree_node(std::move(root_node));

    current_branches.push_back(BranchState{{}, root_node_id, 1.0, std::vector<double>(dimensions, 0.0), std::move(root_slice)});

    std::vector<double> terminal_weights;
    std::vector<std::vector<double>> terminal_abscissas(dimensions);
    std::vector<ConditioningPath> terminal_paths;
    std::vector<std::size_t> terminal_node_ids;
    std::vector<ConditioningPath> invalidated_active_count_prefixes;

    for (std::size_t level = 0; level < dimensions; ++level)
    {
        const std::size_t phase_space_dimension = request.hierarchy.conditioning_order()[level];
        const std::size_t configured_node_count = request.hierarchy.node_count(phase_space_dimension);
        std::vector<BranchState> next_branches;
        if (level + 1 < dimensions)
        {
            next_branches.reserve(current_branches.size() * configured_node_count);
        }

        for (const BranchState& branch : current_branches)
        {
            internal::performance::record_cqmom_branch_solve();
            BranchNodeCountLimit node_count_limit{configured_node_count, std::nullopt};
            try
            {
                node_count_limit = recursive_branch_node_count_limit(active_topology, level, branch.path, phase_space_dimension,
                                                                     configured_node_count, invalidated_active_count_prefixes);
            }
            catch (const std::exception& ex)
            {
                record_failure(result, RecursiveQuadratureStatus::invalid_input, branch.path, level, ex.what());
                return result;
            }
            RecursiveQuadratureBranchDiagnostics branch_diagnostics;
            branch_diagnostics.conditioning_level = level;
            branch_diagnostics.path = branch.path;
            branch_diagnostics.phase_space_dimension = phase_space_dimension;
            branch_diagnostics.requested_node_count = node_count_limit.ceiling;
            branch_diagnostics.adaptive_branch_counts_enabled = request.adaptive_branch_counts;

            std::vector<double> moment_sequence;
            try
            {
                moment_sequence = one_dimensional_moments(branch.moment_slice);
            }
            catch (const std::exception& ex)
            {
                branch_diagnostics.message = ex.what();
                branch_diagnostics.conditional_moment_status = ConditionalMomentStatus::invalid_input;
                result.diagnostics.push_back(branch_diagnostics);
                record_failure(result, RecursiveQuadratureStatus::invalid_input, branch.path, level, ex.what());
                return result;
            }

            Quadrature1D quadrature;
            std::size_t selected_node_count_ceiling = node_count_limit.ceiling;
            bool support_driven_reduction = false;
            const std::optional<BoundedInterval> support_interval = support_interval_for_dimension(request.coordinate_support, phase_space_dimension);
            if (support_interval.has_value())
            {
                BoundedBranchSelection bounded_selection = select_bounded_branch_quadrature(
                    moment_sequence, node_count_limit, configured_node_count, *support_interval, wheeler_options, request.adaptive_branch_counts,
                    request.allow_active_branch_count_growth, request.diagnostic_level == CqmomDiagnosticLevel::detailed, branch_diagnostics);
                quadrature = std::move(bounded_selection.quadrature);
                selected_node_count_ceiling = bounded_selection.selected_node_count_ceiling;
                support_driven_reduction = bounded_selection.reduced;
                if (!bounded_selection.success && !bounded_selection.message.empty())
                {
                    quadrature.message = bounded_selection.message;
                }
            }
            else
            {
                const bool branch_growth_possible = request.allow_active_branch_count_growth && request.adaptive_branch_counts &&
                                                    wheeler_options.adaptive && node_count_limit.active_node_count.has_value() &&
                                                    *node_count_limit.active_node_count < configured_node_count;
                quadrature = branch_growth_possible
                                 ? compute_branch_quadrature(moment_sequence, node_count_limit.ceiling, wheeler_options)
                                 : compute_branch_quadrature(std::move(moment_sequence), node_count_limit.ceiling, wheeler_options);
                if (branch_growth_trial_allowed(request.allow_active_branch_count_growth, request.adaptive_branch_counts, wheeler_options,
                                                node_count_limit, configured_node_count, quadrature))
                {
                    const std::size_t growth_ceiling = *node_count_limit.active_node_count + 1;
                    Quadrature1D growth = compute_branch_quadrature(std::move(moment_sequence), growth_ceiling, wheeler_options);
                    if (growth.ok() && growth.weights.size() == growth_ceiling)
                    {
                        quadrature = std::move(growth);
                        selected_node_count_ceiling = growth_ceiling;
                    }
                }
            }
            branch_diagnostics.requested_node_count = selected_node_count_ceiling;
            branch_diagnostics.wheeler_status = quadrature.status;
            branch_diagnostics.actual_node_count = quadrature.weights.size();
            branch_diagnostics.adaptive_node_count_reduced =
                support_driven_reduction || (quadrature.ok() && quadrature.weights.size() < selected_node_count_ceiling);
            if (!quadrature.ok())
            {
                branch_diagnostics.message = quadrature.message;
                result.diagnostics.push_back(branch_diagnostics);
                record_failure(result, RecursiveQuadratureStatus::wheeler_failed, branch.path, level, quadrature.message);
                return result;
            }
            if (active_topology.constrained && !branch_limit_invalidated(branch.path, invalidated_active_count_prefixes) &&
                node_count_limit.active_node_count.has_value() && quadrature.weights.size() != *node_count_limit.active_node_count)
            {
                for (std::size_t child = 0; child < quadrature.weights.size(); ++child)
                {
                    invalidated_active_count_prefixes.push_back(branch.path.with_child(child));
                }
            }
            if (quadrature.weights.empty() || quadrature.weights.size() != quadrature.abscissas.size() ||
                quadrature.weights.size() > selected_node_count_ceiling ||
                (!request.adaptive_branch_counts && quadrature.weights.size() != selected_node_count_ceiling))
            {
                branch_diagnostics.wheeler_status = quadrature.status;
                branch_diagnostics.message = "Wheeler inversion returned an unexpected node count.";
                result.diagnostics.push_back(branch_diagnostics);
                record_failure(result, RecursiveQuadratureStatus::wheeler_failed, branch.path, level, branch_diagnostics.message);
                return result;
            }

            ConditionalSliceResult conditional_slices;
            if (level + 1 < dimensions)
            {
                conditional_slices =
                    reconstruct_child_slices(branch.moment_slice, quadrature.abscissas, quadrature.weights, request.conditional_options,
                                             request.diagnostic_level, request.conditional_vandermonde_solver);
                branch_diagnostics.conditional_moment_status = conditional_slices.status;
                branch_diagnostics.conditional_moment_diagnostics = conditional_slices.diagnostics;
                branch_diagnostics.conditional_vandermonde_requested_solver = conditional_slices.vandermonde_requested_solver;
                branch_diagnostics.conditional_vandermonde_selected_solver = conditional_slices.vandermonde_selected_solver;
                branch_diagnostics.conditional_vandermonde_algorithm_used = conditional_slices.vandermonde_algorithm_used;
                branch_diagnostics.dense_vandermonde_fallback_occurred = conditional_slices.dense_vandermonde_fallback_occurred;
                if (conditional_slices.status != ConditionalMomentStatus::success)
                {
                    branch_diagnostics.message = conditional_slices.message;
                    result.diagnostics.push_back(branch_diagnostics);
                    record_failure(result, RecursiveQuadratureStatus::conditional_moment_failed, branch.path, level, conditional_slices.message);
                    return result;
                }
            }

            branch_diagnostics.message = "Recursive branch level constructed.";
            if (records_success_branch_diagnostics(request.diagnostic_level))
            {
                result.diagnostics.push_back(branch_diagnostics);
            }

            for (std::size_t node = 0; node < quadrature.weights.size(); ++node)
            {
                BranchState child;
                child.path = branch.path.with_child(node);
                child.tree_node_id = result.quadrature.tree_nodes().size();
                // Terminal weights are products of all conditional local
                // weights encountered along the branch path.
                child.accumulated_weight = branch.accumulated_weight * quadrature.weights[node];
                child.abscissas_by_dimension = branch.abscissas_by_dimension;
                child.abscissas_by_dimension[phase_space_dimension] = quadrature.abscissas[node];

                RecursiveQuadrature::TreeNode tree_node;
                tree_node.parent_id = branch.tree_node_id;
                tree_node.path = child.path;
                tree_node.conditioning_level = level;
                tree_node.phase_space_dimension = phase_space_dimension;
                tree_node.local_node_index = node;
                tree_node.local_weight = quadrature.weights[node];
                tree_node.accumulated_weight = child.accumulated_weight;
                tree_node.abscissas_by_dimension = child.abscissas_by_dimension;
                tree_node.terminal = level + 1 == dimensions;
                child.tree_node_id = result.quadrature.add_tree_node(std::move(tree_node));

                if (level + 1 == dimensions)
                {
                    internal::performance::record_terminal_node_assembled();
                    terminal_weights.push_back(child.accumulated_weight);
                    for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
                    {
                        terminal_abscissas[dimension].push_back(child.abscissas_by_dimension[dimension]);
                    }
                    terminal_paths.push_back(child.path);
                    terminal_node_ids.push_back(child.tree_node_id);

                    if (!request.adaptive_branch_counts)
                    {
                        // Fixed-node hierarchies also populate rectangular
                        // storage indexed directly by the full node multi-index.
                        MultiIndex node_index;
                        try
                        {
                            node_index = recursive_path_to_node_index(request.hierarchy, child.path);
                        }
                        catch (const std::exception& ex)
                        {
                            branch_diagnostics.message = ex.what();
                            if (records_success_branch_diagnostics(request.diagnostic_level) && !result.diagnostics.empty() &&
                                result.diagnostics.back().conditioning_level == branch_diagnostics.conditioning_level &&
                                result.diagnostics.back().path == branch_diagnostics.path)
                            {
                                result.diagnostics.back().message = branch_diagnostics.message;
                            }
                            else
                            {
                                result.diagnostics.push_back(branch_diagnostics);
                            }
                            record_failure(result, RecursiveQuadratureStatus::invalid_input, child.path, level, ex.what());
                            return result;
                        }
                        result.quadrature.set_weight(node_index, child.accumulated_weight);
                        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
                        {
                            result.quadrature.set_abscissa(dimension, node_index, child.abscissas_by_dimension[dimension]);
                        }
                    }
                }
                else
                {
                    child.moment_slice = conditional_slices.child_slices[node];
                    next_branches.push_back(std::move(child));
                }
            }
        }

        current_branches = std::move(next_branches);
    }

    if (request.adaptive_branch_counts)
    {
        try
        {
            result.quadrature.set_terminal_quadrature(std::move(terminal_weights), std::move(terminal_abscissas), std::move(terminal_paths),
                                                      std::move(terminal_node_ids));
        }
        catch (const std::exception& ex)
        {
            record_failure(result, RecursiveQuadratureStatus::invalid_input, {}, dimensions == 0 ? 0 : dimensions - 1, ex.what());
            return result;
        }
    }

    result.status = RecursiveQuadratureStatus::success;
    result.message = "Recursive quadrature constructed.";
    return result;
}

namespace
{

RecursiveQuadratureResult compute_recursive_quadrature_from_validated_moment_view(const RecursiveQuadratureMomentViewRequest& request,
                                                                                  const CqmomMomentLayout& validated_layout)
{
    RecursiveQuadratureMomentViewRequest validated_request = request;
    validated_request.validated_layout = &validated_layout;
    return compute_recursive_quadrature_from_moment_view(validated_request);
}

} // namespace

RecursiveQuadratureResult compute_recursive_quadrature_from_moment_data(const RecursiveQuadratureMomentDataRequest& request)
{
    const WheelerOptions wheeler_options = resolve_cqmom_wheeler_options(request.wheeler_options, request.adaptive_branch_counts);
    RecursiveQuadratureResult result = create_recursive_result(request.hierarchy);
    result.quadrature.set_adaptive_branch_counts_enabled(request.adaptive_branch_counts);

    if (!validate_recursive_request_preflight(request.conditional_options, request.diagnostic_level, wheeler_options, request.adaptive_branch_counts,
                                              result))
    {
        return result;
    }

    IndexedMomentSet indexed_moments;
    try
    {
        validate_cqmom_moment_data_matches_hierarchy(request.moments, request.hierarchy,
                                                     {"CQMOM moment data layout must match the recursive hierarchy.",
                                                      "CQMOM moment data values must match the layout size.",
                                                      "CQMOM moment data values must be finite."});
        indexed_moments = create_indexed_cqmom_moment_set(request.moments);
    }
    catch (const std::exception& ex)
    {
        record_failure(result, RecursiveQuadratureStatus::invalid_input, {}, 0, ex.what());
        return result;
    }

    return compute_recursive_quadrature_from_validated_moment_view({request.hierarchy, CqmomMomentView(indexed_moments), request.conditional_options,
                                                                    wheeler_options, request.adaptive_branch_counts, request.diagnostic_level,
                                                                    request.conditional_vandermonde_solver, nullptr, request.active_branch_counts,
                                                                    request.allow_active_branch_count_growth, request.coordinate_support},
                                                                   request.moments.layout);
}

RecursiveQuadratureResult compute_recursive_quadrature_from_trusted_moment_data(const RecursiveQuadratureMomentDataRequest& request)
{
    const WheelerOptions wheeler_options = resolve_cqmom_wheeler_options(request.wheeler_options, request.adaptive_branch_counts);
    RecursiveQuadratureResult result = create_recursive_result(request.hierarchy);
    result.quadrature.set_adaptive_branch_counts_enabled(request.adaptive_branch_counts);

    if (!validate_recursive_request_preflight(request.conditional_options, request.diagnostic_level, wheeler_options, request.adaptive_branch_counts,
                                              result))
    {
        return result;
    }

    IndexedMomentSet indexed_moments;
    try
    {
        validate_cqmom_moment_values_match_layout(
            request.moments, {"", "CQMOM moment data values must match the layout size.", "CQMOM moment data values must be finite."});
        indexed_moments = create_indexed_cqmom_moment_set(request.moments);
    }
    catch (const std::exception& ex)
    {
        record_failure(result, RecursiveQuadratureStatus::invalid_input, {}, 0, ex.what());
        return result;
    }

    return compute_recursive_quadrature_from_validated_moment_view({request.hierarchy, CqmomMomentView(indexed_moments), request.conditional_options,
                                                                    wheeler_options, request.adaptive_branch_counts, request.diagnostic_level,
                                                                    request.conditional_vandermonde_solver, nullptr, request.active_branch_counts,
                                                                    request.allow_active_branch_count_growth, request.coordinate_support},
                                                                   request.moments.layout);
}

RecursiveQuadratureResult compute_recursive_quadrature(const RecursiveQuadratureRequest& request)
{
    const WheelerOptions wheeler_options = resolve_cqmom_wheeler_options(request.wheeler_options, request.adaptive_branch_counts);
    RecursiveQuadratureResult result = create_recursive_result(request.hierarchy);
    result.quadrature.set_adaptive_branch_counts_enabled(request.adaptive_branch_counts);

    if (!validate_recursive_request_preflight(request.conditional_options, request.diagnostic_level, wheeler_options, request.adaptive_branch_counts,
                                              result))
    {
        return result;
    }
    try
    {
        validate_dense_cqmom_reconstruction_tensor(request.raw_moments, request.hierarchy, "Raw moment extents must match the hierarchy dimension.",
                                                   "Raw moment data size does not match its extents.", "Raw moments must be finite.",
                                                   "Raw moment extents must contain at least 2*N orders per dimension.");
    }
    catch (const std::invalid_argument& ex)
    {
        record_failure(result, RecursiveQuadratureStatus::invalid_input, {}, 0, ex.what());
        return result;
    }

    return compute_recursive_quadrature_from_moment_view({request.hierarchy,
                                                          CqmomMomentView(request.raw_moments),
                                                          request.conditional_options,
                                                          wheeler_options,
                                                          request.adaptive_branch_counts,
                                                          request.diagnostic_level,
                                                          request.conditional_vandermonde_solver,
                                                          nullptr,
                                                          {},
                                                          false,
                                                          request.coordinate_support});
}

} // namespace internal

} // namespace qmomkit
