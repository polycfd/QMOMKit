#include "internal/recursive_quadrature.hpp"

#include "internal/bounded_quadrature.hpp"
#include "internal/checked_size.hpp"
#include "internal/cqmom_validation.hpp"
#include "internal/linear_algebra.hpp"
#include "internal/numeric_utils.hpp"
#include "internal/performance.hpp"
#include "internal/reconstruction_scaling.hpp"
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

namespace
{

[[nodiscard]] std::shared_ptr<const BranchMomentLayout> create_ordered_root_layout(const CqmomMomentLayout& layout, const bool adaptive_branch_counts)
{
    const PhaseSpaceHierarchy& hierarchy = layout.hierarchy;
    const std::vector<std::size_t> moment_extents = hierarchy.moment_extents();
    std::vector<std::size_t> ordered_extents;
    ordered_extents.reserve(hierarchy.number_of_dimensions());
    std::vector<MultiIndex> ordered_indices;
    ordered_indices.reserve(layout.required_indices.size());
    for (const std::size_t dimension : hierarchy.conditioning_order())
    {
        ordered_extents.push_back(moment_extents[dimension]);
    }

    // Branch-local layouts are ordered by conditioning level. Raw/global
    // moment storage remains in phase-space dimension order.
    for (const MultiIndexView phase_space_index : layout.required_indices)
    {
        MultiIndex ordered_index(hierarchy.number_of_dimensions(), 0);
        for (std::size_t axis = 0; axis < hierarchy.number_of_dimensions(); ++axis)
        {
            ordered_index[axis] = phase_space_index[hierarchy.conditioning_order()[axis]];
        }
        ordered_indices.push_back(std::move(ordered_index));
    }
    return create_branch_moment_layout(std::move(ordered_extents), std::move(ordered_indices), adaptive_branch_counts);
}

// Mixed-radix code of a multi-index over `bounds`, or nullopt when the index
// does not fit. Deliberately does not reuse `flatten_multi_index`, which throws
// its own out_of_range for an out-of-range component; callers here need to tell
// "absent" apart from "malformed" and to keep their own error text.
[[nodiscard]] std::optional<std::size_t> multi_index_code(const MultiIndex& index, const std::vector<std::size_t>& bounds)
{
    if (index.size() != bounds.size())
    {
        return std::nullopt;
    }
    std::size_t code = 0;
    for (std::size_t axis = 0; axis < bounds.size(); ++axis)
    {
        if (index[axis] >= bounds[axis])
        {
            return std::nullopt;
        }
        // An index space too large to code is reported as uncodable rather than
        // wrapped: a wrapped code would collide silently and two distinct
        // indices would be treated as one.
        if (bounds[axis] != 0 && code > (std::numeric_limits<std::size_t>::max() - index[axis]) / bounds[axis])
        {
            return std::nullopt;
        }
        code = code * bounds[axis] + index[axis];
    }
    return code;
}

// Insert `code` into a sorted code set, reporting whether it was new. A
// std::unordered_set would allocate a node per element, which measured as a net
// loss here: the sets are small and short-lived and this code is allocation
// bound, so one growing vector beats per-element nodes.
[[nodiscard]] bool insert_unique_code(std::vector<std::size_t>& sorted_codes, const std::size_t code)
{
    const auto position = std::lower_bound(sorted_codes.begin(), sorted_codes.end(), code);
    if (position != sorted_codes.end() && *position == code)
    {
        return false;
    }
    sorted_codes.insert(position, code);
    return true;
}

// Append `index` unless an equal index is already present. `seen_codes` carries
// the codes appended so far, so the duplicate test is an integer lookup rather
// than a scan comparing whole vectors. First-occurrence order is preserved.
void append_unique_index(std::vector<MultiIndex>& indices, std::vector<std::size_t>& seen_codes, const std::vector<std::size_t>& bounds,
                         MultiIndex index)
{
    const std::optional<std::size_t> code = multi_index_code(index, bounds);
    if (!code.has_value())
    {
        // Outside the layout's own extents: fall back to the exhaustive test
        // rather than risk two distinct indices sharing a code.
        if (std::find(indices.begin(), indices.end(), index) == indices.end())
        {
            indices.push_back(std::move(index));
        }
        return;
    }
    if (insert_unique_code(seen_codes, *code))
    {
        indices.push_back(std::move(index));
    }
}

[[nodiscard]] std::size_t local_branch_offset(const BranchMomentLayout& layout, const MultiIndex& index)
{
    // The table is absent only when some index of this layout could not be
    // coded, in which case coding the query would be wasted work and its result
    // unusable. Testing for it first keeps that case on the scan it needs.
    if (!layout.offsets_by_code.empty())
    {
        const std::optional<std::size_t> code = multi_index_code(index, layout.extents);
        if (code.has_value())
        {
            const auto found =
                std::lower_bound(layout.offsets_by_code.begin(), layout.offsets_by_code.end(), *code,
                                 [](const std::pair<std::size_t, std::size_t>& entry, const std::size_t value) { return entry.first < value; });
            if (found != layout.offsets_by_code.end() && found->first == *code)
            {
                return found->second;
            }
            throw std::out_of_range("Branch-local moment layout does not contain the requested moment.");
        }
    }

    const auto found = std::find(layout.indices.begin(), layout.indices.end(), index);
    if (found == layout.indices.end())
    {
        throw std::out_of_range("Branch-local moment layout does not contain the requested moment.");
    }
    return static_cast<std::size_t>(std::distance(layout.indices.begin(), found));
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
    const std::vector<std::size_t> child_bounds(layout.extents.begin() + 1, layout.extents.end());
    std::vector<std::size_t> child_codes;
    for (const MultiIndex& parent_index : layout.indices)
    {
        if (!parent_index.empty() && parent_index.front() < parent_node_count)
        {
            append_unique_index(child_indices, child_codes, child_bounds, MultiIndex(parent_index.begin() + 1, parent_index.end()));
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
            column.parent_offsets_by_order.push_back(local_branch_offset(layout, parent_index));
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

    // Built before the first lookup below, and reused by every later one.
    layout->offsets_by_code.reserve(layout->indices.size());
    for (std::size_t position = 0; position < layout->indices.size(); ++position)
    {
        const std::optional<std::size_t> code = multi_index_code(layout->indices[position], layout->extents);
        if (!code.has_value())
        {
            // An index outside the declared extents cannot be coded uniquely, so
            // abandon the table and let lookups fall back to scanning.
            layout->offsets_by_code.clear();
            break;
        }
        layout->offsets_by_code.emplace_back(*code, position);
    }
    std::sort(layout->offsets_by_code.begin(), layout->offsets_by_code.end());

    // A branch-local layout is ordered by conditioning level. Its leading axis
    // supplies the one-dimensional moment sequence for the next Wheeler solve.
    layout->wheeler_offsets.reserve(layout->extents.front());
    MultiIndex wheeler_index(layout->extents.size(), 0);
    for (std::size_t order = 0; order < layout->extents.front(); ++order)
    {
        wheeler_index.front() = order;
        layout->wheeler_offsets.push_back(local_branch_offset(*layout, wheeler_index));
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

PreparedCqmomMetadata::PreparedCqmomMetadata(CqmomMomentLayout prepared_moment_layout,
                                             std::shared_ptr<const BranchMomentLayout> prepared_root_branch_layout,
                                             const bool prepared_adaptive_branch_counts)
    : moment_layout(std::move(prepared_moment_layout)), root_branch_layout(std::move(prepared_root_branch_layout)),
      adaptive_branch_counts(prepared_adaptive_branch_counts)
{
    if (root_branch_layout == nullptr)
    {
        throw std::invalid_argument("Prepared CQMOM metadata requires a root branch layout.");
    }
}

std::shared_ptr<const PreparedCqmomMetadata> create_prepared_cqmom_metadata(CqmomMomentLayout moment_layout, const bool adaptive_branch_counts)
{
    validate_fixed_cqmom_moment_layout(moment_layout, "Prepared CQMOM metadata requires the canonical CQMOM moment layout.");
    // An adaptive tree is built the same way, only wider: create_branch_moment_layout
    // materializes a child case for every parent node count Wheeler could
    // realize, so the result covers every runtime outcome and depends on nothing
    // but the hierarchy. That is what makes it cacheable rather than rebuilt per
    // reconstruction.
    std::shared_ptr<const BranchMomentLayout> root_branch_layout = create_ordered_root_layout(moment_layout, adaptive_branch_counts);
    std::shared_ptr<const PreparedCqmomMetadata> metadata(
        new PreparedCqmomMetadata(std::move(moment_layout), std::move(root_branch_layout), adaptive_branch_counts));
    internal::performance::record_cqmom_prepared_metadata_build();
    return metadata;
}

bool prepared_cqmom_metadata_is_compatible(const PreparedCqmomMetadata& metadata, const PhaseSpaceHierarchy& hierarchy,
                                           const bool adaptive_branch_counts, const CqmomMomentLayout* const expected_layout)
{
    if (metadata.root_branch_layout == nullptr || metadata.adaptive_branch_counts != adaptive_branch_counts)
    {
        return false;
    }
    const PhaseSpaceHierarchy& metadata_hierarchy = metadata.moment_layout.hierarchy;
    if (metadata_hierarchy.node_counts() != hierarchy.node_counts() || metadata_hierarchy.conditioning_order() != hierarchy.conditioning_order())
    {
        return false;
    }
    return expected_layout == nullptr || metadata.moment_layout.required_indices == expected_layout->required_indices;
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
    const std::optional<std::size_t> product = internal::checked_size_multiply(lhs, rhs);
    if (!product.has_value())
    {
        throw std::invalid_argument("Index extents overflow size_t.");
    }
    return *product;
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
};

[[nodiscard]] bool valid_conditional_options(const ConditionalMomentOptions& options)
{
    return std::isfinite(options.singular_tolerance) && options.singular_tolerance > 0.0 && std::isfinite(options.conditioning_warning_threshold) &&
           options.conditioning_warning_threshold > 0.0 && std::isfinite(options.moore_penrose_rcond) && options.moore_penrose_rcond >= 0.0 &&
           options.moore_penrose_rcond < 1.0;
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

// The root branch layout is built from the canonical CQMOM moment layout, the
// layout both entry points require the moment data to carry, so the root slice
// takes the moment values in their stored order.
[[nodiscard]] BranchMomentSlice create_ordered_root_slice(const CqmomMomentData& moments, const internal::BranchMomentLayout& root_branch_layout)
{
    BranchMomentSlice slice;
    slice.layout = &root_branch_layout;
    slice.values = moments.values;
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
    const internal::VandermondeSolveOptions solve_options{vandermonde_solver, options.singular_tolerance, false, 1.0e-10,
                                                          options.moore_penrose_rcond};
    internal::VandermondeSolveResult solve = internal::solve_conditional_vandermonde_system(parent_abscissas, parent_weights, rhs, solve_options,
                                                                                            "Branch-local weighted Vandermonde system is singular.");
    result.vandermonde_requested_solver = solve.requested_solver;
    result.vandermonde_selected_solver = solve.selected_solver;
    result.vandermonde_algorithm_used = solve.algorithm_used;
    if (!solve.success)
    {
        result.status = ConditionalMomentStatus::singular_system;
        result.message = solve.message;
        return result;
    }

    if (diagnostic_level == CqmomDiagnosticLevel::detailed)
    {
        DenseMatrix system = internal::make_conditional_vandermonde_matrix(parent_abscissas, parent_weights);
        const double residual = internal::residual_norm(system, solve.solution, rhs);
        const double condition = internal::condition_estimate(system, options.singular_tolerance);
        result.diagnostics.residual_norm = residual;
        result.diagnostics.condition_estimate = condition;
        result.diagnostics.conditioning_warning = condition > options.conditioning_warning_threshold;
    }

    const BranchMomentSlice empty_child_slice{child_case->child_layout.get(), std::vector<double>(rhs_cols, 0.0)};
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
    return internal::RecursiveQuadratureResult{
        internal::RecursiveQuadrature(hierarchy), internal::RecursiveQuadratureStatus::success, "", {}, {}, 0, {}};
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
                                                        const CqmomDiagnosticLevel diagnostic_level, internal::RecursiveQuadratureResult& result)
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
    return true;
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
    internal::performance::record_active_count_validation();
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

[[nodiscard]] Quadrature1D compute_branch_quadrature(std::vector<double> moments, const std::size_t node_count, const WheelerOptions& wheeler_options,
                                                     const double coordinate_scale)
{
    return internal::compute_conditional_quadrature(MomentSet(std::move(moments)), node_count, wheeler_options, coordinate_scale);
}

[[nodiscard]] bool validate_coordinate_support(const PhaseSpaceHierarchy& hierarchy, const CqmomAbscissaBounds& support, std::string& message)
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
        const AbscissaBounds interval = *support[dimension];
        if (!std::isfinite(interval.lower) || !std::isfinite(interval.upper) || !(interval.lower < interval.upper))
        {
            message = "CQMOM coordinate_support bounded intervals must have finite lower < upper values.";
            return false;
        }
    }
    return true;
}

[[nodiscard]] std::optional<AbscissaBounds> support_interval_for_dimension(const CqmomAbscissaBounds& support,
                                                                           const std::size_t phase_space_dimension)
{
    if (support.empty())
    {
        return std::nullopt;
    }
    return support.at(phase_space_dimension);
}

void record_bounded_candidate_diagnostics(internal::RecursiveQuadratureBranchDiagnostics& diagnostics,
                                          const internal::BoundedQuadratureCandidate& candidate, const AbscissaBounds interval,
                                          const std::optional<std::size_t> previous_active_node_count, const std::size_t accepted_fallback_node_count)
{
    diagnostics.bounded_support_checked = true;
    diagnostics.bounded_support_interval_present = true;
    diagnostics.bounded_support_interval = interval;
    diagnostics.bounded_previous_active_node_count = previous_active_node_count;
    diagnostics.bounded_candidate_node_count = candidate.node_count;
    diagnostics.bounded_candidate_attempt_type = candidate.attempt_type;
    diagnostics.bounded_accepted_fallback_node_count = accepted_fallback_node_count;
    if (candidate.inversion_attempted)
    {
        diagnostics.bounded_returned_nodes = candidate.quadrature.abscissas;
        diagnostics.bounded_returned_weights = candidate.quadrature.weights;
    }
}

[[nodiscard]] bool branch_growth_trial_allowed(bool allow_active_branch_count_growth, const WheelerOptions& wheeler_options,
                                               const BranchNodeCountLimit& limit, std::size_t configured_node_count, const Quadrature1D& current);

[[nodiscard]] bool branch_growth_trial_allowed(const bool allow_active_branch_count_growth, const WheelerOptions& wheeler_options,
                                               const BranchNodeCountLimit& limit, const std::size_t configured_node_count,
                                               const Quadrature1D& current)
{
    return allow_active_branch_count_growth && wheeler_options.adaptive && limit.active_node_count.has_value() && current.ok() &&
           current.weights.size() == *limit.active_node_count && *limit.active_node_count < configured_node_count;
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

std::size_t flatten_multi_index(const MultiIndexView multi_index, const std::vector<std::size_t>& extents)
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

bool operator==(const MultiIndexView& lhs, const MultiIndexView& rhs) noexcept
{
    return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

bool operator!=(const MultiIndexView& lhs, const MultiIndexView& rhs) noexcept
{
    return !(lhs == rhs);
}

bool operator==(const MultiIndexView& lhs, const MultiIndex& rhs) noexcept
{
    return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

bool operator==(const MultiIndex& lhs, const MultiIndexView& rhs) noexcept
{
    return rhs == lhs;
}

bool operator!=(const MultiIndexView& lhs, const MultiIndex& rhs) noexcept
{
    return !(lhs == rhs);
}

bool operator!=(const MultiIndex& lhs, const MultiIndexView& rhs) noexcept
{
    return !(rhs == lhs);
}

MultiIndexTable::MultiIndexTable(const std::size_t rank, std::vector<std::size_t> values) : values_(std::move(values)), rank_(rank)
{
    if (rank_ == 0)
    {
        if (!values_.empty())
        {
            throw std::invalid_argument("A rank-zero multi-index table must be empty.");
        }
        return;
    }
    if (values_.size() % rank_ != 0)
    {
        throw std::invalid_argument("Multi-index table storage must contain a whole number of equal-rank indices.");
    }
}

MultiIndexTable::MultiIndexTable(const std::vector<MultiIndex>& indices)
{
    if (indices.empty())
    {
        return;
    }
    rank_ = indices.front().size();
    values_.reserve(indices.size() * rank_);
    for (const MultiIndex& index : indices)
    {
        if (index.size() != rank_)
        {
            throw std::invalid_argument("Multi-index table entries must all have the same rank.");
        }
        values_.insert(values_.end(), index.begin(), index.end());
    }
}

MultiIndexView MultiIndexTable::at(const std::size_t row) const
{
    if (row >= size())
    {
        throw std::out_of_range("Multi-index table row is outside the table.");
    }
    return (*this)[row];
}

void MultiIndexTable::push_back(const MultiIndex& index)
{
    push_back(MultiIndexView(index.data(), index.size()));
}

void MultiIndexTable::push_back(const MultiIndexView index)
{
    if (values_.empty() && rank_ == 0)
    {
        rank_ = index.size();
    }
    if (index.size() != rank_)
    {
        throw std::invalid_argument("Multi-index table entries must all have the same rank.");
    }
    values_.insert(values_.end(), index.begin(), index.end());
}

void MultiIndexTable::pop_back()
{
    if (values_.empty())
    {
        throw std::out_of_range("Cannot remove an index from an empty multi-index table.");
    }
    values_.resize(values_.size() - rank_);
}

void MultiIndexTable::reserve(const std::size_t rows)
{
    values_.reserve(rows * rank_);
}

void MultiIndexTable::clear() noexcept
{
    values_.clear();
}

std::vector<MultiIndex> MultiIndexTable::to_multi_indices() const
{
    std::vector<MultiIndex> indices;
    indices.reserve(size());
    for (std::size_t row = 0; row < size(); ++row)
    {
        indices.push_back((*this)[row].to_multi_index());
    }
    return indices;
}

bool operator==(const MultiIndexTable& lhs, const MultiIndexTable& rhs) noexcept
{
    return lhs.rank() == rhs.rank() && lhs.values() == rhs.values();
}

bool operator!=(const MultiIndexTable& lhs, const MultiIndexTable& rhs) noexcept
{
    return !(lhs == rhs);
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

std::size_t PhaseSpaceHierarchy::node_count(const std::size_t dimension) const
{
    if (dimension >= node_counts_.size())
    {
        throw std::out_of_range("Phase-space dimension is out of range.");
    }
    return node_counts_[dimension];
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

std::size_t PhaseSpaceHierarchy::flatten_node_index(const MultiIndexView node_index) const
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

CqmomMomentLayout create_cqmom_moment_layout(const PhaseSpaceHierarchy& hierarchy)
{
    internal::performance::record_moment_layout_build();
    const std::size_t dimensions = hierarchy.number_of_dimensions();
    const std::vector<std::size_t>& order = hierarchy.conditioning_order();

    std::vector<MultiIndex> indices;
    // Every level enumerates a sub-box of the same index space, so the levels
    // overlap and the union has to be deduplicated. The widest any phase-space
    // axis gets is the 2*N Wheeler span it spans while it is the active axis,
    // which bounds the code used to test for duplicates.
    std::vector<std::size_t> index_bounds(dimensions, 0);
    for (std::size_t phase_space_dimension = 0; phase_space_dimension < dimensions; ++phase_space_dimension)
    {
        index_bounds[phase_space_dimension] = checked_multiply(2, hierarchy.node_count(phase_space_dimension));
    }
    std::vector<std::size_t> seen_codes;
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
            const std::optional<std::size_t> code = internal::multi_index_code(phase_space_index, index_bounds);
            if (code.has_value())
            {
                if (internal::insert_unique_code(seen_codes, *code))
                {
                    indices.push_back(phase_space_index);
                }
            }
            else if (std::find(indices.begin(), indices.end(), phase_space_index) == indices.end())
            {
                indices.push_back(phase_space_index);
            }
        }
    }

    return CqmomMomentLayout{hierarchy, MultiIndexTable(std::move(indices))};
}

CqmomMomentData::CqmomMomentData(CqmomMomentLayout moment_layout, std::vector<double> moment_values)
    : layout(std::move(moment_layout)), values(std::move(moment_values))
{
}

bool CqmomMomentData::contains_required_index(const MultiIndexView index) const
{
    return std::find(layout.required_indices.begin(), layout.required_indices.end(), index) != layout.required_indices.end();
}

double CqmomMomentData::moment(const MultiIndexView index) const
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

CqmomMomentData extract_cqmom_moments(const RawMomentTensor& moment_tensor, const CqmomMomentLayout& layout)
{
    if (moment_tensor.extents.size() != layout.hierarchy.number_of_dimensions())
    {
        throw std::invalid_argument("The rank of the full tensor of raw moments must match the CQMOM moment layout dimension count.");
    }

    const std::size_t expected_values = product_extents(moment_tensor.extents);
    if (moment_tensor.values.size() != expected_values)
    {
        throw std::invalid_argument("The value count of the full tensor of raw moments does not match its extents.");
    }

    std::vector<double> values;
    values.reserve(layout.required_indices.size());
    for (const MultiIndexView index : layout.required_indices)
    {
        if (!internal::multi_index_within_extents(index, moment_tensor.extents))
        {
            throw std::out_of_range("The full tensor of raw moments does not contain a required CQMOM moment.");
        }
        values.push_back(moment_tensor.values[flatten_multi_index(index, moment_tensor.extents)]);
    }

    return CqmomMomentData{layout, std::move(values)};
}

CqmomMomentData create_cqmom_moment_data(const CqmomMomentLayout& layout, const SelectedRawMomentResult& result)
{
    if (!result.ok())
    {
        throw std::invalid_argument("Selected raw moment result must be successful before conversion to CQMOM moment data.");
    }
    if (MultiIndexTable(result.indices) != layout.required_indices)
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

[[nodiscard]] ConditioningPathTable rectangular_terminal_paths(const PhaseSpaceHierarchy& hierarchy)
{
    const std::size_t levels = hierarchy.number_of_dimensions();
    std::vector<std::size_t> node_indices;
    node_indices.reserve(hierarchy.total_nodes() * levels);

    // Row-major walk over node indices using one reusable buffer.
    // enumerate_multi_indices would heap-allocate one MultiIndex per terminal
    // node purely to iterate over them.
    const std::vector<std::size_t>& extents = hierarchy.node_counts();
    const std::vector<std::size_t>& conditioning_order = hierarchy.conditioning_order();
    MultiIndex node_index(extents.size(), 0);
    for (std::size_t flat_index = 0; flat_index < hierarchy.total_nodes(); ++flat_index)
    {
        for (std::size_t level = 0; level < levels; ++level)
        {
            node_indices.push_back(node_index[conditioning_order[level]]);
        }
        for (std::size_t reverse_axis = extents.size(); reverse_axis > 0; --reverse_axis)
        {
            if (++node_index[reverse_axis - 1] < extents[reverse_axis - 1])
            {
                break;
            }
            node_index[reverse_axis - 1] = 0;
        }
    }
    return ConditioningPathTable(levels, std::move(node_indices));
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

const ConditioningPathTable& RecursiveQuadrature::terminal_paths() const noexcept
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
                                                  ConditioningPathTable paths, std::vector<std::size_t> terminal_node_ids)
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

namespace
{

// The recursive traversal behind both entry points. They have checked the
// request's options and moment values. The moment data carries the canonical
// CQMOM moment layout of `request.hierarchy`, as does any prepared metadata the
// request holds.
RecursiveQuadratureResult compute_recursive_quadrature_from_validated_moment_data(const RecursiveQuadratureMomentDataRequest& request)
{
    const WheelerOptions& wheeler_options = request.wheeler_options;
    internal::performance::record_cqmom_reconstruction();
    // Recursive CQMOM traversal:
    // 1. adapt root moments into conditioning-order branch storage;
    // 2. run Wheeler on the leading axis of each branch;
    // 3. reconstruct child conditional moments for the remaining axes;
    // 4. multiply local weights along each path to assemble terminal nodes.
    RecursiveQuadratureResult result = create_recursive_result(request.hierarchy);
    result.quadrature.set_adaptive_branch_counts_enabled(wheeler_options.adaptive);
    const std::size_t dimensions = request.hierarchy.number_of_dimensions();

    std::string coordinate_support_message;
    if (!validate_coordinate_support(request.hierarchy, request.coordinate_support, coordinate_support_message))
    {
        record_failure(result, RecursiveQuadratureStatus::invalid_input, {}, 0, coordinate_support_message);
        return result;
    }
    ActiveBranchTopology active_topology;
    try
    {
        active_topology = validate_active_branch_topology(request.hierarchy, wheeler_options.adaptive, request.active_branch_counts);
    }
    catch (const std::exception& ex)
    {
        record_failure(result, RecursiveQuadratureStatus::invalid_input, {}, 0, ex.what());
        return result;
    }
    BranchMomentSlice root_slice;
    std::shared_ptr<const BranchMomentLayout> dynamic_root_branch_layout;
    try
    {
        // Adaptive requests use prepared metadata too: the tree it holds already
        // covers every parent node count, so nothing about it depends on what
        // Wheeler realizes at run time. Compatibility requires the flags to
        // match, so an adaptive request never picks up a fixed tree.
        const bool use_prepared_metadata =
            request.prepared_metadata != nullptr && request.prepared_metadata->adaptive_branch_counts == wheeler_options.adaptive;
        if (use_prepared_metadata)
        {
            const PreparedCqmomMetadata& metadata = *request.prepared_metadata;
            internal::performance::record_cqmom_prepared_metadata_use();
            root_slice = create_ordered_root_slice(request.moments, *metadata.root_branch_layout);
        }
        else
        {
            if (request.require_prepared_metadata)
            {
                throw std::logic_error("Normal CQMOM reconstruction requires prepared metadata.");
            }
            internal::performance::record_cqmom_dynamic_metadata_build();
            dynamic_root_branch_layout = internal::create_ordered_root_layout(request.moments.layout, wheeler_options.adaptive);
            root_slice = create_ordered_root_slice(request.moments, *dynamic_root_branch_layout);
        }
    }
    catch (const std::exception& ex)
    {
        record_failure(result, RecursiveQuadratureStatus::invalid_input, {}, 0, ex.what());
        return result;
    }

    // The root inverts the marginal moments of the first coordinate, which are
    // data. Every later level inverts conditional moments, whose rounding is
    // relative to the scale of their coordinate. Only adaptive reduction uses
    // the scales.
    const std::vector<double> coordinate_scales =
        wheeler_options.adaptive ? cqmom_coordinate_scales(request.moments) : std::vector<double>(dimensions, 0.0);

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
    ConditioningPathTable terminal_paths;
    std::vector<std::size_t> terminal_node_ids;
    std::vector<ConditioningPath> invalidated_active_count_prefixes;

    for (std::size_t level = 0; level < dimensions; ++level)
    {
        const std::size_t phase_space_dimension = request.hierarchy.conditioning_order()[level];
        const std::size_t configured_node_count = request.hierarchy.node_count(phase_space_dimension);
        const double coordinate_scale = level == 0 ? 0.0 : coordinate_scales[phase_space_dimension];
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
            RecursiveQuadratureBranchDiagnostics branch_diagnostics{};
            branch_diagnostics.conditioning_level = level;
            branch_diagnostics.phase_space_dimension = phase_space_dimension;
            branch_diagnostics.requested_node_count = node_count_limit.ceiling;
            branch_diagnostics.adaptive_branch_counts_enabled = wheeler_options.adaptive;
            // Copying the branch path allocates, and at diagnostic level `none`
            // a successful branch publishes no record at all, so the copy would
            // be built and discarded once per branch per reconstruction. The
            // path is instead materialized here, at the single point every
            // record passes through on its way into `result.diagnostics`, which
            // keeps published records identical at every level.
            const auto publish_branch_diagnostics = [&result, &branch](RecursiveQuadratureBranchDiagnostics& record)
            {
                record.path = branch.path;
                result.diagnostics.push_back(record);
            };

            std::vector<double> moment_sequence;
            try
            {
                moment_sequence = one_dimensional_moments(branch.moment_slice);
            }
            catch (const std::exception& ex)
            {
                branch_diagnostics.message = ex.what();
                branch_diagnostics.conditional_moment_status = ConditionalMomentStatus::invalid_input;
                publish_branch_diagnostics(branch_diagnostics);
                record_failure(result, RecursiveQuadratureStatus::invalid_input, branch.path, level, ex.what());
                return result;
            }

            Quadrature1D quadrature;
            std::size_t selected_node_count_ceiling = node_count_limit.ceiling;
            bool support_driven_reduction = false;
            const std::optional<AbscissaBounds> support_interval = support_interval_for_dimension(request.coordinate_support, phase_space_dimension);
            if (support_interval.has_value())
            {
                const internal::BoundedQuadratureSelectionRequest bounded_request{
                    node_count_limit.active_node_count.value_or(node_count_limit.ceiling), configured_node_count, node_count_limit.active_node_count,
                    wheeler_options.adaptive, request.allow_active_branch_count_growth};
                internal::BoundedQuadratureSelection bounded_selection =
                    internal::select_bounded_quadrature(moment_sequence, bounded_request, *support_interval, wheeler_options, coordinate_scale);
                quadrature = std::move(bounded_selection.quadrature);
                selected_node_count_ceiling = bounded_selection.selected_node_count_ceiling;
                support_driven_reduction = bounded_selection.reduced;
                record_bounded_candidate_diagnostics(branch_diagnostics, bounded_selection.diagnostic_candidate, *support_interval,
                                                     node_count_limit.active_node_count, bounded_selection.accepted_fallback_node_count);
            }
            else
            {
                const bool branch_growth_possible = request.allow_active_branch_count_growth && wheeler_options.adaptive &&
                                                    node_count_limit.active_node_count.has_value() &&
                                                    *node_count_limit.active_node_count < configured_node_count;
                quadrature = branch_growth_possible
                                 ? compute_branch_quadrature(moment_sequence, node_count_limit.ceiling, wheeler_options, coordinate_scale)
                                 : compute_branch_quadrature(std::move(moment_sequence), node_count_limit.ceiling, wheeler_options, coordinate_scale);
                if (branch_growth_trial_allowed(request.allow_active_branch_count_growth, wheeler_options, node_count_limit, configured_node_count,
                                                quadrature))
                {
                    const std::size_t growth_ceiling = *node_count_limit.active_node_count + 1;
                    Quadrature1D growth = compute_branch_quadrature(std::move(moment_sequence), growth_ceiling, wheeler_options, coordinate_scale);
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
                publish_branch_diagnostics(branch_diagnostics);
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
                (!wheeler_options.adaptive && quadrature.weights.size() != selected_node_count_ceiling))
            {
                branch_diagnostics.wheeler_status = quadrature.status;
                branch_diagnostics.message = "Wheeler inversion returned an unexpected node count.";
                publish_branch_diagnostics(branch_diagnostics);
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
                if (conditional_slices.status != ConditionalMomentStatus::success)
                {
                    branch_diagnostics.message = conditional_slices.message;
                    publish_branch_diagnostics(branch_diagnostics);
                    record_failure(result, RecursiveQuadratureStatus::conditional_moment_failed, branch.path, level, conditional_slices.message);
                    return result;
                }
            }

            if (records_success_branch_diagnostics(request.diagnostic_level))
            {
                branch_diagnostics.message = "Recursive branch level constructed.";
                publish_branch_diagnostics(branch_diagnostics);
            }
            if (wheeler_options.adaptive && request.record_active_branch_counts)
            {
                result.active_branch_counts.push_back({level, branch.path, phase_space_dimension, configured_node_count, quadrature.weights.size()});
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

                    if (!wheeler_options.adaptive)
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
                                result.diagnostics.back().path == branch.path)
                            {
                                result.diagnostics.back().message = branch_diagnostics.message;
                            }
                            else
                            {
                                publish_branch_diagnostics(branch_diagnostics);
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

    if (wheeler_options.adaptive)
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
        // The traversal visits the branches level by level. Lexicographic
        // parent-path order lists every branch before its descendants, as the
        // terminal paths do.
        std::sort(result.active_branch_counts.begin(), result.active_branch_counts.end(),
                  [](const CqmomActiveBranchCount& lhs, const CqmomActiveBranchCount& rhs) { return lhs.parent_path < rhs.parent_path; });
    }

    result.status = RecursiveQuadratureStatus::success;
    result.message = "Recursive quadrature constructed.";
    return result;
}

} // namespace

RecursiveQuadratureResult compute_recursive_quadrature_from_moment_data(const RecursiveQuadratureMomentDataRequest& request)
{
    RecursiveQuadratureResult result = create_recursive_result(request.hierarchy);
    result.quadrature.set_adaptive_branch_counts_enabled(request.wheeler_options.adaptive);

    if (!validate_recursive_request_preflight(request.conditional_options, request.diagnostic_level, result))
    {
        return result;
    }

    try
    {
        validate_cqmom_moment_data_matches_hierarchy(request.moments, request.hierarchy,
                                                     {"CQMOM moment data layout must match the recursive hierarchy.",
                                                      "CQMOM moment data values must match the layout size.",
                                                      "CQMOM moment data values must be finite."});
    }
    catch (const std::exception& ex)
    {
        record_failure(result, RecursiveQuadratureStatus::invalid_input, {}, 0, ex.what());
        return result;
    }

    return compute_recursive_quadrature_from_validated_moment_data(request);
}

RecursiveQuadratureResult compute_recursive_quadrature_from_trusted_moment_data(const RecursiveQuadratureMomentDataRequest& request)
{
    RecursiveQuadratureResult result = create_recursive_result(request.hierarchy);
    result.quadrature.set_adaptive_branch_counts_enabled(request.wheeler_options.adaptive);

    if (!validate_recursive_request_preflight(request.conditional_options, request.diagnostic_level, result))
    {
        return result;
    }

    try
    {
        validate_cqmom_moment_values_match_layout(
            request.moments, {"", "CQMOM moment data values must match the layout size.", "CQMOM moment data values must be finite."});
    }
    catch (const std::exception& ex)
    {
        record_failure(result, RecursiveQuadratureStatus::invalid_input, {}, 0, ex.what());
        return result;
    }

    return compute_recursive_quadrature_from_validated_moment_data(request);
}

} // namespace internal

} // namespace qmomkit
