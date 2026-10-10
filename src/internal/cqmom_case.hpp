#pragma once

#include "qmomkit/cqmom.hpp"

#include <algorithm>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace qmomkit
{
struct CqmomActiveBranchCount;
} // namespace qmomkit

namespace qmomkit::internal
{

struct PreparedCqmomMetadata;
using PreparedCqmomMetadataPtr = std::shared_ptr<const PreparedCqmomMetadata>;

// A successor case is the same case with new moment values: every other
// setting is carried unchanged. This helper is the one place that carries them,
// so a setting added to the case struct is added here as well.
//
// The settings are listed individually rather than carried by copying the whole
// case and then replacing its moments, because a wholesale copy also copies the
// moment values that are about to be replaced.
[[nodiscard]] inline CqmomCase with_moments(const CqmomCase& cqmom_case, CqmomMomentData moments)
{
    return CqmomCase{std::move(moments),          cqmom_case.conditional_options, cqmom_case.qmom_options,
                     cqmom_case.diagnostic_level, cqmom_case.vandermonde_solver,  cqmom_case.abscissa_bounds};
}

// The marginal moments M[0,...,k_d,...,0] of one coordinate, in ascending
// order of k_d, from values in the order of a CQMOM moment layout.
[[nodiscard]] inline std::vector<double> extract_coordinate_marginal(const std::vector<double>& values, const CqmomMomentLayout& layout,
                                                                     const std::size_t coordinate)
{
    std::vector<std::pair<std::size_t, std::size_t>> orders_and_positions;
    for (std::size_t position = 0; position < layout.required_indices.size(); ++position)
    {
        const MultiIndexView index = layout.required_indices[position];
        bool pure_coordinate = true;
        for (std::size_t dimension = 0; dimension < index.size(); ++dimension)
        {
            if (dimension != coordinate && index[dimension] != 0)
            {
                pure_coordinate = false;
                break;
            }
        }
        if (pure_coordinate)
        {
            orders_and_positions.emplace_back(index[coordinate], position);
        }
    }

    std::sort(orders_and_positions.begin(), orders_and_positions.end());
    std::vector<double> marginal;
    marginal.reserve(orders_and_positions.size());
    for (const auto& order_and_position : orders_and_positions)
    {
        marginal.push_back(values[order_and_position.second]);
    }
    return marginal;
}

class CqmomConstructionAccess
{
public:
    [[nodiscard]] static CqmomPreparedCase create_prepared_case(CqmomCase cqmom_case, PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomCase cqmom_case, PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomCase cqmom_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomCase cqmom_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 std::optional<CqmomQuadrature> accepted_quadrature,
                                                 PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomCase cqmom_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 CqmomQuadrature&& accepted_quadrature, PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static const CqmomCase& case_definition(const CqmomPreparedCase& cqmom_prepared_case);
    [[nodiscard]] static const CqmomCase& state_definition(const CqmomState& state);
    [[nodiscard]] static const std::vector<CqmomActiveBranchCount>& state_active_branch_counts(const CqmomState& state);
    [[nodiscard]] static const std::optional<CqmomQuadrature>& state_accepted_quadrature(const CqmomState& state);
    [[nodiscard]] static const PreparedCqmomMetadataPtr& case_prepared_metadata(const CqmomPreparedCase& cqmom_prepared_case);
    [[nodiscard]] static const PreparedCqmomMetadataPtr& state_prepared_metadata(const CqmomState& state);
    [[nodiscard]] static CqmomQuadrature create_quadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights,
                                                           std::vector<std::vector<double>> abscissas, bool adaptive_branch_counts_enabled,
                                                           bool rectangular_storage, ConditioningPathTable terminal_paths);
};

// Reconstruction of a QMOMKit-owned CQMOM case or state. Every prepared case
// and state carries prepared metadata that matches its hierarchy and
// adaptivity: callers either pass on the metadata of the same lineage or start
// a new one. Internal reconstruction rechecks evolving numerical values only.
// A successful reconstruction stores its active branch counts in
// `accepted_active_branch_counts` when that is given; they are empty for a
// fixed case.
[[nodiscard]] CqmomResult compute_cqmom_internal(const CqmomCase& cqmom_case, std::vector<CqmomActiveBranchCount> active_branch_counts = {},
                                                 bool allow_active_branch_count_growth = false,
                                                 const PreparedCqmomMetadata* prepared_metadata = nullptr, bool require_prepared_metadata = false,
                                                 std::vector<CqmomActiveBranchCount>* accepted_active_branch_counts = nullptr);

} // namespace qmomkit::internal
