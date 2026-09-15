#pragma once

#include "qmomkit/cqmom.hpp"

#include <memory>
#include <optional>
#include <vector>

namespace qmomkit
{
struct CqmomActiveBranchCount;
} // namespace qmomkit

namespace qmomkit::internal
{

struct PreparedCqmomMetadata;
using PreparedCqmomMetadataPtr = std::shared_ptr<const PreparedCqmomMetadata>;

class CqmomConstructionAccess
{
public:
    [[nodiscard]] static CqmomCase create_case(CqmomDenseCase dense_case, PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomCase create_case(CqmomLeanCase lean_case, PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomDenseCase dense_case, PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomDenseCase dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomDenseCase dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 std::optional<CqmomQuadrature> accepted_quadrature,
                                                 PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomDenseCase dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 CqmomQuadrature&& accepted_quadrature, PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case, PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 std::optional<CqmomQuadrature> accepted_quadrature,
                                                 PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 CqmomQuadrature&& accepted_quadrature, PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static const CqmomDenseCase* case_dense_case(const CqmomCase& case_);
    [[nodiscard]] static const CqmomLeanCase* case_lean_case(const CqmomCase& case_);
    [[nodiscard]] static const CqmomDenseCase* state_dense_case(const CqmomState& state);
    [[nodiscard]] static const CqmomLeanCase* state_lean_case(const CqmomState& state);
    [[nodiscard]] static const std::vector<CqmomActiveBranchCount>& state_active_branch_counts(const CqmomState& state);
    [[nodiscard]] static const std::optional<CqmomQuadrature>& state_accepted_quadrature(const CqmomState& state);
    [[nodiscard]] static const PreparedCqmomMetadataPtr& case_prepared_metadata(const CqmomCase& case_);
    [[nodiscard]] static const PreparedCqmomMetadataPtr& state_prepared_metadata(const CqmomState& state);
    [[nodiscard]] static CqmomQuadrature create_quadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights,
                                                           std::vector<std::vector<double>> abscissas, bool adaptive_branch_counts_enabled,
                                                           bool rectangular_storage, ConditioningPathTable terminal_paths);
};

// Internal construction access for QMOMKit-owned CQMOM objects. Fixed normal
// cases and states always receive a compatible eagerly prepared owner: callers
// either propagate one from the same lineage or start a new lineage. Adaptive
// states carry no fixed prepared owner. Internal reconstruction rechecks
// evolving numerical values only.
[[nodiscard]] CqmomResult compute_cqmom_internal(const CqmomDenseCase& dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts = {},
                                                 bool allow_active_branch_count_growth = false,
                                                 const PreparedCqmomMetadata* prepared_metadata = nullptr, bool require_prepared_metadata = false);
[[nodiscard]] CqmomResult compute_cqmom_internal(const CqmomLeanCase& lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts = {},
                                                 bool allow_active_branch_count_growth = false,
                                                 const PreparedCqmomMetadata* prepared_metadata = nullptr, bool require_prepared_metadata = false);

} // namespace qmomkit::internal
