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

struct PreparedFixedCqmomMetadata;
using PreparedFixedCqmomMetadataPtr = std::shared_ptr<const PreparedFixedCqmomMetadata>;

class CqmomConstructionAccess
{
public:
    [[nodiscard]] static CqmomCase create_case(CqmomDenseCase dense_case, PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
    [[nodiscard]] static CqmomCase create_case(CqmomLeanCase lean_case, PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomDenseCase dense_case, PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomDenseCase dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomDenseCase dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 std::optional<CqmomQuadrature> accepted_quadrature,
                                                 PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomDenseCase dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 CqmomQuadrature&& accepted_quadrature,
                                                 PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case, PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 std::optional<CqmomQuadrature> accepted_quadrature,
                                                 PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 CqmomQuadrature&& accepted_quadrature,
                                                 PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
    [[nodiscard]] static const CqmomDenseCase* case_dense_case(const CqmomCase& case_);
    [[nodiscard]] static const CqmomLeanCase* case_lean_case(const CqmomCase& case_);
    [[nodiscard]] static const CqmomDenseCase* state_dense_case(const CqmomState& state);
    [[nodiscard]] static const CqmomLeanCase* state_lean_case(const CqmomState& state);
    [[nodiscard]] static const std::vector<CqmomActiveBranchCount>& state_active_branch_counts(const CqmomState& state);
    [[nodiscard]] static const std::optional<CqmomQuadrature>& state_accepted_quadrature(const CqmomState& state);
    [[nodiscard]] static const PreparedFixedCqmomMetadataPtr& case_prepared_fixed_metadata(const CqmomCase& case_);
    [[nodiscard]] static const PreparedFixedCqmomMetadataPtr& state_prepared_fixed_metadata(const CqmomState& state);
    [[nodiscard]] static CqmomQuadrature create_quadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights,
                                                           std::vector<std::vector<double>> abscissas, bool adaptive_branch_counts_enabled,
                                                           bool rectangular_storage, std::vector<ConditioningPath> terminal_paths);
};

// Internal construction access for QMOMKit-owned CQMOM objects. Fixed normal
// cases and states always receive a compatible eagerly prepared owner: callers
// either propagate one from the same lineage or start a new lineage. Adaptive
// states carry no fixed prepared owner. Trusted reconstruction rechecks
// evolving numerical values only.
[[nodiscard]] CqmomCase create_cqmom_case(CqmomDenseCase dense_case, PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
[[nodiscard]] CqmomCase create_cqmom_case(CqmomLeanCase lean_case, PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
[[nodiscard]] CqmomState create_cqmom_state(CqmomDenseCase dense_case, PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
[[nodiscard]] CqmomState create_cqmom_state(CqmomDenseCase dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                            PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
[[nodiscard]] CqmomState create_cqmom_state(CqmomDenseCase dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                            std::optional<CqmomQuadrature> accepted_quadrature,
                                            PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
[[nodiscard]] CqmomState create_cqmom_state(CqmomDenseCase dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                            CqmomQuadrature&& accepted_quadrature, PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
[[nodiscard]] CqmomState create_cqmom_state(CqmomLeanCase lean_case, PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
[[nodiscard]] CqmomState create_cqmom_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                            PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
[[nodiscard]] CqmomState create_cqmom_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                            std::optional<CqmomQuadrature> accepted_quadrature,
                                            PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
[[nodiscard]] CqmomState create_cqmom_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                            CqmomQuadrature&& accepted_quadrature, PreparedFixedCqmomMetadataPtr prepared_fixed_metadata = nullptr);
[[nodiscard]] const CqmomDenseCase* cqmom_case_dense_case(const CqmomCase& case_);
[[nodiscard]] const CqmomLeanCase* cqmom_case_lean_case(const CqmomCase& case_);
[[nodiscard]] const CqmomDenseCase* cqmom_state_dense_case(const CqmomState& state);
[[nodiscard]] const CqmomLeanCase* cqmom_state_lean_case(const CqmomState& state);
[[nodiscard]] const std::vector<CqmomActiveBranchCount>& cqmom_state_active_branch_counts(const CqmomState& state);
[[nodiscard]] const std::optional<CqmomQuadrature>& cqmom_state_accepted_quadrature(const CqmomState& state);
[[nodiscard]] const PreparedFixedCqmomMetadataPtr& cqmom_case_prepared_fixed_metadata(const CqmomCase& case_);
[[nodiscard]] const PreparedFixedCqmomMetadataPtr& cqmom_state_prepared_fixed_metadata(const CqmomState& state);
[[nodiscard]] CqmomResult compute_trusted_cqmom(const CqmomDenseCase& dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts = {},
                                                bool allow_active_branch_count_growth = false,
                                                const PreparedFixedCqmomMetadata* prepared_fixed_metadata = nullptr,
                                                bool require_prepared_fixed_metadata = false);
[[nodiscard]] CqmomResult compute_trusted_cqmom(const CqmomLeanCase& lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts = {},
                                                bool allow_active_branch_count_growth = false,
                                                const PreparedFixedCqmomMetadata* prepared_fixed_metadata = nullptr,
                                                bool require_prepared_fixed_metadata = false);
[[nodiscard]] CqmomQuadrature create_cqmom_quadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights,
                                                      std::vector<std::vector<double>> abscissas, bool adaptive_branch_counts_enabled,
                                                      bool rectangular_storage, std::vector<ConditioningPath> terminal_paths);

} // namespace qmomkit::internal
