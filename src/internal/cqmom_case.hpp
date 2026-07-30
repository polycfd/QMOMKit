#pragma once

#include "qmomkit/cqmom.hpp"

#include <optional>
#include <vector>

namespace qmomkit
{
struct CqmomActiveBranchCount;
} // namespace qmomkit

namespace qmomkit::internal
{

class CqmomConstructionAccess
{
public:
    [[nodiscard]] static CqmomCase create_case(CqmomDenseCase dense_case);
    [[nodiscard]] static CqmomCase create_case(CqmomLeanCase lean_case);
    [[nodiscard]] static CqmomState create_state(CqmomDenseCase dense_case);
    [[nodiscard]] static CqmomState create_state(CqmomDenseCase dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts);
    [[nodiscard]] static CqmomState create_state(CqmomDenseCase dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 std::optional<CqmomQuadrature> accepted_quadrature);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 std::optional<CqmomQuadrature> accepted_quadrature);
    [[nodiscard]] static const CqmomDenseCase* state_dense_case(const CqmomState& state);
    [[nodiscard]] static const CqmomLeanCase* state_lean_case(const CqmomState& state);
    [[nodiscard]] static const std::vector<CqmomActiveBranchCount>& state_active_branch_counts(const CqmomState& state);
    [[nodiscard]] static const std::optional<CqmomQuadrature>& state_accepted_quadrature(const CqmomState& state);
    [[nodiscard]] static CqmomQuadrature create_quadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights,
                                                           std::vector<std::vector<double>> abscissas, bool adaptive_branch_counts_enabled,
                                                           bool rectangular_storage, std::vector<ConditioningPath> terminal_paths);
};

// Internal construction access for QMOMKit-owned CQMOM objects. Callers must
// pass metadata that has already crossed the normal parser/API validation
// boundary; trusted reconstruction relies on that invariant and rechecks only
// evolving numerical values.
[[nodiscard]] CqmomCase create_cqmom_case(CqmomDenseCase dense_case);
[[nodiscard]] CqmomCase create_cqmom_case(CqmomLeanCase lean_case);
[[nodiscard]] CqmomState create_cqmom_state(CqmomDenseCase dense_case);
[[nodiscard]] CqmomState create_cqmom_state(CqmomDenseCase dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts);
[[nodiscard]] CqmomState create_cqmom_state(CqmomDenseCase dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                            std::optional<CqmomQuadrature> accepted_quadrature);
[[nodiscard]] CqmomState create_cqmom_state(CqmomLeanCase lean_case);
[[nodiscard]] CqmomState create_cqmom_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts);
[[nodiscard]] CqmomState create_cqmom_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                            std::optional<CqmomQuadrature> accepted_quadrature);
[[nodiscard]] const CqmomDenseCase* cqmom_state_dense_case(const CqmomState& state);
[[nodiscard]] const CqmomLeanCase* cqmom_state_lean_case(const CqmomState& state);
[[nodiscard]] const std::vector<CqmomActiveBranchCount>& cqmom_state_active_branch_counts(const CqmomState& state);
[[nodiscard]] const std::optional<CqmomQuadrature>& cqmom_state_accepted_quadrature(const CqmomState& state);
[[nodiscard]] CqmomResult compute_trusted_cqmom(const CqmomDenseCase& dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts = {},
                                                bool allow_active_branch_count_growth = false);
[[nodiscard]] CqmomResult compute_trusted_cqmom(const CqmomLeanCase& lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts = {},
                                                bool allow_active_branch_count_growth = false);
[[nodiscard]] CqmomQuadrature create_cqmom_quadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights,
                                                      std::vector<std::vector<double>> abscissas, bool adaptive_branch_counts_enabled,
                                                      bool rectangular_storage, std::vector<ConditioningPath> terminal_paths);

} // namespace qmomkit::internal
