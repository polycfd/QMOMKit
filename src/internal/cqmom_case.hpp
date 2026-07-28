#pragma once

#include "qmomkit/cqmom.hpp"

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
    [[nodiscard]] static CqmomCase create_case(CqmomLeanCase lean_case);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts);
    [[nodiscard]] static const CqmomLeanCase& state_case(const CqmomState& state);
    [[nodiscard]] static const std::vector<CqmomActiveBranchCount>& state_active_branch_counts(const CqmomState& state);
    [[nodiscard]] static CqmomQuadrature create_quadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights,
                                                           std::vector<std::vector<double>> abscissas, bool adaptive_branch_counts_enabled,
                                                           bool rectangular_storage, std::vector<ConditioningPath> terminal_paths);
};

// Internal construction access for QMOMKit-owned CQMOM objects. Callers must
// pass metadata that has already crossed the normal parser/API validation
// boundary; trusted reconstruction relies on that invariant and rechecks only
// evolving numerical values.
[[nodiscard]] CqmomCase create_cqmom_case(CqmomLeanCase lean_case);
[[nodiscard]] CqmomState create_cqmom_state(CqmomLeanCase lean_case);
[[nodiscard]] CqmomState create_cqmom_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts);
[[nodiscard]] const CqmomLeanCase& cqmom_state_case(const CqmomState& state);
[[nodiscard]] const std::vector<CqmomActiveBranchCount>& cqmom_state_active_branch_counts(const CqmomState& state);
[[nodiscard]] CqmomResult compute_trusted_cqmom(const CqmomDenseCase& dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts = {});
[[nodiscard]] CqmomResult compute_trusted_cqmom(const CqmomLeanCase& lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts = {});
[[nodiscard]] CqmomQuadrature create_cqmom_quadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights,
                                                      std::vector<std::vector<double>> abscissas, bool adaptive_branch_counts_enabled,
                                                      bool rectangular_storage, std::vector<ConditioningPath> terminal_paths);

} // namespace qmomkit::internal
