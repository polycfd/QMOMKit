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

// A successor case is the same case with new moment values: every other
// setting is carried unchanged. This helper is the one place that knows how to
// carry them, so the four sites that build a successor no longer each re-list
// the settings. Adding a setting to the case struct means adding it here as
// well, which is one place to miss rather than four.
//
// The settings are listed individually rather than carried by copying the whole
// case and then replacing its moments: a wholesale copy also duplicates the
// outgoing moment values that are about to be replaced, which measured at
// +12 allocations per phase-space step.
[[nodiscard]] inline CqmomLeanCase with_moments(const CqmomLeanCase& case_, CqmomMomentData moments)
{
    return CqmomLeanCase{std::move(moments),     case_.conditional_options, case_.qmom_options,      case_.adaptive_branch_counts,
                         case_.diagnostic_level, case_.vandermonde_solver,  case_.coordinate_support};
}

class CqmomConstructionAccess
{
public:
    [[nodiscard]] static CqmomPreparedCase create_prepared_case(CqmomLeanCase lean_case, PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case, PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 std::optional<CqmomQuadrature> accepted_quadrature,
                                                 PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static CqmomState create_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 CqmomQuadrature&& accepted_quadrature, PreparedCqmomMetadataPtr prepared_metadata = nullptr);
    [[nodiscard]] static const CqmomLeanCase& case_definition(const CqmomPreparedCase& case_);
    [[nodiscard]] static const CqmomLeanCase& state_definition(const CqmomState& state);
    [[nodiscard]] static const std::vector<CqmomActiveBranchCount>& state_active_branch_counts(const CqmomState& state);
    [[nodiscard]] static const std::optional<CqmomQuadrature>& state_accepted_quadrature(const CqmomState& state);
    [[nodiscard]] static const PreparedCqmomMetadataPtr& case_prepared_metadata(const CqmomPreparedCase& case_);
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
[[nodiscard]] CqmomResult compute_cqmom_internal(const CqmomLeanCase& lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts = {},
                                                 bool allow_active_branch_count_growth = false,
                                                 const PreparedCqmomMetadata* prepared_metadata = nullptr, bool require_prepared_metadata = false);

} // namespace qmomkit::internal
