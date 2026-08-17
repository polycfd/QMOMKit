#pragma once

#include "cqmom_case.hpp"
#include "qmomkit/moment_source.hpp"

namespace qmomkit::internal
{

class CqmomSourceContextAccess
{
public:
    [[nodiscard]] static CqmomSourceContext create(double time, const PhaseSpaceHierarchy& hierarchy, const std::vector<MultiIndex>& indices,
                                                   const std::vector<double>& values, const CqmomLeanCase& reconstruction_case,
                                                   const std::vector<CqmomActiveBranchCount>& active_branch_counts,
                                                   const PreparedFixedCqmomMetadata* prepared_fixed_metadata);
};

// This entry point is for validated CqmomState phase-space evolution only. The caller owns
// the proof that the fixed source layout is immutable and valid; evolving
// values, source controls, callbacks, and callback rates are still validated.
[[nodiscard]] CqmomMomentSourceUpdateResult advance_trusted_cqmom_moments_with_source(const CqmomMomentData& moments, double t, double time_step,
                                                                                      const CqmomMomentSource& source,
                                                                                      const MomentSourceUpdateOptions& options);

[[nodiscard]] CqmomMomentSourceUpdateResult advance_trusted_cqmom_moments_with_source(const CqmomMomentData& moments, double t, double time_step,
                                                                                      const CqmomMomentSource& source,
                                                                                      const MomentSourceUpdateOptions& options,
                                                                                      const CqmomLeanCase& reconstruction_case,
                                                                                      const std::vector<CqmomActiveBranchCount>& active_branch_counts,
                                                                                      const PreparedFixedCqmomMetadata* prepared_fixed_metadata);

[[nodiscard]] CqmomMomentSourceUpdateResult advance_trusted_cqmom_moments_with_source(const CqmomMomentData& moments, double t, double time_step,
                                                                                      const CqmomMomentSource& source,
                                                                                      const CqmomMomentBoundaryContribution& boundary,
                                                                                      const MomentSourceUpdateOptions& options);

[[nodiscard]] CqmomMomentSourceUpdateResult
advance_trusted_cqmom_moments_with_source(const CqmomMomentData& moments, double t, double time_step, const CqmomMomentSource& source,
                                          const CqmomMomentBoundaryContribution& boundary, const MomentSourceUpdateOptions& options,
                                          const CqmomLeanCase& reconstruction_case, const std::vector<CqmomActiveBranchCount>& active_branch_counts,
                                          const PreparedFixedCqmomMetadata* prepared_fixed_metadata);

} // namespace qmomkit::internal
