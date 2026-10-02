#pragma once

#include "cqmom_case.hpp"
#include "qmomkit/moment_source.hpp"

namespace qmomkit::internal
{

class QmomSourceContextAccess
{
public:
    [[nodiscard]] static QmomSourceContext create(double time, const std::vector<double>& values, std::size_t reconstruction_node_count,
                                                  const QmomOptions& reconstruction_options, std::optional<AbscissaBounds> reconstruction_bounds);
};

class CqmomSourceContextAccess
{
public:
    [[nodiscard]] static CqmomSourceContext create(double time, const PhaseSpaceHierarchy& hierarchy, const MultiIndexTable& indices,
                                                   const std::vector<double>& values, const CqmomCase& reconstruction_cqmom_case,
                                                   const std::vector<CqmomActiveBranchCount>& active_branch_counts,
                                                   const PreparedCqmomMetadata* prepared_metadata);
};

// This entry point is for validated QMOM phase-space evolution only. The
// caller owns the source-stage reconstruction metadata; evolving values,
// source controls, callbacks, and callback rates are still validated.
[[nodiscard]] QmomMomentSourceUpdateResult
advance_trusted_qmom_moments_with_source(const std::vector<double>& moments, double t, double time_step, const QmomMomentSource& source,
                                         const MomentSourceUpdateOptions& options, std::size_t reconstruction_node_count,
                                         const QmomOptions& reconstruction_options, std::optional<AbscissaBounds> reconstruction_bounds);

[[nodiscard]] QmomMomentSourceUpdateResult
advance_trusted_qmom_moments_with_source(const std::vector<double>& moments, double t, double time_step, const QmomMomentSource& source,
                                         const QmomMomentBoundaryContribution& boundary, const MomentSourceUpdateOptions& options,
                                         std::size_t reconstruction_node_count, const QmomOptions& reconstruction_options,
                                         std::optional<AbscissaBounds> reconstruction_bounds);

// This entry point is for validated CqmomState phase-space evolution only. The caller owns
// the proof that the fixed source layout is immutable and valid; evolving
// values, source controls, callbacks, and callback rates are still validated.
[[nodiscard]] CqmomMomentSourceUpdateResult advance_trusted_cqmom_moments_with_source(const CqmomMomentData& moments, double t, double time_step,
                                                                                      const CqmomMomentSource& source,
                                                                                      const MomentSourceUpdateOptions& options);

[[nodiscard]] CqmomMomentSourceUpdateResult advance_trusted_cqmom_moments_with_source(const CqmomMomentData& moments, double t, double time_step,
                                                                                      const CqmomMomentSource& source,
                                                                                      const MomentSourceUpdateOptions& options,
                                                                                      const CqmomCase& reconstruction_cqmom_case,
                                                                                      const std::vector<CqmomActiveBranchCount>& active_branch_counts,
                                                                                      const PreparedCqmomMetadata* prepared_metadata);

[[nodiscard]] CqmomMomentSourceUpdateResult advance_trusted_cqmom_moments_with_source(const CqmomMomentData& moments, double t, double time_step,
                                                                                      const CqmomMomentSource& source,
                                                                                      const CqmomMomentBoundaryContribution& boundary,
                                                                                      const MomentSourceUpdateOptions& options);

[[nodiscard]] CqmomMomentSourceUpdateResult
advance_trusted_cqmom_moments_with_source(const CqmomMomentData& moments, double t, double time_step, const CqmomMomentSource& source,
                                          const CqmomMomentBoundaryContribution& boundary, const MomentSourceUpdateOptions& options,
                                          const CqmomCase& reconstruction_cqmom_case, const std::vector<CqmomActiveBranchCount>& active_branch_counts,
                                          const PreparedCqmomMetadata* prepared_metadata);

} // namespace qmomkit::internal
