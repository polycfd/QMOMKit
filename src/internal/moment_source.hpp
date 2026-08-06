#pragma once

#include "qmomkit/moment_source.hpp"

namespace qmomkit::internal
{

// This entry point is for validated CqmomState phase-space evolution only. The caller owns
// the proof that the fixed source layout is immutable and valid; evolving
// values, source controls, callbacks, and callback rates are still validated.
[[nodiscard]] CqmomMomentSourceUpdateResult advance_trusted_cqmom_moments_with_source(const CqmomMomentData& moments, double t, double time_step,
                                                                                      const CqmomMomentSource& source,
                                                                                      const MomentSourceUpdateOptions& options);

[[nodiscard]] CqmomMomentSourceUpdateResult advance_trusted_cqmom_moments_with_source(const CqmomMomentData& moments, double t, double time_step,
                                                                                      const CqmomMomentSource& source,
                                                                                      const CqmomMomentBoundaryContribution& boundary,
                                                                                      const MomentSourceUpdateOptions& options);

} // namespace qmomkit::internal
