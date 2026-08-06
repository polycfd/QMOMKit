#pragma once

#include "qmomkit/moment_rates.hpp"

#include <vector>

namespace qmomkit::internal
{

[[nodiscard]] QmomMomentSourceRateResult evaluate_qmom_moment_source_rates(double time, const std::vector<double>& moments,
                                                                           const MomentSource& source);

[[nodiscard]] QmomMomentSourceRateResult evaluate_qmom_moment_source_rates(double time, const std::vector<double>& moments,
                                                                           const MomentSource& source, const MomentBoundaryContribution& boundary);

[[nodiscard]] CqmomMomentSourceRateResult evaluate_cqmom_moment_source_rates(double time, const CqmomMomentData& moments,
                                                                             const CqmomMomentSource& source);

} // namespace qmomkit::internal
