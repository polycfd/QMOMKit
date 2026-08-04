#pragma once

#include "qmomkit/moment_rates.hpp"

#include <vector>

namespace qmomkit::internal
{

[[nodiscard]] QmomDiscontinuousMomentRateResult evaluate_qmom_discontinuous_moment_rates(double time, const std::vector<double>& moments,
                                                                                         const MomentSource& source);

[[nodiscard]] QmomDiscontinuousMomentRateResult evaluate_qmom_discontinuous_moment_rates(double time, const std::vector<double>& moments,
                                                                                         const MomentSource& source,
                                                                                         const MomentBoundaryContribution& boundary);

[[nodiscard]] CqmomDiscontinuousMomentRateResult evaluate_cqmom_discontinuous_moment_rates(double time, const CqmomMomentData& moments,
                                                                                           const CqmomMomentSource& source);

} // namespace qmomkit::internal
