#pragma once

#include "qmomkit/moment_rates.hpp"

#include <vector>

namespace qmomkit::internal
{

[[nodiscard]] std::vector<std::vector<double>> evaluate_terminal_velocities(const PhaseSpaceVelocity& velocity, double time,
                                                                            const std::vector<std::vector<double>>& terminal_abscissas);

[[nodiscard]] QmomAdvectionMomentRateResult evaluate_qmom_advection_moment_rates(const std::vector<double>& moments, std::size_t node_count,
                                                                                 double time, const PhaseSpaceVelocity& velocity,
                                                                                 const QmomOptions& qmom_options = {});

[[nodiscard]] CqmomAdvectionMomentRateResult evaluate_cqmom_advection_moment_rates(const CqmomLeanCase& lean_case, double time,
                                                                                   const PhaseSpaceVelocity& velocity);

} // namespace qmomkit::internal
