#pragma once

#include "qmomkit/cqmom.hpp"

#include <vector>

namespace qmomkit::internal
{

// Reconstruct moments for QMOMKit-owned CQMOM state transport. The supplied
// layout is assumed to be the already validated fixed layout for the quadrature
// hierarchy; the implementation still checks quadrature storage and generated
// finite values.
[[nodiscard]] CqmomMomentData reconstruct_trusted_cqmom_moments(const CqmomQuadrature& quadrature, const CqmomMomentLayout& layout);

// Reconstruct raw-moment rates from terminal quadrature storage and already
// evaluated terminal-node velocities. Abscissas are stored by phase-space
// dimension, velocities by terminal node.
[[nodiscard]] std::vector<double> reconstruct_continuous_moment_rates(const std::vector<double>& weights,
                                                                      const std::vector<std::vector<double>>& abscissas,
                                                                      const std::vector<std::vector<double>>& velocities,
                                                                      const std::vector<MultiIndex>& indices);

} // namespace qmomkit::internal
