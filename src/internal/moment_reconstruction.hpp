#pragma once

#include "qmomkit/cqmom.hpp"

namespace qmomkit::internal
{

// Reconstruct moments for QMOMKit-owned CQMOM state transport. The supplied
// layout is assumed to be the already validated fixed layout for the quadrature
// hierarchy; the implementation still checks quadrature storage and generated
// finite values.
[[nodiscard]] CqmomMomentData reconstruct_trusted_cqmom_moments(const CqmomQuadrature& quadrature, const CqmomMomentLayout& layout);

} // namespace qmomkit::internal
