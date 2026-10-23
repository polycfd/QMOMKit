#include <qmomkit/wheeler.hpp>

namespace public_header_pins
{

// Each pin takes the address of one entry point through its exact signature, as
// in include_cqmom.cpp: a changed parameter or return type fails to compile, and
// because the pins have external linkage, a missing definition fails to link.

using QuadratureInversion = qmomkit::Quadrature1D (*)(const qmomkit::MomentSet&, std::size_t, const qmomkit::WheelerOptions&);

QuadratureInversion quadrature_inversion = qmomkit::compute_quadrature;

} // namespace public_header_pins
