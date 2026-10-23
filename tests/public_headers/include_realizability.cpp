#include <qmomkit/realizability.hpp>

namespace public_header_pins
{

// Each pin takes the address of one entry point through its exact signature, as
// in include_cqmom.cpp: a changed parameter or return type fails to compile, and
// because the pins have external linkage, a missing definition fails to link.

using SupportRealizabilityCheck = qmomkit::RealizabilityResult (*)(const qmomkit::MomentSet&, qmomkit::MomentSupport,
                                                                   const qmomkit::RealizabilityOptions&);

using IntervalRealizabilityCheck = qmomkit::RealizabilityResult (*)(const qmomkit::MomentSet&, qmomkit::BoundedInterval,
                                                                    const qmomkit::RealizabilityOptions&);

SupportRealizabilityCheck support_realizability_check = qmomkit::check_realizability;

IntervalRealizabilityCheck interval_realizability_check = qmomkit::check_realizability;

} // namespace public_header_pins
