#include <qmomkit/cqmom.hpp>
#include <type_traits>

// The open case is filled positionally by brace initialization throughout the
// library, examples and tests, so it must stay an aggregate.
static_assert(std::is_aggregate_v<qmomkit::CqmomLeanCase>);

// A prepared case is built from an open case only on request. Allowing the
// conversion implicitly would make overloads that accept either case ambiguous.
static_assert(std::is_constructible_v<qmomkit::CqmomPreparedCase, qmomkit::CqmomLeanCase>);
static_assert(!std::is_convertible_v<qmomkit::CqmomLeanCase, qmomkit::CqmomPreparedCase>);

namespace public_header_pins
{

// Each pin takes the address of one entry point through its exact signature. A
// changed parameter or return type fails to compile here. The pins have
// external linkage, so each one is emitted and its function must be defined
// when this test is linked. In an unnamed namespace the compiler may drop them,
// as Clang does even without optimization, and the definitions would go
// unchecked.

using CqmomLeanCaseInversion = qmomkit::CqmomResult (*)(const qmomkit::CqmomLeanCase&);

using CqmomPreparedCaseInversion = qmomkit::CqmomResult (*)(const qmomkit::CqmomPreparedCase&);

using CqmomMomentExtraction = qmomkit::CqmomMomentData (*)(const qmomkit::RawMomentTensor&, const qmomkit::CqmomMomentLayout&);

using CqmomPreparedCaseInitialState = qmomkit::CqmomState (qmomkit::CqmomPreparedCase::*)() const;

CqmomLeanCaseInversion cqmom_lean_case_inversion = qmomkit::compute_cqmom;

CqmomPreparedCaseInversion cqmom_prepared_case_inversion = qmomkit::compute_cqmom;

CqmomMomentExtraction cqmom_moment_extraction = qmomkit::extract_cqmom_moments;

CqmomPreparedCaseInitialState cqmom_prepared_case_initial_state = &qmomkit::CqmomPreparedCase::initial_state;

} // namespace public_header_pins
