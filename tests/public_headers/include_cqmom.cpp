#include <qmomkit/cqmom.hpp>
#include <type_traits>

// The open case is filled positionally by brace initialization throughout the
// library, examples and tests, so it must stay an aggregate.
static_assert(std::is_aggregate_v<qmomkit::CqmomCase>);

// A prepared case is built from an open case only on request. Allowing the
// conversion implicitly would make overloads that accept either case ambiguous.
static_assert(std::is_constructible_v<qmomkit::CqmomPreparedCase, qmomkit::CqmomCase>);
static_assert(!std::is_convertible_v<qmomkit::CqmomCase, qmomkit::CqmomPreparedCase>);

namespace public_header_pins
{

// Each pin takes the address of one entry point through its exact signature. A
// changed parameter or return type fails to compile here. The pins have
// external linkage, so each one is emitted and its function must be defined
// when this test is linked. In an unnamed namespace the compiler may drop them,
// as Clang does even without optimization, and the definitions would go
// unchecked.

using CqmomCaseInversion = qmomkit::CqmomResult (*)(const qmomkit::CqmomCase&);

using CqmomPreparedCaseInversion = qmomkit::CqmomResult (*)(const qmomkit::CqmomPreparedCase&);

using CqmomMomentExtraction = qmomkit::CqmomMomentData (*)(const qmomkit::RawMomentTensor&, const qmomkit::CqmomMomentLayout&);

using CqmomCaseInitialState = qmomkit::CqmomState (qmomkit::CqmomCase::*)() const;

using CqmomPreparedCaseInitialState = qmomkit::CqmomState (qmomkit::CqmomPreparedCase::*)() const;

CqmomCaseInversion cqmom_case_inversion = qmomkit::compute_cqmom;

CqmomPreparedCaseInversion cqmom_prepared_case_inversion = qmomkit::compute_cqmom;

CqmomMomentExtraction cqmom_moment_extraction = qmomkit::extract_cqmom_moments;

CqmomCaseInitialState cqmom_case_initial_state = &qmomkit::CqmomCase::initial_state;

CqmomPreparedCaseInitialState cqmom_prepared_case_initial_state = &qmomkit::CqmomPreparedCase::initial_state;

using CqmomStateCandidateCreation = qmomkit::CqmomStateCandidate (qmomkit::CqmomPreparedCase::*)(const qmomkit::CqmomMomentData&) const;

using CqmomCheckedStateCandidateCreation = qmomkit::CqmomStateCandidate (qmomkit::CqmomPreparedCase::*)(
    const qmomkit::CqmomMomentData&, const std::vector<std::optional<qmomkit::RealizabilityPolicy>>&) const;

CqmomStateCandidateCreation cqmom_state_candidate_creation = &qmomkit::CqmomPreparedCase::state_candidate;

CqmomCheckedStateCandidateCreation cqmom_checked_state_candidate_creation = &qmomkit::CqmomPreparedCase::state_candidate;

} // namespace public_header_pins
