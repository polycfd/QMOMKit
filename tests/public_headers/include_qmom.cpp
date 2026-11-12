#include <qmomkit/qmom.hpp>
#include <type_traits>
#include <utility>
#include <vector>

static_assert(std::is_same_v<decltype(std::declval<const qmomkit::QmomState&>().maximum_node_count()), std::size_t>);
static_assert(std::is_same_v<decltype(std::declval<const qmomkit::QmomState&>().qmom_options()), const qmomkit::QmomOptions&>);

// The quadrature family is an option of its own, with the classical family as
// its default.
static_assert(std::is_same_v<decltype(qmomkit::QmomOptions::family), qmomkit::QmomFamily>);
static_assert(qmomkit::QmomOptions{}.family == qmomkit::QmomFamily::classical);

// The settings of the generalized family: N_G has no default, and the
// continuation and nu default to normal and 1.
static_assert(std::is_same_v<decltype(qmomkit::QmomOptions::generalized), qmomkit::GeneralizedOptions>);
static_assert(std::is_same_v<decltype(qmomkit::GeneralizedOptions::node_count), std::size_t>);
static_assert(qmomkit::GeneralizedOptions{}.node_count == 0);
static_assert(qmomkit::GeneralizedOptions{}.continuation == qmomkit::GeneralizedContinuation::normal);
static_assert(qmomkit::GeneralizedOptions{}.nu == 1.0);

// The status of a state candidate distinguishes success from failure.
static_assert(std::is_same_v<decltype(std::declval<const qmomkit::QmomStateCandidate&>().ok()), bool>);
static_assert(std::is_same_v<decltype(qmomkit::QmomStateCandidate::status), qmomkit::StateCandidateStatus>);

// The characteristic speeds carry a QMOM status.
static_assert(std::is_same_v<decltype(std::declval<const qmomkit::HyqmomCharacteristicSpeeds&>().ok()), bool>);
static_assert(std::is_same_v<decltype(qmomkit::HyqmomCharacteristicSpeeds::speeds), std::vector<double>>);
static_assert(std::is_same_v<decltype(qmomkit::HyqmomCharacteristicSpeeds::status), qmomkit::QmomStatus>);

namespace public_header_pins
{

// Each pin takes the address of one entry point through its exact signature, as
// in include_cqmom.cpp. The policy overload compiles with RealizabilityPolicy
// only declared, as this header declares it.

using QmomStateCandidateCreation = qmomkit::QmomStateCandidate (qmomkit::QmomCase::*)(const qmomkit::MomentSet&) const;

using QmomCheckedStateCandidateCreation = qmomkit::QmomStateCandidate (qmomkit::QmomCase::*)(const qmomkit::MomentSet&,
                                                                                             const qmomkit::RealizabilityPolicy&) const;

using QmomCaseInitialState = qmomkit::QmomState (qmomkit::QmomCase::*)() const;

using QmomInversion = qmomkit::QmomResult (*)(const qmomkit::MomentSet&, std::size_t, const qmomkit::QmomOptions&);

using QmomBoundedInversion = qmomkit::QmomResult (*)(const qmomkit::MomentSet&, std::size_t, const qmomkit::QmomOptions&,
                                                     std::optional<qmomkit::AbscissaBounds>);

using QmomCaseInversion = qmomkit::QmomResult (*)(const qmomkit::QmomCase&);

using HyqmomCharacteristicSpeedsOfState = qmomkit::HyqmomCharacteristicSpeeds (*)(const qmomkit::QmomState&);

QmomInversion qmom_inversion = qmomkit::compute_qmom;

QmomBoundedInversion qmom_bounded_inversion = qmomkit::compute_qmom;

QmomCaseInversion qmom_case_inversion = qmomkit::compute_qmom;

QmomStateCandidateCreation qmom_state_candidate_creation = &qmomkit::QmomCase::state_candidate;

QmomCheckedStateCandidateCreation qmom_checked_state_candidate_creation = &qmomkit::QmomCase::state_candidate;

QmomCaseInitialState qmom_case_initial_state = &qmomkit::QmomCase::initial_state;

HyqmomCharacteristicSpeedsOfState hyqmom_characteristic_speeds_of_state = qmomkit::compute_hyqmom_characteristic_speeds;

} // namespace public_header_pins
