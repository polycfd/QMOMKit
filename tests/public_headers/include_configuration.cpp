#include <qmomkit/configuration.hpp>

namespace
{

using CqmomInitialMomentReplacement = qmomkit::config::ConfigurationResult<qmomkit::CqmomPreparedCase> (*)(const qmomkit::CqmomPreparedCase&,
                                                                                                           qmomkit::RawMomentTensor);

using CqmomLeanInitialMomentReplacement = qmomkit::config::ConfigurationResult<qmomkit::CqmomPreparedCase> (*)(const qmomkit::CqmomPreparedCase&,
                                                                                                               qmomkit::CqmomMomentData);

[[maybe_unused]] CqmomInitialMomentReplacement cqmom_initial_moment_replacement = qmomkit::config::with_initial_moments;

[[maybe_unused]] CqmomLeanInitialMomentReplacement cqmom_lean_initial_moment_replacement = qmomkit::config::with_initial_moments;

} // namespace
