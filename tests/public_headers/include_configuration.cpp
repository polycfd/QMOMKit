#include <qmomkit/configuration.hpp>

namespace
{

using CqmomInitialMomentReplacement = qmomkit::config::ConfigurationResult<qmomkit::CqmomCase> (*)(const qmomkit::CqmomCase&,
                                                                                                   qmomkit::RawMomentTensor);

[[maybe_unused]] CqmomInitialMomentReplacement cqmom_initial_moment_replacement = qmomkit::config::with_initial_moments;

} // namespace
