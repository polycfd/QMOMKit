#pragma once

#include "qmomkit/realizability_policy.hpp"

namespace qmomkit::internal
{

// Why `policy` cannot be applied, or nullptr when it can. A policy whose
// check_mode is None is never applied and therefore always valid. The rules are
// those check_realizability applies to its options and interval.
[[nodiscard]] const char* realizability_policy_error(const RealizabilityPolicy& policy);

} // namespace qmomkit::internal
