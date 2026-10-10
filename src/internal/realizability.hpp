#pragma once

#include "qmomkit/realizability_policy.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace qmomkit::internal
{

// Why `policy` cannot be applied, or nullptr when it can. A policy whose
// check_mode is None is never applied and therefore always valid. The rules are
// those check_realizability applies to its options and interval.
[[nodiscard]] const char* realizability_policy_error(const RealizabilityPolicy& policy);

// Why CQMOM coordinate policies cannot be applied, or an empty string when they
// can. A policy beyond `dimensions` is invalid whatever its mode, and is
// reported before any policy within them, which realizability_policy_error
// checks.
[[nodiscard]] std::string coordinate_realizability_policies_error(const std::vector<std::optional<RealizabilityPolicy>>& policies,
                                                                  std::size_t dimensions);

} // namespace qmomkit::internal
