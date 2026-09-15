#pragma once

#include <optional>

namespace qmomkit
{

enum class ConditionalMomentStatus
{
    success,
    invalid_input,
    singular_system
};

struct ConditionalMomentOptions
{
    double singular_tolerance = 1.0e-14;
    double conditioning_warning_threshold = 1.0e12;
    double moore_penrose_rcond = 1.0e-14;
};

struct ConditionalMomentDiagnostics
{
    /**
     * Expensive diagnostic quantities are reported when the caller requested
     * them. They are empty when diagnostic computation is disabled, so callers
     * should check availability before reading the values.
     */
    std::optional<double> residual_norm;
    std::optional<double> condition_estimate;
    bool conditioning_warning = false;
};

} // namespace qmomkit
