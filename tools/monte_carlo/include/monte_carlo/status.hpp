#pragma once

namespace qmomkit::monte_carlo
{

/**
 * Shared status for comparison-oriented utilities.
 *
 * The status distinguishes malformed inputs from failures in the deterministic
 * QMOMKit path, the Monte Carlo sampled path, or the final moment comparison.
 */
enum class ComparisonStatus
{
    success,
    invalid_input,
    sampling_failed,
    qmomkit_failed,
    monte_carlo_failed,
    moment_failed
};

} // namespace qmomkit::monte_carlo
