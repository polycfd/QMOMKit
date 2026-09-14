#pragma once

#include "numeric_utils.hpp"
#include "qmomkit/runge_kutta.hpp"

#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace qmomkit::internal
{

#ifndef NDEBUG
[[nodiscard]] inline std::string ode_rhs_unfilled_message(const double time, const std::size_t component)
{
    std::ostringstream output;
    output << "ODE right-hand side did not fill every derivative entry:\n"
           << "  time: " << time << '\n'
           << "  component: " << component;
    return output.str();
}
#endif

/**
 * Evaluate the right-hand side into `derivative` and check what it wrote.
 *
 * The buffer is sized to the state and handed to the callback, so a buffer that
 * has already served an earlier evaluation is reused rather than reallocated.
 * Both Runge-Kutta drivers validate the same three things in the same order:
 * that the callback did not resize the buffer, that it wrote every entry, and
 * that every value is finite.
 */
inline void rhs_checked(const OdeRhs& rhs, const double time, const std::vector<double>& state, std::vector<double>& derivative)
{
    // Zeroed, not merely resized: a callback that violates the fill-every-entry
    // contract then degrades to a zero derivative rather than silently
    // integrating whatever the previous stage left in the reused buffer. The
    // extra pass is negligible beside the callback, and debug builds overwrite
    // it immediately with the unwritten-entry marking.
    derivative.assign(state.size(), 0.0);
    poison_output_buffer(derivative);
    rhs(time, state, derivative);

    if (derivative.size() != state.size())
    {
        throw std::invalid_argument("RHS returned a vector with the wrong size.");
    }
#ifndef NDEBUG
    // Ahead of the finiteness scan so that an unwritten entry is reported as the
    // missing write it is, rather than as a non-finite value the callback never
    // produced.
    for (std::size_t component = 0; component < derivative.size(); ++component)
    {
        if (is_output_buffer_poison(derivative[component]))
        {
            throw std::invalid_argument(ode_rhs_unfilled_message(time, component));
        }
    }
#endif
    if (!all_finite(derivative))
    {
        throw std::invalid_argument("RHS returned a non-finite value.");
    }
}

} // namespace qmomkit::internal
