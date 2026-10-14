#include "monte_carlo/evolution.hpp"

#include "ensemble_checks.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace qmomkit::monte_carlo
{
namespace
{

[[nodiscard]] bool all_finite(const std::vector<double>& values)
{
    return std::all_of(values.begin(), values.end(), [](const double value) { return std::isfinite(value); });
}

[[nodiscard]] std::string format_diagnostic_vector(const std::vector<double>& values)
{
    std::ostringstream output;
    output << std::setprecision(std::numeric_limits<double>::max_digits10) << '[';
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (i > 0)
        {
            output << ", ";
        }
        output << values[i];
    }
    output << ']';
    return output.str();
}

[[nodiscard]] std::string phase_space_velocity_wrong_dimension_message(const std::size_t expected_size, const std::size_t returned_size)
{
    std::ostringstream output;
    output << "Phase-space velocity callback returned the wrong dimension:\n"
           << "  expected size: " << expected_size << '\n'
           << "  returned size: " << returned_size;
    return output.str();
}

[[nodiscard]] std::string phase_space_velocity_non_finite_message(const double time, const std::size_t dimension, const double value,
                                                                  const std::vector<double>& state, const std::vector<double>& velocity)
{
    std::ostringstream output;
    output << std::setprecision(std::numeric_limits<double>::max_digits10) << "Phase-space velocity callback returned a non-finite value:\n"
           << "  time: " << time << '\n'
           << "  dimension: " << dimension << '\n'
           << "  value: " << value << '\n'
           << "  state: " << format_diagnostic_vector(state) << '\n'
           << "  velocity: " << format_diagnostic_vector(velocity);
    return output.str();
}

// Mirrors the library's debug-only under-fill check. Samples are advanced
// through the same packed callback, which fills a buffer supplied here rather
// than returning one, so a callback that writes only part of it would leave the
// remaining coordinates advancing on whatever was already in the buffer.
// Neither the dimension check nor the finiteness check below can see that: the
// buffer is the right size and zero is a finite velocity.
//
// Debug builds fill it with a recognizable quiet-NaN payload first and report
// any entry that still holds it. Release builds compile this away and keep the
// plain zero initialization. Comparison is by bit pattern because NaN never
// compares equal to itself. Diagnostics stay per sample, matching the
// surrounding messages; `evolve_ensemble` reports which sample failed.
#ifndef NDEBUG
inline constexpr std::uint64_t phase_space_velocity_poison_bits = 0x7ff8000051434b50ULL;

[[nodiscard]] double phase_space_velocity_poison() noexcept
{
    double value = 0.0;
    std::memcpy(&value, &phase_space_velocity_poison_bits, sizeof(value));
    return value;
}

[[nodiscard]] bool is_phase_space_velocity_poison(const double value) noexcept
{
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits == phase_space_velocity_poison_bits;
}

[[nodiscard]] std::string phase_space_velocity_unfilled_message(const double time, const std::size_t dimension, const std::vector<double>& coordinate)
{
    std::ostringstream output;
    output << std::setprecision(std::numeric_limits<double>::max_digits10) << "Phase-space velocity callback did not fill every velocity entry:\n"
           << "  time: " << time << '\n'
           << "  dimension: " << dimension << '\n'
           << "  coordinate: " << format_diagnostic_vector(coordinate);
    return output.str();
}
#endif

void poison_velocity_buffer([[maybe_unused]] std::vector<double>& velocities) noexcept
{
#ifndef NDEBUG
    std::fill(velocities.begin(), velocities.end(), phase_space_velocity_poison());
#endif
}

[[nodiscard]] EnsembleEvolutionResult invalid_result(std::string message, const std::size_t failed_sample_index = invalid_sample_index)
{
    EnsembleEvolutionResult result;
    result.status = EnsembleEvolutionStatus::invalid_input;
    result.failed_sample_index = failed_sample_index;
    result.rk23_status = qmomkit::AdaptiveRkStatus::invalid_input;
    result.message = std::move(message);
    return result;
}

[[nodiscard]] EnsembleEvolutionResult validate_inputs(const Ensemble& ensemble, const qmomkit::PhaseSpaceVelocity& velocity, const double t0,
                                                      const double t1)
{
    if (ensemble.empty())
    {
        return invalid_result(internal::ensemble_problem_message(internal::EnsembleProblem::empty));
    }
    if (!velocity)
    {
        return invalid_result("Phase-space velocity callback must be provided.");
    }
    if (!std::isfinite(t0) || !std::isfinite(t1) || !std::isfinite(t1 - t0))
    {
        return invalid_result("Evolution time input contains non-finite values.");
    }
    if (t0 == t1)
    {
        return invalid_result("Evolution time interval must be nonzero.");
    }

    // Every problem with the samples, a non-finite coordinate included, is
    // reported as invalid input for the sample it was found in.
    const internal::EnsembleCheck check = internal::check_ensemble(ensemble);
    if (check.problem != internal::EnsembleProblem::none)
    {
        return invalid_result(internal::ensemble_problem_message(check.problem), check.sample_index);
    }
    return {};
}

void checked_velocity(const qmomkit::PhaseSpaceVelocity& velocity, const std::size_t dimensions, const double t,
                      const std::vector<double>& coordinates, std::vector<double>& derivative)
{
    poison_velocity_buffer(derivative);
    velocity(t, coordinates, derivative, dimensions);
    if (derivative.size() != dimensions)
    {
        throw std::invalid_argument(phase_space_velocity_wrong_dimension_message(dimensions, derivative.size()));
    }
#ifndef NDEBUG
    // Ahead of the finiteness scan so that an unwritten entry is reported as
    // the missing write it is, rather than as a non-finite value the callback
    // never produced.
    for (std::size_t dimension = 0; dimension < derivative.size(); ++dimension)
    {
        if (is_phase_space_velocity_poison(derivative[dimension]))
        {
            throw std::invalid_argument(phase_space_velocity_unfilled_message(t, dimension, coordinates));
        }
    }
#endif
    for (std::size_t dimension = 0; dimension < derivative.size(); ++dimension)
    {
        if (!std::isfinite(derivative[dimension]))
        {
            throw std::invalid_argument(phase_space_velocity_non_finite_message(t, dimension, derivative[dimension], coordinates, derivative));
        }
    }
}

[[nodiscard]] EnsembleEvolutionResult evolve_samples(const Ensemble& ensemble, const std::size_t dimensions,
                                                     const qmomkit::PhaseSpaceVelocity& velocity, const double t0, const double t1,
                                                     const qmomkit::AdaptiveRkOptions& options)
{
    Ensemble evolved;
    evolved.reserve(ensemble.size());

    EnsembleEvolutionResult result;
    result.final_time = t0;

    for (std::size_t sample_index = 0; sample_index < ensemble.size(); ++sample_index)
    {
        qmomkit::AdaptiveRkResult rk_result;
        try
        {
            rk_result = qmomkit::integrate_embedded_rk(
                ensemble[sample_index], t0, t1,
                [&velocity, dimensions](const double rhs_time, const std::vector<double>& coordinates, std::vector<double>& derivative)
                { checked_velocity(velocity, dimensions, rhs_time, coordinates, derivative); }, qmomkit::EmbeddedRkMethod::rk23, options);
        }
        catch (const std::exception& error)
        {
            result.ensemble.clear();
            result.status = EnsembleEvolutionStatus::velocity_failed;
            result.failed_sample_index = sample_index;
            result.rk23_status = qmomkit::AdaptiveRkStatus::invalid_input;
            result.message = error.what();
            result.final_time = t0;
            return result;
        }

        // Evolution is all-or-nothing for the returned ensemble. Diagnostics
        // keep the failed index and RK23 counters, but partial sample states
        // are not exposed as a successful realization.
        result.total_accepted_steps += rk_result.accepted_steps;
        result.total_rejected_steps += rk_result.rejected_steps;
        result.total_function_evaluations += rk_result.function_evaluations;
        result.maximum_rejected_steps = std::max(result.maximum_rejected_steps, rk_result.rejected_steps);
        result.rk23_status = rk_result.status;

        if (rk_result.status != qmomkit::AdaptiveRkStatus::success)
        {
            result.ensemble.clear();
            result.status = EnsembleEvolutionStatus::rk23_failed;
            result.failed_sample_index = sample_index;
            result.message = rk_result.message;
            result.final_time = rk_result.final_time;
            return result;
        }
        if (rk_result.state.size() != dimensions || !all_finite(rk_result.state))
        {
            result.ensemble.clear();
            result.status = EnsembleEvolutionStatus::rk23_failed;
            result.failed_sample_index = sample_index;
            result.rk23_status = qmomkit::AdaptiveRkStatus::invalid_input;
            result.message = "RK23 returned an invalid sample state.";
            result.final_time = rk_result.final_time;
            return result;
        }

        evolved.push_back(std::move(rk_result.state));
    }

    result.ensemble = std::move(evolved);
    result.status = EnsembleEvolutionStatus::success;
    result.failed_sample_index = invalid_sample_index;
    result.rk23_status = qmomkit::AdaptiveRkStatus::success;
    result.message.clear();
    result.final_time = t1;
    return result;
}

} // namespace

EnsembleEvolutionResult evolve_ensemble(const Ensemble& ensemble, const qmomkit::PhaseSpaceVelocity& velocity, const double t0, const double t1,
                                        const qmomkit::AdaptiveRkOptions& options)
{
    EnsembleEvolutionResult validation = validate_inputs(ensemble, velocity, t0, t1);
    if (!validation.ok())
    {
        validation.final_time = t0;
        return validation;
    }

    return evolve_samples(ensemble, ensemble.front().size(), velocity, t0, t1, options);
}

} // namespace qmomkit::monte_carlo
