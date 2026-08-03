#include "monte_carlo/evolution.hpp"

#include "explicit_distribution_detail.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
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
        return invalid_result("Ensemble must contain at least one sample.");
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

    const std::size_t dimensions = ensemble.front().size();
    if (dimensions == 0)
    {
        return invalid_result("Samples must contain at least one phase-space coordinate.", 0);
    }
    for (std::size_t sample_index = 0; sample_index < ensemble.size(); ++sample_index)
    {
        if (ensemble[sample_index].size() != dimensions)
        {
            return invalid_result("All samples must have the same phase-space dimension.", sample_index);
        }
        if (!all_finite(ensemble[sample_index]))
        {
            return invalid_result("Sample coordinate is non-finite.", sample_index);
        }
    }

    return {};
}

struct StateEvolutionResult
{
    Ensemble states;
    EnsembleEvolutionStatus status = EnsembleEvolutionStatus::success;
    std::size_t failed_state_index = invalid_sample_index;
    qmomkit::AdaptiveRkStatus rk23_status = qmomkit::AdaptiveRkStatus::success;
    std::string message;
    std::size_t total_accepted_steps = 0;
    std::size_t total_rejected_steps = 0;
    std::size_t total_function_evaluations = 0;
    std::size_t maximum_rejected_steps = 0;
    double final_time = 0.0;
};

[[nodiscard]] std::vector<double> checked_velocity(const qmomkit::PhaseSpaceVelocity& velocity, const std::size_t dimensions, const double t,
                                                   const std::vector<double>& coordinates)
{
    const std::vector<double> derivative = velocity(t, coordinates);
    if (derivative.size() != dimensions)
    {
        throw std::invalid_argument(phase_space_velocity_wrong_dimension_message(dimensions, derivative.size()));
    }
    for (std::size_t dimension = 0; dimension < derivative.size(); ++dimension)
    {
        if (!std::isfinite(derivative[dimension]))
        {
            throw std::invalid_argument(phase_space_velocity_non_finite_message(t, dimension, derivative[dimension], coordinates, derivative));
        }
    }
    return derivative;
}

[[nodiscard]] StateEvolutionResult evolve_states(const Ensemble& states, const std::size_t dimensions, const qmomkit::PhaseSpaceVelocity& velocity,
                                                 const double t0, const double t1, const qmomkit::AdaptiveRkOptions& options,
                                                 std::string invalid_state_message)
{
    Ensemble evolved;
    evolved.reserve(states.size());

    StateEvolutionResult result;
    result.final_time = t0;

    for (std::size_t state_index = 0; state_index < states.size(); ++state_index)
    {
        qmomkit::AdaptiveRkResult rk_result;
        try
        {
            rk_result = qmomkit::integrate_rk23(
                states[state_index], t0, t1, [&velocity, dimensions](const double rhs_time, const std::vector<double>& coordinates)
                { return checked_velocity(velocity, dimensions, rhs_time, coordinates); }, options);
        }
        catch (const std::exception& error)
        {
            result.states.clear();
            result.status = EnsembleEvolutionStatus::velocity_failed;
            result.failed_state_index = state_index;
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
            result.states.clear();
            result.status = EnsembleEvolutionStatus::rk23_failed;
            result.failed_state_index = state_index;
            result.message = rk_result.message;
            result.final_time = rk_result.t;
            return result;
        }
        if (rk_result.state.size() != dimensions || !all_finite(rk_result.state))
        {
            result.states.clear();
            result.status = EnsembleEvolutionStatus::rk23_failed;
            result.failed_state_index = state_index;
            result.rk23_status = qmomkit::AdaptiveRkStatus::invalid_input;
            result.message = std::move(invalid_state_message);
            result.final_time = rk_result.t;
            return result;
        }

        evolved.push_back(std::move(rk_result.state));
    }

    result.states = std::move(evolved);
    result.status = EnsembleEvolutionStatus::success;
    result.failed_state_index = invalid_sample_index;
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

    const std::size_t dimensions = ensemble.front().size();
    StateEvolutionResult states = evolve_states(ensemble, dimensions, velocity, t0, t1, options, "RK23 returned an invalid sample state.");

    EnsembleEvolutionResult result;
    result.ensemble = std::move(states.states);
    result.status = states.status;
    result.failed_sample_index = states.failed_state_index;
    result.rk23_status = states.rk23_status;
    result.message = states.message;
    result.total_accepted_steps = states.total_accepted_steps;
    result.total_rejected_steps = states.total_rejected_steps;
    result.total_function_evaluations = states.total_function_evaluations;
    result.maximum_rejected_steps = states.maximum_rejected_steps;
    result.final_time = states.final_time;
    return result;
}

qmomkit::QmomNodeAdvectionResult evolve_quadrature_nodes(const qmomkit::Quadrature1D& quadrature, const qmomkit::PhaseSpaceVelocity& velocity,
                                                         const double t0, const double t1, const qmomkit::AdaptiveRkOptions& options)
{
    qmomkit::NodeAdvectionOptions node_options;
    node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    node_options.adaptive_options = options;
    return qmomkit::advance_nodes(quadrature, t0, t1 - t0, velocity, node_options);
}

ExplicitQuadratureEvolutionResult evolve_explicit_quadrature(const ExplicitDiscreteDistribution& distribution,
                                                             const qmomkit::PhaseSpaceVelocity& velocity, const double t0, const double t1,
                                                             const qmomkit::AdaptiveRkOptions& options)
{
    ExplicitQuadratureEvolutionResult result;
    result.initial_distribution = distribution;
    result.final_time = t0;

    const detail::ExplicitDistributionValidation validation = detail::validate_explicit_distribution(distribution);
    if (!validation.ok)
    {
        result.status = ComparisonStatus::invalid_input;
        result.rk23_status = qmomkit::AdaptiveRkStatus::invalid_input;
        result.message = validation.message;
        return result;
    }
    if (!velocity)
    {
        result.status = ComparisonStatus::invalid_input;
        result.rk23_status = qmomkit::AdaptiveRkStatus::invalid_input;
        result.message = "Phase-space velocity callback must be provided.";
        return result;
    }
    if (!std::isfinite(t0) || !std::isfinite(t1) || t0 == t1)
    {
        result.status = ComparisonStatus::invalid_input;
        result.rk23_status = qmomkit::AdaptiveRkStatus::invalid_input;
        result.message = "Evolution time interval is invalid.";
        return result;
    }

    StateEvolutionResult states =
        evolve_states(distribution.nodes, validation.dimensions, velocity, t0, t1, options, "RK23 returned an invalid explicit node state.");
    result.total_accepted_steps = states.total_accepted_steps;
    result.total_rejected_steps = states.total_rejected_steps;
    result.total_function_evaluations = states.total_function_evaluations;
    result.maximum_rejected_steps = states.maximum_rejected_steps;
    result.rk23_status = states.rk23_status;
    result.final_time = states.final_time;

    if (states.status != EnsembleEvolutionStatus::success)
    {
        result.status = ComparisonStatus::qmomkit_failed;
        result.failed_node_index = states.failed_state_index;
        result.message = states.message;
        result.evolved_distribution = {};
        return result;
    }

    result.evolved_distribution = ExplicitDiscreteDistribution{std::move(states.states), distribution.weights, validation.number_density};
    result.status = ComparisonStatus::success;
    result.failed_node_index = static_cast<std::size_t>(-1);
    result.rk23_status = qmomkit::AdaptiveRkStatus::success;
    result.message.clear();
    result.final_time = t1;
    return result;
}

} // namespace qmomkit::monte_carlo
