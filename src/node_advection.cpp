#include "qmomkit/node_advection.hpp"

#include "internal/cqmom_case.hpp"
#include "internal/numeric_utils.hpp"
#include "internal/performance.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace qmomkit
{
namespace
{

[[nodiscard]] bool all_finite(const std::vector<std::vector<double>>& values)
{
    return std::all_of(values.begin(), values.end(), [](const std::vector<double>& row) { return internal::all_finite(row); });
}

[[nodiscard]] NodeAdvectionDiagnostics diagnostics_from_rk23(const AdaptiveRk23Result& result)
{
    NodeAdvectionDiagnostics diagnostics;
    diagnostics.scheme = NodeAdvectionScheme::adaptive_rk23;
    diagnostics.step_count = result.accepted_steps;
    diagnostics.rejected_step_count = result.rejected_steps;
    diagnostics.rhs_evaluations = result.function_evaluations;
    diagnostics.final_time = result.t;
    diagnostics.last_step = result.last_step;
    diagnostics.adaptive_rk23_status = result.status;
    diagnostics.message = result.message;
    return diagnostics;
}

[[nodiscard]] std::vector<double> pack_nodes(const std::vector<std::vector<double>>& abscissas)
{
    if (abscissas.empty())
    {
        return {};
    }

    const std::size_t dimensions = abscissas.size();
    const std::size_t nodes = abscissas.front().size();
    // Continuous transport integrates node coordinates only, packed node-major
    // as [node0_dim0, node0_dim1, ...]; weights are not part of the RK state.
    std::vector<double> state(nodes * dimensions, 0.0);
    for (std::size_t node = 0; node < nodes; ++node)
    {
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
        {
            state[node * dimensions + dimension] = abscissas[dimension][node];
        }
    }
    return state;
}

[[nodiscard]] std::vector<std::vector<double>> unpack_nodes(const std::vector<double>& state, const std::size_t dimensions)
{
    const std::size_t nodes = state.size() / dimensions;
    std::vector<std::vector<double>> abscissas(dimensions, std::vector<double>(nodes, 0.0));
    for (std::size_t node = 0; node < nodes; ++node)
    {
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
        {
            abscissas[dimension][node] = state[node * dimensions + dimension];
        }
    }
    return abscissas;
}

void validate_weights(const std::vector<double>& weights)
{
    if (weights.empty())
    {
        throw std::invalid_argument("Quadrature must contain at least one node.");
    }
    if (!internal::all_finite(weights))
    {
        throw std::invalid_argument("Quadrature weights must be finite.");
    }
}

void validate_abscissas(const std::vector<std::vector<double>>& abscissas, const std::size_t dimensions, const std::size_t nodes)
{
    if (dimensions == 0)
    {
        throw std::invalid_argument("Phase-space dimension count must be positive.");
    }
    if (abscissas.size() != dimensions)
    {
        throw std::invalid_argument("Quadrature abscissas must contain one row per phase-space dimension.");
    }
    for (const std::vector<double>& dimension_abscissas : abscissas)
    {
        if (dimension_abscissas.size() != nodes)
        {
            throw std::invalid_argument("Every abscissa dimension must contain one value per quadrature node.");
        }
    }
    if (!all_finite(abscissas))
    {
        throw std::invalid_argument("Quadrature abscissas must be finite.");
    }
}

[[nodiscard]] std::vector<double> node_rhs_from_velocity(const PhaseSpaceVelocity& velocity, const std::size_t dimensions, const double t,
                                                         const std::vector<double>& state)
{
    const std::size_t nodes = state.size() / dimensions;
    std::vector<double> derivative(state.size(), 0.0);
    std::vector<double> node_coordinate(dimensions, 0.0);

    // The velocity callback is evaluated independently at each quadrature node
    // and returns derivatives in phase-space dimension order.
    for (std::size_t node = 0; node < nodes; ++node)
    {
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
        {
            node_coordinate[dimension] = state[node * dimensions + dimension];
        }

        const std::vector<double> node_velocity = velocity(t, node_coordinate);
        internal::performance::record_node_rhs_evaluation();
        if (node_velocity.size() != dimensions)
        {
            throw std::invalid_argument("Phase-space velocity callback returned the wrong dimension.");
        }
        if (!internal::all_finite(node_velocity))
        {
            throw std::invalid_argument("Phase-space velocity callback returned a non-finite value.");
        }

        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
        {
            derivative[node * dimensions + dimension] = node_velocity[dimension];
        }
    }

    return derivative;
}

struct IntegratedNodes
{
    std::vector<std::vector<double>> abscissas;
    NodeAdvectionStatus status = NodeAdvectionStatus::success;
    NodeAdvectionDiagnostics diagnostics;
    std::string message;
};

[[nodiscard]] std::vector<double> add_scaled(const std::vector<double>& state, const double scale, const std::vector<double>& derivative)
{
    if (state.size() != derivative.size())
    {
        throw std::invalid_argument("Vector sizes do not match.");
    }
    std::vector<double> result(state.size(), 0.0);
    for (std::size_t i = 0; i < state.size(); ++i)
    {
        result[i] = state[i] + scale * derivative[i];
    }
    return result;
}

[[nodiscard]] IntegratedNodes integrate_nodes_rk23(const std::vector<std::vector<double>>& abscissas, const double t, const double time_step,
                                                   const PhaseSpaceVelocity& velocity, const NodeAdvectionOptions& options,
                                                   const std::size_t dimensions, const std::vector<double>& initial_state)
{
    AdaptiveRk23Result rk_result;
    try
    {
        rk_result = integrate_rk23(
            initial_state, t, t + time_step, [&velocity, dimensions](const double rhs_time, const std::vector<double>& state)
            { return node_rhs_from_velocity(velocity, dimensions, rhs_time, state); }, options.rk23_options);
    }
    catch (const std::exception& error)
    {
        NodeAdvectionDiagnostics diagnostics;
        diagnostics.scheme = NodeAdvectionScheme::adaptive_rk23;
        diagnostics.adaptive_rk23_status = AdaptiveRk23Status::invalid_input;
        diagnostics.final_time = t;
        diagnostics.message = error.what();
        return IntegratedNodes{abscissas, NodeAdvectionStatus::velocity_failed, diagnostics, error.what()};
    }

    NodeAdvectionDiagnostics diagnostics = diagnostics_from_rk23(rk_result);
    if (rk_result.status != AdaptiveRk23Status::success)
    {
        return IntegratedNodes{abscissas, NodeAdvectionStatus::integration_failed, diagnostics, rk_result.message};
    }
    if (rk_result.state.size() != initial_state.size() || !internal::all_finite(rk_result.state))
    {
        diagnostics.message = "RK23 returned an invalid node state.";
        return IntegratedNodes{abscissas, NodeAdvectionStatus::integration_failed, diagnostics, diagnostics.message};
    }

    return IntegratedNodes{unpack_nodes(rk_result.state, dimensions), NodeAdvectionStatus::success, diagnostics, ""};
}

[[nodiscard]] IntegratedNodes integrate_nodes_ssp_rk3(const std::vector<std::vector<double>>& abscissas, const double t, const double time_step,
                                                      const PhaseSpaceVelocity& velocity, const std::size_t dimensions,
                                                      const std::vector<double>& initial_state)
{
    NodeAdvectionDiagnostics diagnostics;
    diagnostics.scheme = NodeAdvectionScheme::ssp_rk3;
    diagnostics.final_time = t;

    if (time_step == 0.0)
    {
        return IntegratedNodes{abscissas, NodeAdvectionStatus::success, diagnostics, ""};
    }

    try
    {
        // This is the SSP-RK3 trajectory formula only. It deliberately avoids
        // the source-update wrapper, including source predicates, retries, and
        // source-specific performance counters.
        std::size_t rhs_evaluations = 0;
        const auto rhs = [&velocity, dimensions, &rhs_evaluations](const double rhs_time, const std::vector<double>& state)
        {
            ++rhs_evaluations;
            return node_rhs_from_velocity(velocity, dimensions, rhs_time, state);
        };

        const std::vector<double> k1 = rhs(t, initial_state);
        const std::vector<double> stage_1 = add_scaled(initial_state, time_step, k1);

        const std::vector<double> k2 = rhs(t + time_step, stage_1);
        std::vector<double> stage_2(initial_state.size(), 0.0);
        for (std::size_t i = 0; i < initial_state.size(); ++i)
        {
            stage_2[i] = 0.75 * initial_state[i] + 0.25 * (stage_1[i] + time_step * k2[i]);
        }

        const std::vector<double> k3 = rhs(t + 0.5 * time_step, stage_2);
        std::vector<double> updated_state(initial_state.size(), 0.0);
        for (std::size_t i = 0; i < initial_state.size(); ++i)
        {
            updated_state[i] = (initial_state[i] + 2.0 * (stage_2[i] + time_step * k3[i])) / 3.0;
        }

        diagnostics.step_count = 1;
        diagnostics.rhs_evaluations = rhs_evaluations;
        diagnostics.final_time = t + time_step;
        diagnostics.last_step = time_step;
        if (!internal::all_finite(updated_state))
        {
            diagnostics.message = "SSP-RK3 returned an invalid node state.";
            return IntegratedNodes{abscissas, NodeAdvectionStatus::integration_failed, diagnostics, diagnostics.message};
        }

        return IntegratedNodes{unpack_nodes(updated_state, dimensions), NodeAdvectionStatus::success, diagnostics, ""};
    }
    catch (const std::exception& error)
    {
        diagnostics.message = error.what();
        return IntegratedNodes{abscissas, NodeAdvectionStatus::velocity_failed, diagnostics, error.what()};
    }
}

[[nodiscard]] IntegratedNodes integrate_nodes(const std::vector<std::vector<double>>& abscissas, const double t, const double time_step,
                                              const PhaseSpaceVelocity& velocity, const NodeAdvectionOptions& options)
{
    if (!velocity)
    {
        return IntegratedNodes{abscissas, NodeAdvectionStatus::invalid_input, {}, "Phase-space velocity callback must be provided."};
    }
    if (!std::isfinite(t) || !std::isfinite(time_step) || !std::isfinite(t + time_step))
    {
        return IntegratedNodes{abscissas, NodeAdvectionStatus::invalid_input, {}, "Node advection time input contains non-finite values."};
    }

    const std::size_t dimensions = abscissas.size();
    const std::vector<double> initial_state = pack_nodes(abscissas);
    switch (options.scheme)
    {
    case NodeAdvectionScheme::ssp_rk3:
        return integrate_nodes_ssp_rk3(abscissas, t, time_step, velocity, dimensions, initial_state);
    case NodeAdvectionScheme::adaptive_rk23:
        return integrate_nodes_rk23(abscissas, t, time_step, velocity, options, dimensions, initial_state);
    }

    NodeAdvectionDiagnostics diagnostics;
    diagnostics.message = "Unsupported node-advection scheme.";
    return IntegratedNodes{abscissas, NodeAdvectionStatus::invalid_input, diagnostics, diagnostics.message};
}

[[nodiscard]] const char* scheme_name(const NodeAdvectionScheme scheme)
{
    switch (scheme)
    {
    case NodeAdvectionScheme::ssp_rk3:
        return "ssp_rk3";
    case NodeAdvectionScheme::adaptive_rk23:
        return "adaptive_rk23";
    }
    return "unknown";
}

} // namespace

QmomNodeAdvectionResult advance_nodes(const Quadrature1D& quadrature, const double t, const double time_step, const PhaseSpaceVelocity& velocity,
                                      const NodeAdvectionOptions& options)
{
    try
    {
        if (!quadrature.ok())
        {
            return QmomNodeAdvectionResult{quadrature, NodeAdvectionStatus::invalid_input, {}, "Input QMOM quadrature is not successful."};
        }
        validate_weights(quadrature.weights);
        validate_abscissas({quadrature.abscissas}, 1, quadrature.weights.size());
    }
    catch (const std::exception& error)
    {
        return QmomNodeAdvectionResult{quadrature, NodeAdvectionStatus::invalid_input, {}, error.what()};
    }

    const IntegratedNodes integrated = integrate_nodes({quadrature.abscissas}, t, time_step, velocity, options);
    if (integrated.status != NodeAdvectionStatus::success)
    {
        return QmomNodeAdvectionResult{quadrature, integrated.status, integrated.diagnostics, integrated.message};
    }

    Quadrature1D updated = quadrature;
    updated.abscissas = integrated.abscissas.front();
    return QmomNodeAdvectionResult{std::move(updated), NodeAdvectionStatus::success, integrated.diagnostics, ""};
}

QmomNodeAdvectionResult advance_nodes(const QmomResult& quadrature, const double t, const double time_step, const PhaseSpaceVelocity& velocity,
                                      const NodeAdvectionOptions& options)
{
    if (!quadrature.ok())
    {
        return QmomNodeAdvectionResult{{quadrature.abscissas, quadrature.weights, QuadratureStatus::invalid_input, quadrature.message},
                                       NodeAdvectionStatus::invalid_input,
                                       {},
                                       "Input QMOM result is not successful."};
    }

    return advance_nodes(Quadrature1D{quadrature.abscissas, quadrature.weights, QuadratureStatus::success, quadrature.message}, t, time_step,
                         velocity, options);
}

CqmomNodeAdvectionResult advance_nodes(const CqmomQuadrature& quadrature, const double t, const double time_step, const PhaseSpaceVelocity& velocity,
                                       const NodeAdvectionOptions& options)
{
    try
    {
        validate_weights(quadrature.weights());
        validate_abscissas(quadrature.abscissas(), quadrature.number_of_dimensions(), quadrature.weights().size());
        if (!quadrature.has_rectangular_storage() && quadrature.terminal_paths().size() != quadrature.weights().size())
        {
            throw std::invalid_argument("Adaptive CQMOM quadrature must provide one terminal path per terminal node.");
        }
    }
    catch (const std::exception& error)
    {
        return CqmomNodeAdvectionResult{internal::create_cqmom_quadrature(quadrature.hierarchy(), quadrature.weights(), quadrature.abscissas(),
                                                                          quadrature.adaptive_branch_counts_enabled(),
                                                                          quadrature.has_rectangular_storage(), quadrature.terminal_paths()),
                                        NodeAdvectionStatus::invalid_input,
                                        {},
                                        error.what()};
    }

    const IntegratedNodes integrated = integrate_nodes(quadrature.abscissas(), t, time_step, velocity, options);
    if (integrated.status != NodeAdvectionStatus::success)
    {
        return CqmomNodeAdvectionResult{internal::create_cqmom_quadrature(quadrature.hierarchy(), quadrature.weights(), quadrature.abscissas(),
                                                                          quadrature.adaptive_branch_counts_enabled(),
                                                                          quadrature.has_rectangular_storage(), quadrature.terminal_paths()),
                                        integrated.status, integrated.diagnostics, integrated.message};
    }

    return CqmomNodeAdvectionResult{internal::create_cqmom_quadrature(quadrature.hierarchy(), quadrature.weights(), integrated.abscissas,
                                                                      quadrature.adaptive_branch_counts_enabled(),
                                                                      quadrature.has_rectangular_storage(), quadrature.terminal_paths()),
                                    NodeAdvectionStatus::success, integrated.diagnostics, ""};
}

std::string format_node_advection_options(const NodeAdvectionOptions& options)
{
    std::ostringstream output;
    output << std::setprecision(17);
    output << "Node advection settings:\n";
    output << "  scheme = " << scheme_name(options.scheme) << '\n';

    if (options.scheme != NodeAdvectionScheme::adaptive_rk23)
    {
        return output.str();
    }

    output << "  rk23.relative_tolerance = " << options.rk23_options.relative_tolerance << '\n';
    output << "  rk23.absolute_tolerance = " << options.rk23_options.absolute_tolerance << '\n';
    output << "  rk23.initial_step = ";
    if (options.rk23_options.initial_step.has_value())
    {
        output << *options.rk23_options.initial_step;
    }
    else
    {
        output << "automatic";
    }
    output << '\n';
    output << "  rk23.max_step = ";
    if (options.rk23_options.max_step > 0.0)
    {
        output << options.rk23_options.max_step;
    }
    else
    {
        output << "full integration interval";
    }
    output << '\n';
    output << "  rk23.min_step = ";
    if (options.rk23_options.min_step > 0.0)
    {
        output << options.rk23_options.min_step;
    }
    else
    {
        output << "machine-dependent automatic minimum";
    }
    output << '\n';
    output << "  rk23.max_steps = " << options.rk23_options.max_steps << '\n';
    return output.str();
}

} // namespace qmomkit
