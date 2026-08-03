#include "internal/continuous_moment_rates.hpp"

#include "internal/cqmom_case.hpp"
#include "internal/moment_reconstruction.hpp"
#include "internal/numeric_utils.hpp"
#include "qmomkit/moment_rates.hpp"
#include "qmomkit/moment_set.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace qmomkit::internal
{
namespace
{

[[nodiscard]] std::vector<MultiIndex> qmom_moment_indices(const std::size_t moment_count)
{
    std::vector<MultiIndex> indices;
    indices.reserve(moment_count);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        indices.push_back(MultiIndex{order});
    }
    return indices;
}

[[nodiscard]] QmomContinuousMomentRateResult qmom_invalid_input(std::string message)
{
    return QmomContinuousMomentRateResult{{}, ContinuousMomentRateStatus::invalid_input, QmomStatus::invalid_input, std::move(message)};
}

[[nodiscard]] QmomContinuousMomentRateResult qmom_velocity_failed(std::string message)
{
    return QmomContinuousMomentRateResult{{}, ContinuousMomentRateStatus::velocity_failed, QmomStatus::success, std::move(message)};
}

[[nodiscard]] ContinuousMomentRateStatus rate_status_from_qmom_status(const QmomStatus status)
{
    return status == QmomStatus::invalid_input ? ContinuousMomentRateStatus::invalid_input : ContinuousMomentRateStatus::inversion_failed;
}

[[nodiscard]] CqmomContinuousMomentRateResult cqmom_invalid_input(std::string message)
{
    return CqmomContinuousMomentRateResult{
        std::nullopt, ContinuousMomentRateStatus::invalid_input, CqmomStatus::invalid_input, std::move(message), {}};
}

[[nodiscard]] CqmomContinuousMomentRateResult cqmom_velocity_failed(std::string message)
{
    return CqmomContinuousMomentRateResult{std::nullopt, ContinuousMomentRateStatus::velocity_failed, CqmomStatus::success, std::move(message), {}};
}

[[nodiscard]] ContinuousMomentRateStatus rate_status_from_cqmom_status(const CqmomStatus status)
{
    return status == CqmomStatus::invalid_input ? ContinuousMomentRateStatus::invalid_input : ContinuousMomentRateStatus::inversion_failed;
}

void validate_terminal_abscissas_for_velocity(const std::vector<std::vector<double>>& terminal_abscissas)
{
    if (terminal_abscissas.empty())
    {
        throw std::invalid_argument("Terminal velocity evaluation requires at least one phase-space dimension.");
    }

    const std::size_t node_count = terminal_abscissas.front().size();
    if (node_count == 0)
    {
        throw std::invalid_argument("Terminal velocity evaluation requires at least one terminal node.");
    }
    for (const std::vector<double>& dimension_abscissas : terminal_abscissas)
    {
        if (dimension_abscissas.size() != node_count)
        {
            throw std::invalid_argument("Every terminal abscissa row must contain one value per terminal node.");
        }
        if (!all_finite(dimension_abscissas))
        {
            throw std::invalid_argument("Terminal abscissas must be finite.");
        }
    }
}

void validate_qmom_rate_request(const std::vector<double>& moments, const double time, const PhaseSpaceVelocity& velocity)
{
    if (moments.empty())
    {
        throw std::invalid_argument("QMOM continuous moment-rate evaluation requires at least one moment.");
    }
    if (!std::isfinite(time))
    {
        throw std::invalid_argument("QMOM continuous moment-rate evaluation time must be finite.");
    }
    if (!velocity)
    {
        throw std::invalid_argument("Phase-space velocity callback must be provided.");
    }
}

void validate_cqmom_rate_request(const CqmomLeanCase& lean_case, const double time, const PhaseSpaceVelocity& velocity)
{
    if (!std::isfinite(time))
    {
        throw std::invalid_argument("CQMOM continuous moment-rate evaluation time must be finite.");
    }
    if (!velocity)
    {
        throw std::invalid_argument("Phase-space velocity callback must be provided.");
    }
    if (lean_case.moments.values.size() != lean_case.moments.layout.required_indices.size())
    {
        throw std::invalid_argument("CQMOM moment data values must match the layout size.");
    }
}

} // namespace

std::vector<std::vector<double>> evaluate_terminal_velocities(const PhaseSpaceVelocity& velocity, const double time,
                                                              const std::vector<std::vector<double>>& terminal_abscissas)
{
    if (!velocity)
    {
        throw std::invalid_argument("Phase-space velocity callback must be provided.");
    }
    if (!std::isfinite(time))
    {
        throw std::invalid_argument("Terminal velocity evaluation time must be finite.");
    }
    validate_terminal_abscissas_for_velocity(terminal_abscissas);

    const std::size_t dimensions = terminal_abscissas.size();
    const std::size_t node_count = terminal_abscissas.front().size();
    std::vector<std::vector<double>> velocities(node_count, std::vector<double>(dimensions, 0.0));
    std::vector<double> node_coordinate(dimensions, 0.0);

    for (std::size_t node = 0; node < node_count; ++node)
    {
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
        {
            node_coordinate[dimension] = terminal_abscissas[dimension][node];
        }

        std::vector<double> node_velocity;
        try
        {
            node_velocity = velocity(time, node_coordinate);
        }
        catch (const std::exception& error)
        {
            throw std::invalid_argument(error.what());
        }
        validate_phase_space_velocity_result(node_velocity, dimensions, time, node_coordinate);

        velocities[node] = std::move(node_velocity);
    }

    return velocities;
}

QmomContinuousMomentRateResult evaluate_qmom_continuous_moment_rates(const std::vector<double>& moments, const std::size_t node_count,
                                                                     const double time, const PhaseSpaceVelocity& velocity,
                                                                     const QmomOptions& qmom_options)
{
    try
    {
        validate_qmom_rate_request(moments, time, velocity);
    }
    catch (const std::exception& error)
    {
        return qmom_invalid_input(error.what());
    }

    const QmomResult quadrature = compute_qmom(MomentSet(moments), node_count, qmom_options);
    if (!quadrature.ok())
    {
        return QmomContinuousMomentRateResult{{}, rate_status_from_qmom_status(quadrature.status), quadrature.status, quadrature.message};
    }

    std::vector<std::vector<double>> velocities;
    try
    {
        velocities = evaluate_terminal_velocities(velocity, time, {quadrature.abscissas});
    }
    catch (const std::exception& error)
    {
        return qmom_velocity_failed(error.what());
    }

    return QmomContinuousMomentRateResult{
        reconstruct_continuous_moment_rates(quadrature.weights, {quadrature.abscissas}, velocities, qmom_moment_indices(moments.size())),
        ContinuousMomentRateStatus::success, quadrature.status, ""};
}

CqmomContinuousMomentRateResult evaluate_cqmom_continuous_moment_rates(const CqmomLeanCase& lean_case, const double time,
                                                                       const PhaseSpaceVelocity& velocity)
{
    try
    {
        validate_cqmom_rate_request(lean_case, time, velocity);
    }
    catch (const std::exception& error)
    {
        return cqmom_invalid_input(error.what());
    }

    const CqmomResult quadrature = compute_cqmom_internal(lean_case);
    if (!quadrature.ok())
    {
        return CqmomContinuousMomentRateResult{std::nullopt, rate_status_from_cqmom_status(quadrature.status), quadrature.status, quadrature.message,
                                               quadrature.diagnostics};
    }

    std::vector<std::vector<double>> velocities;
    try
    {
        velocities = evaluate_terminal_velocities(velocity, time, quadrature.quadrature.abscissas());
    }
    catch (const std::exception& error)
    {
        return cqmom_velocity_failed(error.what());
    }

    std::vector<double> rates = reconstruct_continuous_moment_rates(quadrature.quadrature.weights(), quadrature.quadrature.abscissas(), velocities,
                                                                    lean_case.moments.layout.required_indices);
    return CqmomContinuousMomentRateResult{CqmomMomentData{lean_case.moments.layout, std::move(rates)}, ContinuousMomentRateStatus::success,
                                           quadrature.status, "", quadrature.diagnostics};
}

} // namespace qmomkit::internal

namespace qmomkit
{
namespace
{

[[nodiscard]] ContinuousMomentRateStatus to_public_status(const internal::ContinuousMomentRateStatus status)
{
    switch (status)
    {
    case internal::ContinuousMomentRateStatus::success:
        return ContinuousMomentRateStatus::success;
    case internal::ContinuousMomentRateStatus::invalid_input:
        return ContinuousMomentRateStatus::invalid_input;
    case internal::ContinuousMomentRateStatus::inversion_failed:
        return ContinuousMomentRateStatus::inversion_failed;
    case internal::ContinuousMomentRateStatus::velocity_failed:
        return ContinuousMomentRateStatus::velocity_failed;
    }
    return ContinuousMomentRateStatus::invalid_input;
}

} // namespace

QmomContinuousMomentRateResult evaluate_qmom_continuous_moment_rates(const std::vector<double>& moments, const std::size_t node_count,
                                                                     const double time, const PhaseSpaceVelocity& velocity,
                                                                     const QmomOptions& options)
{
    internal::QmomContinuousMomentRateResult result = internal::evaluate_qmom_continuous_moment_rates(moments, node_count, time, velocity, options);
    return QmomContinuousMomentRateResult{std::move(result.rates), to_public_status(result.status), result.inversion_status,
                                          std::move(result.message)};
}

CqmomContinuousMomentRateResult evaluate_cqmom_continuous_moment_rates(const CqmomLeanCase& case_, const double time,
                                                                       const PhaseSpaceVelocity& velocity)
{
    internal::CqmomContinuousMomentRateResult result = internal::evaluate_cqmom_continuous_moment_rates(case_, time, velocity);
    return CqmomContinuousMomentRateResult{std::move(result.rates), to_public_status(result.status), result.inversion_status,
                                           std::move(result.message), std::move(result.diagnostics)};
}

} // namespace qmomkit
