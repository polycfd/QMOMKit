#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <qmomkit/moment_rates.hpp>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace moment_rate_element_operator
{

struct Element
{
    std::vector<double> nodes;
    std::vector<double> quadrature_points;
    std::vector<double> quadrature_weights;
    std::vector<std::vector<double>> basis_values;
    std::vector<std::vector<double>> inverse_mass;

    [[nodiscard]] std::size_t basis_size() const noexcept
    {
        return nodes.size();
    }

    [[nodiscard]] std::size_t quadrature_count() const noexcept
    {
        return quadrature_points.size();
    }
};

struct ElementRhsResult
{
    bool ok = false;
    std::vector<double> coefficient_rates;
    std::size_t failed_quadrature_point = std::numeric_limits<std::size_t>::max();
    std::string status;
    std::string inversion_status;
    std::string message;
};

[[nodiscard]] inline std::string to_string(const qmomkit::ContinuousMomentRateStatus status)
{
    switch (status)
    {
    case qmomkit::ContinuousMomentRateStatus::success:
        return "success";
    case qmomkit::ContinuousMomentRateStatus::invalid_input:
        return "invalid_input";
    case qmomkit::ContinuousMomentRateStatus::inversion_failed:
        return "inversion_failed";
    case qmomkit::ContinuousMomentRateStatus::velocity_failed:
        return "velocity_failed";
    }
    return "unknown";
}

[[nodiscard]] inline std::string to_string(const qmomkit::DiscontinuousMomentRateStatus status)
{
    switch (status)
    {
    case qmomkit::DiscontinuousMomentRateStatus::success:
        return "success";
    case qmomkit::DiscontinuousMomentRateStatus::invalid_input:
        return "invalid_input";
    case qmomkit::DiscontinuousMomentRateStatus::source_failed:
        return "source_failed";
    case qmomkit::DiscontinuousMomentRateStatus::invalid_rate:
        return "invalid_rate";
    }
    return "unknown";
}

[[nodiscard]] inline std::string to_string(const qmomkit::QmomStatus status)
{
    switch (status)
    {
    case qmomkit::QmomStatus::success:
        return "success";
    case qmomkit::QmomStatus::invalid_input:
        return "invalid_input";
    case qmomkit::QmomStatus::not_realizable:
        return "not_realizable";
    }
    return "unknown";
}

[[nodiscard]] inline std::string to_string(const qmomkit::CqmomStatus status)
{
    switch (status)
    {
    case qmomkit::CqmomStatus::success:
        return "success";
    case qmomkit::CqmomStatus::invalid_input:
        return "invalid_input";
    case qmomkit::CqmomStatus::qmom_failed:
        return "qmom_failed";
    case qmomkit::CqmomStatus::conditional_moment_failed:
        return "conditional_moment_failed";
    }
    return "unknown";
}

[[nodiscard]] inline double integer_power(const double base, const std::size_t exponent)
{
    double value = 1.0;
    for (std::size_t i = 0; i < exponent; ++i)
    {
        value *= base;
    }
    return value;
}

[[nodiscard]] inline std::vector<double> qmom_moments_from_terminals(const std::vector<double>& weights, const std::vector<double>& abscissas,
                                                                     const std::size_t moment_count)
{
    std::vector<double> moments(moment_count, 0.0);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            moments[order] += weights[node] * integer_power(abscissas[node], order);
        }
    }
    return moments;
}

[[nodiscard]] inline std::vector<double> selected_terminal_moments(const std::vector<double>& weights,
                                                                   const std::vector<std::vector<double>>& abscissas,
                                                                   const std::vector<qmomkit::MultiIndex>& indices)
{
    std::vector<double> moments;
    moments.reserve(indices.size());
    for (const qmomkit::MultiIndex& index : indices)
    {
        double value = 0.0;
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            double contribution = weights[node];
            for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
            {
                contribution *= integer_power(abscissas[dimension][node], index[dimension]);
            }
            value += contribution;
        }
        moments.push_back(value);
    }
    return moments;
}

[[nodiscard]] inline qmomkit::CqmomMomentData cqmom_moment_data_from_terminals(const qmomkit::CqmomMomentLayout& layout,
                                                                               const std::vector<double>& weights,
                                                                               const std::vector<std::vector<double>>& abscissas)
{
    return qmomkit::CqmomMomentData{layout, selected_terminal_moments(weights, abscissas, layout.required_indices)};
}

[[nodiscard]] inline std::vector<double> evaluate_lagrange_basis(const std::vector<double>& nodes, const double x)
{
    std::vector<double> values(nodes.size(), 1.0);
    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
        for (std::size_t j = 0; j < nodes.size(); ++j)
        {
            if (i != j)
            {
                values[i] *= (x - nodes[j]) / (nodes[i] - nodes[j]);
            }
        }
    }
    return values;
}

[[nodiscard]] inline std::vector<std::vector<double>> invert_matrix(const std::vector<std::vector<double>>& matrix)
{
    const std::size_t size = matrix.size();
    std::vector<std::vector<double>> augmented(size, std::vector<double>(2 * size, 0.0));
    for (std::size_t row = 0; row < size; ++row)
    {
        if (matrix[row].size() != size)
        {
            throw std::invalid_argument("Mass matrix must be square.");
        }
        for (std::size_t col = 0; col < size; ++col)
        {
            augmented[row][col] = matrix[row][col];
        }
        augmented[row][size + row] = 1.0;
    }

    for (std::size_t pivot = 0; pivot < size; ++pivot)
    {
        std::size_t best = pivot;
        for (std::size_t row = pivot + 1; row < size; ++row)
        {
            if (std::abs(augmented[row][pivot]) > std::abs(augmented[best][pivot]))
            {
                best = row;
            }
        }
        if (std::abs(augmented[best][pivot]) <= 1.0e-14)
        {
            throw std::runtime_error("Mass matrix is singular.");
        }
        if (best != pivot)
        {
            std::swap(augmented[best], augmented[pivot]);
        }

        const double pivot_value = augmented[pivot][pivot];
        for (double& value : augmented[pivot])
        {
            value /= pivot_value;
        }
        for (std::size_t row = 0; row < size; ++row)
        {
            if (row == pivot)
            {
                continue;
            }
            const double factor = augmented[row][pivot];
            for (std::size_t col = 0; col < 2 * size; ++col)
            {
                augmented[row][col] -= factor * augmented[pivot][col];
            }
        }
    }

    std::vector<std::vector<double>> inverse(size, std::vector<double>(size, 0.0));
    for (std::size_t row = 0; row < size; ++row)
    {
        for (std::size_t col = 0; col < size; ++col)
        {
            inverse[row][col] = augmented[row][size + col];
        }
    }
    return inverse;
}

[[nodiscard]] inline Element make_lagrange_element(const std::size_t degree)
{
    Element element;
    if (degree == 2)
    {
        element.nodes = {-1.0, 0.0, 1.0};
        const double a = std::sqrt(3.0 / 5.0);
        element.quadrature_points = {-a, 0.0, a};
        element.quadrature_weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    }
    else if (degree == 3)
    {
        element.nodes = {-1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0};
        const double a = std::sqrt((3.0 - 2.0 * std::sqrt(6.0 / 5.0)) / 7.0);
        const double b = std::sqrt((3.0 + 2.0 * std::sqrt(6.0 / 5.0)) / 7.0);
        element.quadrature_points = {-b, -a, a, b};
        element.quadrature_weights = {(18.0 - std::sqrt(30.0)) / 36.0, (18.0 + std::sqrt(30.0)) / 36.0, (18.0 + std::sqrt(30.0)) / 36.0,
                                      (18.0 - std::sqrt(30.0)) / 36.0};
    }
    else
    {
        throw std::invalid_argument("Prototype supports polynomial degree 2 or 3.");
    }

    element.basis_values.reserve(element.quadrature_points.size());
    for (const double x : element.quadrature_points)
    {
        element.basis_values.push_back(evaluate_lagrange_basis(element.nodes, x));
    }

    std::vector<std::vector<double>> mass(element.basis_size(), std::vector<double>(element.basis_size(), 0.0));
    for (std::size_t q = 0; q < element.quadrature_count(); ++q)
    {
        for (std::size_t i = 0; i < element.basis_size(); ++i)
        {
            for (std::size_t j = 0; j < element.basis_size(); ++j)
            {
                mass[i][j] += element.quadrature_weights[q] * element.basis_values[q][i] * element.basis_values[q][j];
            }
        }
    }
    element.inverse_mass = invert_matrix(mass);
    return element;
}

[[nodiscard]] inline std::size_t coefficient_index(const std::size_t basis, const std::size_t moment, const std::size_t moment_count)
{
    return basis * moment_count + moment;
}

[[nodiscard]] inline std::vector<double> constant_coefficients(const Element& element, const std::vector<double>& moments)
{
    std::vector<double> coefficients(element.basis_size() * moments.size(), 0.0);
    for (std::size_t basis = 0; basis < element.basis_size(); ++basis)
    {
        for (std::size_t moment = 0; moment < moments.size(); ++moment)
        {
            coefficients[coefficient_index(basis, moment, moments.size())] = moments[moment];
        }
    }
    return coefficients;
}

[[nodiscard]] inline std::vector<double> interpolate_element_state(const Element& element, const std::vector<double>& coefficients,
                                                                   const std::size_t moment_count)
{
    if (coefficients.size() != element.basis_size() * moment_count)
    {
        throw std::invalid_argument("Coefficient vector size does not match element shape.");
    }

    std::vector<double> local(element.quadrature_count() * moment_count, 0.0);
    for (std::size_t q = 0; q < element.quadrature_count(); ++q)
    {
        for (std::size_t basis = 0; basis < element.basis_size(); ++basis)
        {
            const double basis_value = element.basis_values[q][basis];
            for (std::size_t moment = 0; moment < moment_count; ++moment)
            {
                local[q * moment_count + moment] += basis_value * coefficients[coefficient_index(basis, moment, moment_count)];
            }
        }
    }
    return local;
}

[[nodiscard]] inline std::vector<double> project_rates_to_coefficients(const Element& element, const std::vector<double>& local_rates,
                                                                       const std::size_t moment_count)
{
    if (local_rates.size() != element.quadrature_count() * moment_count)
    {
        throw std::invalid_argument("Local rate vector size does not match element quadrature shape.");
    }

    std::vector<double> projected(element.basis_size() * moment_count, 0.0);
    for (std::size_t q = 0; q < element.quadrature_count(); ++q)
    {
        for (std::size_t basis = 0; basis < element.basis_size(); ++basis)
        {
            const double weighted_basis = element.quadrature_weights[q] * element.basis_values[q][basis];
            for (std::size_t moment = 0; moment < moment_count; ++moment)
            {
                projected[coefficient_index(basis, moment, moment_count)] += weighted_basis * local_rates[q * moment_count + moment];
            }
        }
    }

    std::vector<double> coefficient_rates(element.basis_size() * moment_count, 0.0);
    for (std::size_t moment = 0; moment < moment_count; ++moment)
    {
        for (std::size_t row = 0; row < element.basis_size(); ++row)
        {
            for (std::size_t col = 0; col < element.basis_size(); ++col)
            {
                coefficient_rates[coefficient_index(row, moment, moment_count)] +=
                    element.inverse_mass[row][col] * projected[coefficient_index(col, moment, moment_count)];
            }
        }
    }
    return coefficient_rates;
}

[[nodiscard]] inline ElementRhsResult success_result(std::vector<double> coefficient_rates)
{
    return ElementRhsResult{true, std::move(coefficient_rates), std::numeric_limits<std::size_t>::max(), "success", "success", ""};
}

[[nodiscard]] inline ElementRhsResult continuous_failure(const std::size_t point, const qmomkit::ContinuousMomentRateStatus status,
                                                         const std::string& inversion_status, std::string message)
{
    return ElementRhsResult{false, {}, point, to_string(status), inversion_status, std::move(message)};
}

[[nodiscard]] inline ElementRhsResult discontinuous_failure(const std::size_t point, const qmomkit::DiscontinuousMomentRateStatus status,
                                                            std::string message)
{
    return ElementRhsResult{false, {}, point, to_string(status), "not_applicable", std::move(message)};
}

[[nodiscard]] inline ElementRhsResult evaluate_qmom_continuous_rhs(const Element& element, const std::vector<double>& coefficients,
                                                                   const std::size_t moment_count, const std::size_t node_count, const double time,
                                                                   const qmomkit::PhaseSpaceVelocity& velocity,
                                                                   const qmomkit::QmomOptions& options = {})
{
    const std::vector<double> local_moments = interpolate_element_state(element, coefficients, moment_count);
    std::vector<double> local_rates(element.quadrature_count() * moment_count, 0.0);
    for (std::size_t q = 0; q < element.quadrature_count(); ++q)
    {
        const std::vector<double> moments(local_moments.begin() + static_cast<std::ptrdiff_t>(q * moment_count),
                                          local_moments.begin() + static_cast<std::ptrdiff_t>((q + 1) * moment_count));
        const qmomkit::QmomContinuousMomentRateResult result =
            qmomkit::evaluate_qmom_continuous_moment_rates(moments, node_count, time, velocity, options);
        if (!result.ok())
        {
            return continuous_failure(q, result.status, to_string(result.inversion_status), result.message);
        }
        for (std::size_t moment = 0; moment < moment_count; ++moment)
        {
            local_rates[q * moment_count + moment] = result.rates[moment];
        }
    }
    return success_result(project_rates_to_coefficients(element, local_rates, moment_count));
}

[[nodiscard]] inline ElementRhsResult evaluate_qmom_discontinuous_rhs(const Element& element, const std::vector<double>& coefficients,
                                                                      const std::size_t moment_count, const double time,
                                                                      const qmomkit::MomentSource& source)
{
    const std::vector<double> local_moments = interpolate_element_state(element, coefficients, moment_count);
    std::vector<double> local_rates(element.quadrature_count() * moment_count, 0.0);
    for (std::size_t q = 0; q < element.quadrature_count(); ++q)
    {
        const std::vector<double> moments(local_moments.begin() + static_cast<std::ptrdiff_t>(q * moment_count),
                                          local_moments.begin() + static_cast<std::ptrdiff_t>((q + 1) * moment_count));
        const qmomkit::QmomDiscontinuousMomentRateResult result = qmomkit::evaluate_qmom_discontinuous_moment_rates(time, moments, source);
        if (!result.ok())
        {
            return discontinuous_failure(q, result.status, result.message);
        }
        for (std::size_t moment = 0; moment < moment_count; ++moment)
        {
            local_rates[q * moment_count + moment] = result.rates[moment];
        }
    }
    return success_result(project_rates_to_coefficients(element, local_rates, moment_count));
}

[[nodiscard]] inline ElementRhsResult evaluate_cqmom_continuous_rhs(const Element& element, const std::vector<double>& coefficients,
                                                                    const qmomkit::CqmomLeanCase& template_case, const double time,
                                                                    const qmomkit::PhaseSpaceVelocity& velocity)
{
    const std::size_t moment_count = template_case.moments.values.size();
    const std::vector<double> local_moments = interpolate_element_state(element, coefficients, moment_count);
    std::vector<double> local_rates(element.quadrature_count() * moment_count, 0.0);
    for (std::size_t q = 0; q < element.quadrature_count(); ++q)
    {
        qmomkit::CqmomLeanCase local_case = template_case;
        local_case.moments.values.assign(local_moments.begin() + static_cast<std::ptrdiff_t>(q * moment_count),
                                         local_moments.begin() + static_cast<std::ptrdiff_t>((q + 1) * moment_count));
        const qmomkit::CqmomContinuousMomentRateResult result = qmomkit::evaluate_cqmom_continuous_moment_rates(local_case, time, velocity);
        if (!result.ok())
        {
            return continuous_failure(q, result.status, to_string(result.inversion_status), result.message);
        }
        for (std::size_t moment = 0; moment < moment_count; ++moment)
        {
            local_rates[q * moment_count + moment] = result.rates->values[moment];
        }
    }
    return success_result(project_rates_to_coefficients(element, local_rates, moment_count));
}

[[nodiscard]] inline ElementRhsResult evaluate_cqmom_discontinuous_rhs(const Element& element, const std::vector<double>& coefficients,
                                                                       const qmomkit::CqmomMomentData& template_moments, const double time,
                                                                       const qmomkit::CqmomMomentSource& source)
{
    const std::size_t moment_count = template_moments.values.size();
    const std::vector<double> local_moments = interpolate_element_state(element, coefficients, moment_count);
    std::vector<double> local_rates(element.quadrature_count() * moment_count, 0.0);
    for (std::size_t q = 0; q < element.quadrature_count(); ++q)
    {
        qmomkit::CqmomMomentData moments = template_moments;
        moments.values.assign(local_moments.begin() + static_cast<std::ptrdiff_t>(q * moment_count),
                              local_moments.begin() + static_cast<std::ptrdiff_t>((q + 1) * moment_count));
        const qmomkit::CqmomDiscontinuousMomentRateResult result = qmomkit::evaluate_cqmom_discontinuous_moment_rates(time, moments, source);
        if (!result.ok())
        {
            return discontinuous_failure(q, result.status, result.message);
        }
        for (std::size_t moment = 0; moment < moment_count; ++moment)
        {
            local_rates[q * moment_count + moment] = result.rates->values[moment];
        }
    }
    return success_result(project_rates_to_coefficients(element, local_rates, moment_count));
}

template <typename Rhs>
[[nodiscard]] inline ElementRhsResult ssp_rk3_coefficient_step(const std::vector<double>& coefficients, const double time, const double time_step,
                                                               Rhs rhs)
{
    ElementRhsResult k1 = rhs(time, coefficients);
    if (!k1.ok)
    {
        return k1;
    }

    std::vector<double> stage1(coefficients.size(), 0.0);
    for (std::size_t i = 0; i < coefficients.size(); ++i)
    {
        stage1[i] = coefficients[i] + time_step * k1.coefficient_rates[i];
    }

    ElementRhsResult k2 = rhs(time + time_step, stage1);
    if (!k2.ok)
    {
        return k2;
    }

    std::vector<double> stage2(coefficients.size(), 0.0);
    for (std::size_t i = 0; i < coefficients.size(); ++i)
    {
        stage2[i] = 0.75 * coefficients[i] + 0.25 * (stage1[i] + time_step * k2.coefficient_rates[i]);
    }

    ElementRhsResult k3 = rhs(time + 0.5 * time_step, stage2);
    if (!k3.ok)
    {
        return k3;
    }

    std::vector<double> updated(coefficients.size(), 0.0);
    for (std::size_t i = 0; i < coefficients.size(); ++i)
    {
        updated[i] = (coefficients[i] / 3.0) + (2.0 / 3.0) * (stage2[i] + time_step * k3.coefficient_rates[i]);
    }
    return success_result(std::move(updated));
}

[[nodiscard]] inline std::string describe_failure(const ElementRhsResult& result)
{
    std::ostringstream stream;
    stream << "quadrature_point=" << result.failed_quadrature_point << " status=" << result.status << " inversion_status=" << result.inversion_status
           << " message=" << result.message;
    return stream.str();
}

} // namespace moment_rate_element_operator
