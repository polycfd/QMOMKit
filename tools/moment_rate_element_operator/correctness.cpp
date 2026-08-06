#include "element_operator.hpp"

#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using namespace moment_rate_element_operator;

void require(const bool condition, const std::string& message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

void require_near(const double actual, const double expected, const double tolerance, const std::string& message)
{
    if (std::abs(actual - expected) > tolerance)
    {
        throw std::runtime_error(message + ": actual=" + std::to_string(actual) + " expected=" + std::to_string(expected));
    }
}

void require_vector_near(const std::vector<double>& actual, const std::vector<double>& expected, const double tolerance, const std::string& message)
{
    require(actual.size() == expected.size(), message + ": size mismatch");
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        require_near(actual[i], expected[i], tolerance, message + " at index " + std::to_string(i));
    }
}

std::vector<double> qmom_single_node_advection_rates(const std::vector<double>& moments, const double velocity)
{
    const double weight = moments[0];
    const double abscissa = moments[1] / moments[0];
    std::vector<double> rates(moments.size(), 0.0);
    for (std::size_t order = 1; order < moments.size(); ++order)
    {
        rates[order] = static_cast<double>(order) * weight * velocity * integer_power(abscissa, order - 1);
    }
    return rates;
}

std::vector<double> manual_qmom_moment_source_reference(const Element& element, const std::vector<double>& coefficients,
                                                        const std::size_t moment_count, const double time)
{
    const std::vector<double> local_moments = interpolate_element_state(element, coefficients, moment_count);
    std::vector<double> local_rates(element.quadrature_count() * moment_count, 0.0);
    for (std::size_t q = 0; q < element.quadrature_count(); ++q)
    {
        for (std::size_t moment = 0; moment < moment_count; ++moment)
        {
            const double value = local_moments[q * moment_count + moment];
            local_rates[q * moment_count + moment] = 0.25 + 0.1 * static_cast<double>(moment + 1) + 0.05 * time + 0.02 * value;
        }
    }
    return project_rates_to_coefficients(element, local_rates, moment_count);
}

std::vector<double> manual_cqmom_moment_source_reference(const Element& element, const std::vector<double>& coefficients,
                                                         const qmomkit::CqmomMomentLayout& layout, const double time)
{
    const std::size_t moment_count = layout.required_indices.size();
    const std::vector<double> local_moments = interpolate_element_state(element, coefficients, moment_count);
    std::vector<double> local_rates(element.quadrature_count() * moment_count, 0.0);
    for (std::size_t q = 0; q < element.quadrature_count(); ++q)
    {
        for (std::size_t moment = 0; moment < moment_count; ++moment)
        {
            const qmomkit::MultiIndex& index = layout.required_indices[moment];
            const double value = local_moments[q * moment_count + moment];
            local_rates[q * moment_count + moment] =
                0.1 + time + 0.25 * static_cast<double>(index[0]) - 0.15 * static_cast<double>(index[1]) + 0.01 * value;
        }
    }
    return project_rates_to_coefficients(element, local_rates, moment_count);
}

void check_qmom_advection()
{
    const Element element = make_lagrange_element(2);
    const std::vector<double> moments = qmom_moments_from_terminals({2.0}, {3.0}, 4);
    const std::vector<double> coefficients = constant_coefficients(element, moments);
    const double velocity_value = 0.5;
    const qmomkit::PhaseSpaceVelocity velocity = [velocity_value](double, const std::vector<double>&) { return std::vector<double>{velocity_value}; };

    const ElementRhsResult result = evaluate_qmom_advection_rhs(element, coefficients, moments.size(), 1, 0.0, velocity);
    require(result.ok, "QMOM advection RHS failed: " + describe_failure(result));

    const std::vector<double> local_expected = qmom_single_node_advection_rates(moments, velocity_value);
    require_vector_near(result.coefficient_rates, constant_coefficients(element, local_expected), 1.0e-12, "QMOM advection coefficient rate");
}

void check_qmom_moment_source()
{
    const Element element = make_lagrange_element(2);
    std::vector<double> coefficients = constant_coefficients(element, qmom_moments_from_terminals({0.4, 0.6}, {-0.5, 1.2}, 4));
    for (std::size_t basis = 0; basis < element.basis_size(); ++basis)
    {
        coefficients[coefficient_index(basis, 1, 4)] += 0.03 * static_cast<double>(basis);
    }
    const qmomkit::MomentSource source = [](const double time, const std::vector<double>& state)
    {
        std::vector<double> rates(state.size(), 0.0);
        for (std::size_t moment = 0; moment < state.size(); ++moment)
        {
            rates[moment] = 0.25 + 0.1 * static_cast<double>(moment + 1) + 0.05 * time + 0.02 * state[moment];
        }
        return rates;
    };

    const ElementRhsResult result = evaluate_qmom_moment_source_rhs(element, coefficients, 4, 0.25, source);
    require(result.ok, "QMOM moment-source RHS failed: " + describe_failure(result));
    require_vector_near(result.coefficient_rates, manual_qmom_moment_source_reference(element, coefficients, 4, 0.25), 1.0e-12,
                        "QMOM moment-source coefficient rate");
}

void check_cqmom_advection()
{
    const Element element = make_lagrange_element(2);
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::CqmomLeanCase lean_case{cqmom_moment_data_from_terminals(layout, {0.2, 0.3, 0.1, 0.4}, {{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}}),
                                     {},
                                     {},
                                     false,
                                     qmomkit::CqmomDiagnosticLevel::none,
                                     qmomkit::VandermondeSolver::automatic};
    const std::vector<double> coefficients = constant_coefficients(element, lean_case.moments.values);
    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>& coordinate)
    { return std::vector<double>{0.1 + 0.05 * coordinate[0], -0.2 * coordinate[1]}; };

    const ElementRhsResult result = evaluate_cqmom_advection_rhs(element, coefficients, lean_case, 0.0, velocity);
    require(result.ok, "CQMOM advection RHS failed: " + describe_failure(result));

    const std::vector<double> local_moments = interpolate_element_state(element, coefficients, lean_case.moments.values.size());
    std::vector<double> local_rates(element.quadrature_count() * lean_case.moments.values.size(), 0.0);
    for (std::size_t q = 0; q < element.quadrature_count(); ++q)
    {
        qmomkit::CqmomLeanCase local_case = lean_case;
        local_case.moments.values.assign(local_moments.begin() + static_cast<std::ptrdiff_t>(q * lean_case.moments.values.size()),
                                         local_moments.begin() + static_cast<std::ptrdiff_t>((q + 1) * lean_case.moments.values.size()));
        const qmomkit::CqmomAdvectionMomentRateResult local = qmomkit::evaluate_cqmom_advection_moment_rates(local_case, 0.0, velocity);
        require(local.ok(), "manual CQMOM advection local call failed");
        for (std::size_t moment = 0; moment < lean_case.moments.values.size(); ++moment)
        {
            local_rates[q * lean_case.moments.values.size() + moment] = local.rates->values[moment];
        }
    }
    require_vector_near(result.coefficient_rates, project_rates_to_coefficients(element, local_rates, lean_case.moments.values.size()), 1.0e-12,
                        "CQMOM advection coefficient rate");
}

void check_cqmom_moment_source()
{
    const Element element = make_lagrange_element(2);
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData moments =
        cqmom_moment_data_from_terminals(layout, {0.2, 0.3, 0.1, 0.4}, {{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}});
    const std::vector<double> coefficients = constant_coefficients(element, moments.values);
    bool inspected_context = false;
    const qmomkit::CqmomMomentSource source = [&layout, &inspected_context](const double time, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rates(context.size(), 0.0);
        for (std::size_t moment = 0; moment < context.size(); ++moment)
        {
            require(context.index(moment) == layout.required_indices[moment], "CQMOM source context ordering changed");
            const qmomkit::MultiIndex& index = context.index(moment);
            rates[moment] = 0.1 + time + 0.25 * static_cast<double>(index[0]) - 0.15 * static_cast<double>(index[1]) + 0.01 * context.value(moment);
        }
        inspected_context = true;
        return qmomkit::CqmomMomentRate{rates};
    };

    const ElementRhsResult result = evaluate_cqmom_moment_source_rhs(element, coefficients, moments, 0.25, source);
    require(result.ok, "CQMOM moment-source RHS failed: " + describe_failure(result));
    require(inspected_context, "CQMOM moment source was not invoked");
    require_vector_near(result.coefficient_rates, manual_cqmom_moment_source_reference(element, coefficients, layout, 0.25), 1.0e-12,
                        "CQMOM moment-source coefficient rate");
}

void check_external_ssp_rk3()
{
    const Element element = make_lagrange_element(2);
    const std::vector<double> qmom_moments = qmom_moments_from_terminals({2.0}, {3.0}, 4);
    const std::vector<double> qmom_coefficients = constant_coefficients(element, qmom_moments);
    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>&) { return std::vector<double>{0.02}; };
    const ElementRhsResult advection_step =
        ssp_rk3_coefficient_step(qmom_coefficients, 0.0, 0.05, [&element, &velocity](const double time, const std::vector<double>& coefficients)
                                 { return evaluate_qmom_advection_rhs(element, coefficients, 4, 1, time, velocity); });
    require(advection_step.ok, "external advection SSP-RK3 step failed: " + describe_failure(advection_step));

    const qmomkit::MomentSource source = [](double, const std::vector<double>& moments) { return std::vector<double>(moments.size(), 0.01); };
    const ElementRhsResult moment_source_step =
        ssp_rk3_coefficient_step(qmom_coefficients, 0.0, 0.05, [&element, &source](const double time, const std::vector<double>& coefficients)
                                 { return evaluate_qmom_moment_source_rhs(element, coefficients, 4, time, source); });
    require(moment_source_step.ok, "external moment-source SSP-RK3 stage failed: " + describe_failure(moment_source_step));
}

void check_failure_diagnostics()
{
    const Element element = make_lagrange_element(2);
    const std::vector<double> nonrealizable{1.0, 0.0, -1.0, 0.0};
    const std::vector<double> coefficients = constant_coefficients(element, nonrealizable);
    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>&) { return std::vector<double>{0.0}; };
    const ElementRhsResult result = evaluate_qmom_advection_rhs(element, coefficients, nonrealizable.size(), 2, 0.0, velocity);
    require(!result.ok, "expected failed local advection rate call");
    require(result.failed_quadrature_point == 0, "failed quadrature-point index was not reported");
    require(result.status == "inversion_failed", "public rate status was not reported");
    require(result.inversion_status == "not_realizable", "underlying inversion status was not reported");
}

} // namespace

int main()
{
    try
    {
        check_qmom_advection();
        check_qmom_moment_source();
        check_cqmom_advection();
        check_cqmom_moment_source();
        check_external_ssp_rk3();
        check_failure_diagnostics();
    }
    catch (const std::exception& error)
    {
        std::cerr << "moment-rate element-operator check failed: " << error.what() << '\n';
        return 1;
    }

    std::cout << "moment-rate element-operator checks passed\n";
    return 0;
}
