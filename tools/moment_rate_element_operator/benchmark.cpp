#include "element_operator.hpp"

#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using namespace moment_rate_element_operator;

volatile double benchmark_sink = 0.0;

struct CaseMetadata
{
    std::string case_name;
    std::string category;
    std::size_t basis_size = 0;
    std::size_t quadrature_count = 0;
    std::size_t moment_count = 0;
    std::size_t phase_space_dimension = 0;
    std::string node_configuration;
    std::string diagnostic_level;
};

struct Timing
{
    double seconds = 0.0;
    std::size_t repetitions = 0;
    std::size_t local_calls_per_repetition = 0;
};

void consume(const std::vector<double>& values)
{
    double sum = 0.0;
    for (const double value : values)
    {
        sum += value;
    }
    benchmark_sink += sum;
}

template <typename Work> Timing measure(const std::size_t repetitions, const std::size_t local_calls_per_repetition, Work work)
{
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t rep = 0; rep < repetitions; ++rep)
    {
        consume(work());
    }
    const auto stop = std::chrono::steady_clock::now();
    return Timing{std::chrono::duration<double>(stop - start).count(), repetitions, local_calls_per_repetition};
}

void print_row(const CaseMetadata& metadata, const Timing& timing)
{
    const double average_ns = timing.seconds * 1.0e9 / static_cast<double>(timing.repetitions);
    const double average_local_call_ns = timing.local_calls_per_repetition == 0
                                             ? 0.0
                                             : timing.seconds * 1.0e9 / static_cast<double>(timing.repetitions * timing.local_calls_per_repetition);
    std::cout << metadata.case_name << ',' << metadata.category << ',' << metadata.basis_size << ',' << metadata.quadrature_count << ','
              << metadata.moment_count << ',' << metadata.phase_space_dimension << ',' << metadata.node_configuration << ','
              << metadata.diagnostic_level << ',' << timing.repetitions << ',' << timing.seconds << ',' << average_ns << ',' << average_local_call_ns
              << '\n';
}

std::size_t parse_repetitions(const int argc, char** argv)
{
    std::size_t repetitions = 1000;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--short")
        {
            repetitions = 25;
        }
        else if (arg == "--repetitions" && i + 1 < argc)
        {
            repetitions = static_cast<std::size_t>(std::stoull(argv[++i]));
        }
        else
        {
            throw std::invalid_argument("Usage: moment_rate_element_operator_benchmark [--short] [--repetitions N]");
        }
    }
    return repetitions;
}

void run_qmom_continuous_case(const Element& element, const std::size_t repetitions, const std::string& case_name, const std::vector<double>& moments,
                              const std::size_t node_count, const qmomkit::PhaseSpaceVelocity& velocity)
{
    const std::vector<double> coefficients = constant_coefficients(element, moments);
    const std::size_t moment_count = moments.size();
    const std::vector<double> local_moments = interpolate_element_state(element, coefficients, moment_count);
    const ElementRhsResult probe = evaluate_qmom_continuous_rhs(element, coefficients, moment_count, node_count, 0.0, velocity);
    if (!probe.ok)
    {
        throw std::runtime_error("benchmark setup failed: " + describe_failure(probe));
    }

    const CaseMetadata base{
        case_name,       "", element.basis_size(), element.quadrature_count(), moment_count, 1, std::to_string(node_count) + " QMOM nodes",
        "not_applicable"};
    print_row(CaseMetadata{base.case_name, "spatial_interpolation_only", base.basis_size, base.quadrature_count, base.moment_count,
                           base.phase_space_dimension, base.node_configuration, base.diagnostic_level},
              measure(repetitions, 0, [&] { return interpolate_element_state(element, coefficients, moment_count); }));
    print_row(CaseMetadata{base.case_name, "qmom_continuous_local_rate_calls_only", base.basis_size, base.quadrature_count, base.moment_count,
                           base.phase_space_dimension, base.node_configuration, base.diagnostic_level},
              measure(repetitions, element.quadrature_count(),
                      [&]
                      {
                          std::vector<double> rates(element.quadrature_count() * moment_count, 0.0);
                          for (std::size_t q = 0; q < element.quadrature_count(); ++q)
                          {
                              const std::vector<double> local(local_moments.begin() + static_cast<std::ptrdiff_t>(q * moment_count),
                                                              local_moments.begin() + static_cast<std::ptrdiff_t>((q + 1) * moment_count));
                              const qmomkit::QmomContinuousMomentRateResult result =
                                  qmomkit::evaluate_qmom_continuous_moment_rates(local, node_count, 0.0, velocity);
                              for (std::size_t moment = 0; moment < moment_count; ++moment)
                              {
                                  rates[q * moment_count + moment] = result.rates[moment];
                              }
                          }
                          return rates;
                      }));
    const std::vector<double> local_rates = interpolate_element_state(element, probe.coefficient_rates, moment_count);
    print_row(CaseMetadata{base.case_name, "projection_and_inverse_mass_only", base.basis_size, base.quadrature_count, base.moment_count,
                           base.phase_space_dimension, base.node_configuration, base.diagnostic_level},
              measure(repetitions, 0, [&] { return project_rates_to_coefficients(element, local_rates, moment_count); }));
    print_row(CaseMetadata{base.case_name, "complete_element_rhs", base.basis_size, base.quadrature_count, base.moment_count,
                           base.phase_space_dimension, base.node_configuration, base.diagnostic_level},
              measure(repetitions, element.quadrature_count(),
                      [&]
                      {
                          ElementRhsResult result = evaluate_qmom_continuous_rhs(element, coefficients, moment_count, node_count, 0.0, velocity);
                          return result.coefficient_rates;
                      }));
    print_row(CaseMetadata{base.case_name, "complete_external_ssp_rk3_step", base.basis_size, base.quadrature_count, base.moment_count,
                           base.phase_space_dimension, base.node_configuration, base.diagnostic_level},
              measure(repetitions, 3 * element.quadrature_count(),
                      [&]
                      {
                          ElementRhsResult result = ssp_rk3_coefficient_step(
                              coefficients, 0.0, 0.01, [&](const double time, const std::vector<double>& state)
                              { return evaluate_qmom_continuous_rhs(element, state, moment_count, node_count, time, velocity); });
                          return result.coefficient_rates;
                      }));
}

void run_qmom_discontinuous_case(const Element& element, const std::size_t repetitions)
{
    const std::size_t moment_count = 4;
    const std::vector<double> coefficients = constant_coefficients(element, qmom_moments_from_terminals({0.4, 0.6}, {-0.5, 1.2}, moment_count));
    const qmomkit::MomentSource source = [](const double time, const std::vector<double>& state)
    {
        std::vector<double> rates(state.size(), 0.0);
        for (std::size_t moment = 0; moment < state.size(); ++moment)
        {
            rates[moment] = 0.1 + 0.01 * time + 0.02 * static_cast<double>(moment) + 0.005 * state[moment];
        }
        return rates;
    };
    const std::vector<double> local_moments = interpolate_element_state(element, coefficients, moment_count);
    const ElementRhsResult probe = evaluate_qmom_discontinuous_rhs(element, coefficients, moment_count, 0.0, source);
    if (!probe.ok)
    {
        throw std::runtime_error("benchmark setup failed: " + describe_failure(probe));
    }
    const CaseMetadata base{
        "qmom_discontinuous_affine_source", "", element.basis_size(), element.quadrature_count(), moment_count, 1, "2 QMOM nodes", "not_applicable"};
    print_row(CaseMetadata{base.case_name, "qmom_discontinuous_local_rate_calls_only", base.basis_size, base.quadrature_count, base.moment_count,
                           base.phase_space_dimension, base.node_configuration, base.diagnostic_level},
              measure(repetitions, element.quadrature_count(),
                      [&]
                      {
                          std::vector<double> rates(element.quadrature_count() * moment_count, 0.0);
                          for (std::size_t q = 0; q < element.quadrature_count(); ++q)
                          {
                              const std::vector<double> local(local_moments.begin() + static_cast<std::ptrdiff_t>(q * moment_count),
                                                              local_moments.begin() + static_cast<std::ptrdiff_t>((q + 1) * moment_count));
                              const qmomkit::QmomDiscontinuousMomentRateResult result =
                                  qmomkit::evaluate_qmom_discontinuous_moment_rates(0.0, local, source);
                              for (std::size_t moment = 0; moment < moment_count; ++moment)
                              {
                                  rates[q * moment_count + moment] = result.rates[moment];
                              }
                          }
                          return rates;
                      }));
    print_row(CaseMetadata{base.case_name, "complete_element_rhs", base.basis_size, base.quadrature_count, base.moment_count,
                           base.phase_space_dimension, base.node_configuration, base.diagnostic_level},
              measure(repetitions, element.quadrature_count(),
                      [&]
                      {
                          ElementRhsResult result = evaluate_qmom_discontinuous_rhs(element, coefficients, moment_count, 0.0, source);
                          return result.coefficient_rates;
                      }));
    print_row(CaseMetadata{base.case_name, "complete_external_ssp_rk3_step", base.basis_size, base.quadrature_count, base.moment_count,
                           base.phase_space_dimension, base.node_configuration, base.diagnostic_level},
              measure(repetitions, 3 * element.quadrature_count(),
                      [&]
                      {
                          ElementRhsResult result =
                              ssp_rk3_coefficient_step(coefficients, 0.0, 0.01, [&](const double time, const std::vector<double>& state)
                                                       { return evaluate_qmom_discontinuous_rhs(element, state, moment_count, time, source); });
                          return result.coefficient_rates;
                      }));
}

void run_cqmom_continuous_case(const Element& element, const std::size_t repetitions, const std::string& case_name,
                               qmomkit::CqmomDiagnosticLevel diagnostic_level, const std::string& order_name)
{
    const qmomkit::PhaseSpaceHierarchy hierarchy =
        order_name == "identity" ? qmomkit::PhaseSpaceHierarchy({2, 2}) : qmomkit::PhaseSpaceHierarchy({2, 2}, {1, 0});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::CqmomLeanCase lean_case{cqmom_moment_data_from_terminals(layout, {0.2, 0.3, 0.1, 0.4}, {{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}}),
                                     {},
                                     {},
                                     false,
                                     diagnostic_level,
                                     qmomkit::VandermondeSolver::automatic};
    const std::vector<double> coefficients = constant_coefficients(element, lean_case.moments.values);
    const qmomkit::PhaseSpaceVelocity velocity = [](double, const std::vector<double>& coordinate)
    { return std::vector<double>{0.05 + 0.02 * coordinate[0] * coordinate[0], -0.03 * coordinate[1]}; };
    const std::vector<double> local_moments = interpolate_element_state(element, coefficients, lean_case.moments.values.size());
    const ElementRhsResult probe = evaluate_cqmom_continuous_rhs(element, coefficients, lean_case, 0.0, velocity);
    if (!probe.ok)
    {
        throw std::runtime_error("benchmark setup failed: " + describe_failure(probe));
    }

    const std::string diagnostic_name = diagnostic_level == qmomkit::CqmomDiagnosticLevel::none ? "none" : "detailed";
    const CaseMetadata base{
        case_name, "", element.basis_size(), element.quadrature_count(), lean_case.moments.values.size(), 2, order_name + " 2x2", diagnostic_name};
    print_row(CaseMetadata{base.case_name, "cqmom_continuous_local_rate_calls_only", base.basis_size, base.quadrature_count, base.moment_count,
                           base.phase_space_dimension, base.node_configuration, base.diagnostic_level},
              measure(repetitions, element.quadrature_count(),
                      [&]
                      {
                          std::vector<double> rates(element.quadrature_count() * lean_case.moments.values.size(), 0.0);
                          for (std::size_t q = 0; q < element.quadrature_count(); ++q)
                          {
                              qmomkit::CqmomLeanCase local_case = lean_case;
                              local_case.moments.values.assign(
                                  local_moments.begin() + static_cast<std::ptrdiff_t>(q * lean_case.moments.values.size()),
                                  local_moments.begin() + static_cast<std::ptrdiff_t>((q + 1) * lean_case.moments.values.size()));
                              const qmomkit::CqmomContinuousMomentRateResult result =
                                  qmomkit::evaluate_cqmom_continuous_moment_rates(local_case, 0.0, velocity);
                              for (std::size_t moment = 0; moment < lean_case.moments.values.size(); ++moment)
                              {
                                  rates[q * lean_case.moments.values.size() + moment] = result.rates->values[moment];
                              }
                          }
                          return rates;
                      }));
    print_row(CaseMetadata{base.case_name, "complete_element_rhs", base.basis_size, base.quadrature_count, base.moment_count,
                           base.phase_space_dimension, base.node_configuration, base.diagnostic_level},
              measure(repetitions, element.quadrature_count(),
                      [&]
                      {
                          ElementRhsResult result = evaluate_cqmom_continuous_rhs(element, coefficients, lean_case, 0.0, velocity);
                          return result.coefficient_rates;
                      }));
}

void run_cqmom_discontinuous_case(const Element& element, const std::size_t repetitions)
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData moments =
        cqmom_moment_data_from_terminals(layout, {0.2, 0.3, 0.1, 0.4}, {{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}});
    const std::vector<double> coefficients = constant_coefficients(element, moments.values);
    const qmomkit::CqmomMomentSource source = [](const double time, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rates(context.size(), 0.0);
        for (std::size_t moment = 0; moment < context.size(); ++moment)
        {
            const qmomkit::MultiIndex& index = context.index(moment);
            rates[moment] =
                0.02 + 0.01 * time + 0.03 * static_cast<double>(index[0]) - 0.01 * static_cast<double>(index[1]) + 0.002 * context.value(moment);
        }
        return qmomkit::CqmomMomentRate{rates};
    };
    const std::vector<double> local_moments = interpolate_element_state(element, coefficients, moments.values.size());
    const ElementRhsResult probe = evaluate_cqmom_discontinuous_rhs(element, coefficients, moments, 0.0, source);
    if (!probe.ok)
    {
        throw std::runtime_error("benchmark setup failed: " + describe_failure(probe));
    }
    const CaseMetadata base{"cqmom_discontinuous_affine_source",
                            "",
                            element.basis_size(),
                            element.quadrature_count(),
                            moments.values.size(),
                            2,
                            "nonidentity 2x2",
                            "not_applicable"};
    print_row(CaseMetadata{base.case_name, "cqmom_discontinuous_local_rate_calls_only", base.basis_size, base.quadrature_count, base.moment_count,
                           base.phase_space_dimension, base.node_configuration, base.diagnostic_level},
              measure(repetitions, element.quadrature_count(),
                      [&]
                      {
                          std::vector<double> rates(element.quadrature_count() * moments.values.size(), 0.0);
                          for (std::size_t q = 0; q < element.quadrature_count(); ++q)
                          {
                              qmomkit::CqmomMomentData local_moments_data = moments;
                              local_moments_data.values.assign(local_moments.begin() + static_cast<std::ptrdiff_t>(q * moments.values.size()),
                                                               local_moments.begin() + static_cast<std::ptrdiff_t>((q + 1) * moments.values.size()));
                              const qmomkit::CqmomDiscontinuousMomentRateResult result =
                                  qmomkit::evaluate_cqmom_discontinuous_moment_rates(0.0, local_moments_data, source);
                              for (std::size_t moment = 0; moment < moments.values.size(); ++moment)
                              {
                                  rates[q * moments.values.size() + moment] = result.rates->values[moment];
                              }
                          }
                          return rates;
                      }));
    print_row(CaseMetadata{base.case_name, "complete_element_rhs", base.basis_size, base.quadrature_count, base.moment_count,
                           base.phase_space_dimension, base.node_configuration, base.diagnostic_level},
              measure(repetitions, element.quadrature_count(),
                      [&]
                      {
                          ElementRhsResult result = evaluate_cqmom_discontinuous_rhs(element, coefficients, moments, 0.0, source);
                          return result.coefficient_rates;
                      }));
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const std::size_t repetitions = parse_repetitions(argc, argv);
        const Element quadratic = make_lagrange_element(2);
        const Element cubic = make_lagrange_element(3);
        std::cout << "case_name,category,basis_size,spatial_quadrature_count,transported_moment_count,phase_space_dimension,node_configuration,"
                     "diagnostic_level,repetitions,total_elapsed_seconds,average_time_per_element_rhs_ns,average_time_per_local_rate_call_ns\n";
        run_qmom_continuous_case(quadratic, repetitions, "qmom_continuous_one_node_affine_callback", qmom_moments_from_terminals({2.0}, {3.0}, 4), 1,
                                 [](double, const std::vector<double>& coordinate) { return std::vector<double>{0.05 + 0.01 * coordinate[0]}; });
        run_qmom_continuous_case(cubic, repetitions, "qmom_continuous_two_node_nonlinear_callback",
                                 qmom_moments_from_terminals({0.4, 0.6}, {-0.5, 1.2}, 4), 2, [](double, const std::vector<double>& coordinate)
                                 { return std::vector<double>{0.05 + 0.02 * coordinate[0] * coordinate[0]}; });
        run_qmom_discontinuous_case(quadratic, repetitions);
        run_cqmom_continuous_case(quadratic, repetitions, "cqmom_continuous_identity_none", qmomkit::CqmomDiagnosticLevel::none, "identity");
        run_cqmom_continuous_case(quadratic, repetitions, "cqmom_continuous_nonidentity_detailed", qmomkit::CqmomDiagnosticLevel::detailed,
                                  "nonidentity");
        run_cqmom_discontinuous_case(quadratic, repetitions);
        if (benchmark_sink == 123456789.0)
        {
            std::cerr << "unreachable sink value: " << benchmark_sink << '\n';
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << "moment-rate element-operator benchmark failed: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
