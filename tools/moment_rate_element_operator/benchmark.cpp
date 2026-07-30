#include "element_operator.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef MOMENT_RATE_ELEMENT_OPERATOR_COMPILER_ID
#define MOMENT_RATE_ELEMENT_OPERATOR_COMPILER_ID "unknown"
#endif

#ifndef MOMENT_RATE_ELEMENT_OPERATOR_COMPILER_VERSION
#define MOMENT_RATE_ELEMENT_OPERATOR_COMPILER_VERSION "unknown"
#endif

#ifndef MOMENT_RATE_ELEMENT_OPERATOR_BUILD_TYPE
#define MOMENT_RATE_ELEMENT_OPERATOR_BUILD_TYPE "unknown"
#endif

#ifndef MOMENT_RATE_ELEMENT_OPERATOR_SYSTEM_NAME
#define MOMENT_RATE_ELEMENT_OPERATOR_SYSTEM_NAME "unknown"
#endif

namespace
{

using namespace moment_rate_element_operator;

volatile double benchmark_sink = 0.0;

struct BenchmarkOptions
{
    bool csv = false;
    bool short_run = false;
    std::size_t samples = 7;
    std::size_t warmup_iterations = 4;
    std::size_t min_iterations = 1;
    std::size_t max_iterations = 1U << 20U;
    double target_sample_seconds = 0.020;
    std::vector<std::size_t> batch_sizes{1, 100, 1000};
    std::string filter;
};

struct CaseMetadata
{
    std::string case_name;
    std::string category;
    std::size_t basis_size = 0;
    std::size_t quadrature_count = 0;
    std::size_t moment_count = 0;
    std::size_t phase_space_dimension = 0;
    std::string node_configuration;
    std::string conditioning_order;
    std::string diagnostic_level;
    std::string callback_label;
    std::size_t element_batch_size = 1;
    std::size_t local_calls_per_element = 0;
    std::size_t local_calls_per_ssp_rk3 = 0;
};

struct BenchmarkCase
{
    CaseMetadata metadata;
    std::function<void()> validate;
    std::function<std::vector<double>()> work;
};

struct TimingResult
{
    std::size_t sample_count = 0;
    std::size_t iterations_per_sample = 0;
    double median_sample_seconds = 0.0;
    double min_sample_seconds = 0.0;
    double max_sample_seconds = 0.0;
    double checksum = 0.0;
};

[[nodiscard]] std::string assertion_label()
{
#ifdef NDEBUG
    return "disabled";
#else
    return "enabled";
#endif
}

[[nodiscard]] double checksum_values(const std::vector<double>& values)
{
    double checksum = 0.0;
    double scale = 1.0;
    for (const double value : values)
    {
        checksum += scale * value;
        scale += 0.001;
    }
    return checksum;
}

void require_finite(const std::vector<double>& values, const std::string& context)
{
    if (values.empty())
    {
        throw std::runtime_error(context + " produced an empty result.");
    }
    for (const double value : values)
    {
        if (!std::isfinite(value))
        {
            throw std::runtime_error(context + " produced a non-finite value.");
        }
    }
}

[[nodiscard]] std::vector<double> local_moment_slice(const std::vector<double>& local_moments, const std::size_t q, const std::size_t moment_count)
{
    return std::vector<double>(local_moments.begin() + static_cast<std::ptrdiff_t>(q * moment_count),
                               local_moments.begin() + static_cast<std::ptrdiff_t>((q + 1) * moment_count));
}

[[nodiscard]] std::vector<double> repeat_batch(const std::size_t batch_size, const std::function<std::vector<double>()>& single_element_work)
{
    std::vector<double> accumulated;
    for (std::size_t element = 0; element < batch_size; ++element)
    {
        const std::vector<double> values = single_element_work();
        if (accumulated.empty())
        {
            accumulated.assign(values.size(), 0.0);
        }
        if (values.size() != accumulated.size())
        {
            throw std::runtime_error("Benchmark work returned inconsistent result sizes.");
        }
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            accumulated[i] += values[i];
        }
    }
    return accumulated;
}

template <typename Work> [[nodiscard]] std::pair<double, double> time_iterations(const std::size_t iterations, Work work)
{
    double checksum = 0.0;
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t iteration = 0; iteration < iterations; ++iteration)
    {
        checksum += checksum_values(work());
    }
    const auto stop = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(stop - start).count();
    if (!std::isfinite(checksum))
    {
        throw std::runtime_error("Benchmark checksum is non-finite.");
    }
    return {seconds, checksum};
}

[[nodiscard]] std::size_t calibrated_iterations(const BenchmarkOptions& options, const std::function<std::vector<double>()>& work)
{
    std::size_t iterations = std::max<std::size_t>(1, options.min_iterations);
    while (true)
    {
        const auto measured = time_iterations(iterations, work);
        benchmark_sink += measured.second;
        if (measured.first >= options.target_sample_seconds || iterations >= options.max_iterations)
        {
            return iterations;
        }
        const double safe_elapsed = std::max(measured.first, 1.0e-9);
        const double raw_factor = options.target_sample_seconds / safe_elapsed;
        const double factor = std::max(2.0, std::min(20.0, raw_factor * 1.20));
        const std::size_t next = static_cast<std::size_t>(std::ceil(static_cast<double>(iterations) * factor));
        iterations = std::min(options.max_iterations, std::max(iterations + 1, next));
    }
}

[[nodiscard]] TimingResult measure_case(const BenchmarkOptions& options, const BenchmarkCase& benchmark_case)
{
    if (benchmark_case.validate)
    {
        benchmark_case.validate();
    }
    require_finite(benchmark_case.work(), benchmark_case.metadata.case_name);

    for (std::size_t i = 0; i < options.warmup_iterations; ++i)
    {
        benchmark_sink += checksum_values(benchmark_case.work());
    }

    const std::size_t iterations = calibrated_iterations(options, benchmark_case.work);
    std::vector<double> samples;
    samples.reserve(options.samples);
    double checksum = 0.0;
    for (std::size_t sample = 0; sample < options.samples; ++sample)
    {
        const auto measured = time_iterations(iterations, benchmark_case.work);
        samples.push_back(measured.first);
        checksum += measured.second;
    }
    benchmark_sink += checksum;

    std::vector<double> sorted = samples;
    std::sort(sorted.begin(), sorted.end());
    return TimingResult{samples.size(),
                        iterations,
                        sorted[sorted.size() / 2],
                        *std::min_element(samples.begin(), samples.end()),
                        *std::max_element(samples.begin(), samples.end()),
                        checksum};
}

[[nodiscard]] double ns(const double seconds)
{
    return seconds * 1.0e9;
}

[[nodiscard]] double per_element_seconds(const TimingResult& timing, const CaseMetadata& metadata)
{
    return timing.median_sample_seconds / static_cast<double>(timing.iterations_per_sample * metadata.element_batch_size);
}

[[nodiscard]] double per_local_call_seconds(const TimingResult& timing, const CaseMetadata& metadata)
{
    const std::size_t local_calls_per_iteration = metadata.local_calls_per_element * metadata.element_batch_size;
    if (local_calls_per_iteration == 0)
    {
        return 0.0;
    }
    return timing.median_sample_seconds / static_cast<double>(timing.iterations_per_sample * local_calls_per_iteration);
}

void print_csv_header()
{
    std::cout << "compiler_id,compiler_version,build_type,assertions,os,target_sample_ms,case_name,operation_category,basis_size,"
                 "spatial_quadrature_count,transported_moment_count,phase_space_dimension,"
                 "node_configuration,conditioning_order,diagnostic_level,callback_label,element_batch_size,local_calls_per_element_rhs,"
                 "local_calls_per_ssp_rk3,sample_count,iterations_per_sample,median_total_time_ns,median_time_per_element_ns,"
                 "median_time_per_local_call_ns,min_sample_time_ns,max_sample_time_ns,checksum\n";
}

void print_csv_row(const BenchmarkOptions& options, const CaseMetadata& metadata, const TimingResult& timing)
{
    std::cout << MOMENT_RATE_ELEMENT_OPERATOR_COMPILER_ID << ',' << MOMENT_RATE_ELEMENT_OPERATOR_COMPILER_VERSION << ','
              << MOMENT_RATE_ELEMENT_OPERATOR_BUILD_TYPE << ',' << assertion_label() << ',' << MOMENT_RATE_ELEMENT_OPERATOR_SYSTEM_NAME << ','
              << std::fixed << std::setprecision(3) << options.target_sample_seconds * 1.0e3 << ',' << metadata.case_name << ',' << metadata.category
              << ',' << metadata.basis_size << ',' << metadata.quadrature_count << ',' << metadata.moment_count << ','
              << metadata.phase_space_dimension << ',' << metadata.node_configuration << ',' << metadata.conditioning_order << ','
              << metadata.diagnostic_level << ',' << metadata.callback_label << ',' << metadata.element_batch_size << ','
              << metadata.local_calls_per_element << ',' << metadata.local_calls_per_ssp_rk3 << ',' << timing.sample_count << ','
              << timing.iterations_per_sample << ',' << std::setprecision(1) << ns(timing.median_sample_seconds) << ','
              << ns(per_element_seconds(timing, metadata)) << ',' << ns(per_local_call_seconds(timing, metadata)) << ','
              << ns(timing.min_sample_seconds) << ',' << ns(timing.max_sample_seconds) << ',' << std::setprecision(12) << timing.checksum << '\n';
}

void print_human_metadata(const BenchmarkOptions& options)
{
    std::cout << "Moment-rate element-operator benchmark\n";
    std::cout << "Results are serial, machine-specific, and build-specific.\n";
    std::cout << "compiler: " << MOMENT_RATE_ELEMENT_OPERATOR_COMPILER_ID << ' ' << MOMENT_RATE_ELEMENT_OPERATOR_COMPILER_VERSION << '\n';
    std::cout << "build_type: " << MOMENT_RATE_ELEMENT_OPERATOR_BUILD_TYPE << '\n';
    std::cout << "assertions: " << assertion_label() << '\n';
    std::cout << "os: " << MOMENT_RATE_ELEMENT_OPERATOR_SYSTEM_NAME << '\n';
    std::cout << "samples: " << options.samples << '\n';
    std::cout << "target_sample_ms: " << options.target_sample_seconds * 1.0e3 << '\n';
    std::cout << "warmup_iterations: " << options.warmup_iterations << '\n';
    std::cout << '\n';
    std::cout << "case | operation | batch | iterations/sample | median_batch_ns | median_element_ns | median_local_call_ns | "
                 "min_batch_ns | max_batch_ns\n";
}

void print_human_row(const CaseMetadata& metadata, const TimingResult& timing)
{
    const double total_ns = ns(timing.median_sample_seconds / static_cast<double>(timing.iterations_per_sample));
    const double element_ns = ns(per_element_seconds(timing, metadata));
    const double call_ns = ns(per_local_call_seconds(timing, metadata));
    const double min_ns = ns(timing.min_sample_seconds / static_cast<double>(timing.iterations_per_sample));
    const double max_ns = ns(timing.max_sample_seconds / static_cast<double>(timing.iterations_per_sample));
    std::cout << metadata.case_name << " | " << metadata.category << " | " << metadata.element_batch_size << " | " << timing.iterations_per_sample
              << " | " << std::fixed << std::setprecision(1) << total_ns << " | " << element_ns << " | " << call_ns << " | " << min_ns << " | "
              << max_ns << '\n';
}

[[nodiscard]] qmomkit::PhaseSpaceVelocity zero_velocity()
{
    return [](double, const std::vector<double>& coordinate) { return std::vector<double>(coordinate.size(), 0.0); };
}

[[nodiscard]] qmomkit::PhaseSpaceVelocity affine_velocity()
{
    return [](const double time, const std::vector<double>& coordinate)
    {
        std::vector<double> velocity(coordinate.size(), 0.0);
        for (std::size_t dimension = 0; dimension < coordinate.size(); ++dimension)
        {
            velocity[dimension] = 0.03 + 0.01 * time + 0.02 * static_cast<double>(dimension + 1) * coordinate[dimension];
        }
        return velocity;
    };
}

[[nodiscard]] qmomkit::PhaseSpaceVelocity nonlinear_velocity()
{
    return [](const double time, const std::vector<double>& coordinate)
    {
        std::vector<double> velocity(coordinate.size(), 0.0);
        for (std::size_t dimension = 0; dimension < coordinate.size(); ++dimension)
        {
            const double x = coordinate[dimension];
            velocity[dimension] = 0.04 + 0.01 * time + 0.015 * x * x + 0.002 * x * static_cast<double>(dimension + 1);
        }
        return velocity;
    };
}

[[nodiscard]] qmomkit::MomentSource free_qmom_source()
{
    return [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.025); };
}

[[nodiscard]] qmomkit::MomentSource arithmetic_qmom_source()
{
    return [](const double time, const std::vector<double>& state)
    {
        std::vector<double> rates(state.size(), 0.0);
        for (std::size_t moment = 0; moment < state.size(); ++moment)
        {
            rates[moment] = 0.1 + 0.01 * time + 0.02 * static_cast<double>(moment) + 0.005 * state[moment];
        }
        return rates;
    };
}

[[nodiscard]] qmomkit::CqmomMomentSource free_cqmom_source()
{
    return [](double, const qmomkit::CqmomSourceContext& context) { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.015)}; };
}

[[nodiscard]] qmomkit::CqmomMomentSource arithmetic_cqmom_source()
{
    return [](const double time, const qmomkit::CqmomSourceContext& context)
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
}

[[nodiscard]] std::vector<double> flatten_qmom_result(const qmomkit::QmomResult& result)
{
    std::vector<double> values = result.weights;
    values.insert(values.end(), result.abscissas.begin(), result.abscissas.end());
    return values;
}

[[nodiscard]] std::vector<double> flatten_cqmom_result(const qmomkit::CqmomResult& result)
{
    std::vector<double> values = result.quadrature.weights();
    const std::vector<std::vector<double>>& abscissas = result.quadrature.abscissas();
    for (const std::vector<double>& dimension_abscissas : abscissas)
    {
        values.insert(values.end(), dimension_abscissas.begin(), dimension_abscissas.end());
    }
    return values;
}

void validate_qmom_inversion_fixture(const std::vector<double>& local_moments, const std::size_t moment_count, const std::size_t node_count,
                                     const qmomkit::QmomOptions& options)
{
    for (std::size_t q = 0; q * moment_count < local_moments.size(); ++q)
    {
        const qmomkit::QmomResult result =
            qmomkit::compute_qmom(qmomkit::MomentSet(local_moment_slice(local_moments, q, moment_count)), node_count, options);
        if (!result.ok())
        {
            throw std::runtime_error("QMOM inversion fixture failed: " + result.message);
        }
    }
}

void validate_cqmom_inversion_fixture(const std::vector<double>& local_moments, const qmomkit::CqmomLeanCase& lean_case)
{
    const std::size_t moment_count = lean_case.moments.values.size();
    for (std::size_t q = 0; q * moment_count < local_moments.size(); ++q)
    {
        qmomkit::CqmomLeanCase local_case = lean_case;
        local_case.moments.values = local_moment_slice(local_moments, q, moment_count);
        const qmomkit::CqmomResult result = qmomkit::compute_cqmom(local_case);
        if (!result.ok())
        {
            throw std::runtime_error("CQMOM inversion fixture failed: " + result.message);
        }
    }
}

[[nodiscard]] BenchmarkCase make_case(CaseMetadata metadata, std::function<void()> validate, std::function<std::vector<double>()> work)
{
    return BenchmarkCase{std::move(metadata), std::move(validate), std::move(work)};
}

void add_interpolation_projection_rows(std::vector<BenchmarkCase>& cases, const CaseMetadata& base, const Element& element,
                                       const std::vector<double>& coefficients, const std::vector<double>& local_rates,
                                       const std::vector<std::size_t>& batch_sizes)
{
    for (const std::size_t batch_size : batch_sizes)
    {
        CaseMetadata interpolation = base;
        interpolation.category = "interpolation_only";
        interpolation.element_batch_size = batch_size;
        interpolation.local_calls_per_element = 0;
        cases.push_back(
            make_case(interpolation, {},
                      [=] { return repeat_batch(batch_size, [&] { return interpolate_element_state(element, coefficients, base.moment_count); }); }));

        CaseMetadata projection = base;
        projection.category = "projection_inverse_mass_only";
        projection.element_batch_size = batch_size;
        projection.local_calls_per_element = 0;
        cases.push_back(
            make_case(projection, {}, [=]
                      { return repeat_batch(batch_size, [&] { return project_rates_to_coefficients(element, local_rates, base.moment_count); }); }));
    }
}

void add_qmom_continuous_case(std::vector<BenchmarkCase>& cases, const Element& element, const std::vector<std::size_t>& batch_sizes,
                              const std::string& case_name, const std::vector<double>& moments, const std::size_t node_count,
                              const qmomkit::PhaseSpaceVelocity& velocity, const std::string& callback_label)
{
    const std::size_t moment_count = moments.size();
    const std::vector<double> coefficients = constant_coefficients(element, moments);
    const std::vector<double> local_moments = interpolate_element_state(element, coefficients, moment_count);
    const ElementRhsResult probe = evaluate_qmom_continuous_rhs(element, coefficients, moment_count, node_count, 0.0, velocity);
    if (!probe.ok)
    {
        throw std::runtime_error("benchmark setup failed: " + describe_failure(probe));
    }
    const std::vector<double> local_rates = interpolate_element_state(element, probe.coefficient_rates, moment_count);

    const CaseMetadata base{case_name,
                            "",
                            element.basis_size(),
                            element.quadrature_count(),
                            moment_count,
                            1,
                            std::to_string(node_count) + "_qmom_nodes",
                            "not_applicable",
                            "not_applicable",
                            callback_label,
                            1,
                            element.quadrature_count(),
                            3 * element.quadrature_count()};
    add_interpolation_projection_rows(cases, base, element, coefficients, local_rates, batch_sizes);

    for (const std::size_t batch_size : batch_sizes)
    {
        CaseMetadata inversion = base;
        inversion.category = "qmom_inversion_only";
        inversion.element_batch_size = batch_size;
        cases.push_back(make_case(
            inversion, [=] { validate_qmom_inversion_fixture(local_moments, moment_count, node_count, {}); },
            [=]
            {
                return repeat_batch(batch_size,
                                    [&]
                                    {
                                        std::vector<double> values;
                                        for (std::size_t q = 0; q < element.quadrature_count(); ++q)
                                        {
                                            const qmomkit::QmomResult result = qmomkit::compute_qmom(
                                                qmomkit::MomentSet(local_moment_slice(local_moments, q, moment_count)), node_count);
                                            const std::vector<double> flattened = flatten_qmom_result(result);
                                            values.insert(values.end(), flattened.begin(), flattened.end());
                                        }
                                        return values;
                                    });
            }));

        CaseMetadata local = base;
        local.category = "qmom_continuous_local_rate";
        local.element_batch_size = batch_size;
        cases.push_back(make_case(
            local,
            [=]
            {
                const ElementRhsResult result = evaluate_qmom_continuous_rhs(element, coefficients, moment_count, node_count, 0.0, velocity);
                if (!result.ok)
                {
                    throw std::runtime_error("benchmark setup failed: " + describe_failure(result));
                }
            },
            [=]
            {
                return repeat_batch(batch_size,
                                    [&]
                                    {
                                        std::vector<double> rates(element.quadrature_count() * moment_count, 0.0);
                                        for (std::size_t q = 0; q < element.quadrature_count(); ++q)
                                        {
                                            const qmomkit::QmomContinuousMomentRateResult result = qmomkit::evaluate_qmom_continuous_moment_rates(
                                                local_moment_slice(local_moments, q, moment_count), node_count, 0.0, velocity);
                                            for (std::size_t moment = 0; moment < moment_count; ++moment)
                                            {
                                                rates[q * moment_count + moment] = result.rates[moment];
                                            }
                                        }
                                        return rates;
                                    });
            }));

        CaseMetadata rhs = base;
        rhs.category = "complete_element_rhs";
        rhs.element_batch_size = batch_size;
        cases.push_back(make_case(
            rhs,
            [=]
            {
                const ElementRhsResult result = evaluate_qmom_continuous_rhs(element, coefficients, moment_count, node_count, 0.0, velocity);
                if (!result.ok)
                {
                    throw std::runtime_error("benchmark setup failed: " + describe_failure(result));
                }
            },
            [=]
            {
                return repeat_batch(batch_size,
                                    [&]
                                    {
                                        const ElementRhsResult result =
                                            evaluate_qmom_continuous_rhs(element, coefficients, moment_count, node_count, 0.0, velocity);
                                        return result.coefficient_rates;
                                    });
            }));

        CaseMetadata rk = base;
        rk.category = "complete_external_ssp_rk3_step";
        rk.element_batch_size = batch_size;
        rk.local_calls_per_element = 3 * element.quadrature_count();
        cases.push_back(make_case(
            rk,
            [=]
            {
                const ElementRhsResult result =
                    ssp_rk3_coefficient_step(coefficients, 0.0, 0.01, [&](const double time, const std::vector<double>& state)
                                             { return evaluate_qmom_continuous_rhs(element, state, moment_count, node_count, time, velocity); });
                if (!result.ok)
                {
                    throw std::runtime_error("benchmark setup failed: " + describe_failure(result));
                }
            },
            [=]
            {
                return repeat_batch(batch_size,
                                    [&]
                                    {
                                        const ElementRhsResult result = ssp_rk3_coefficient_step(
                                            coefficients, 0.0, 0.01, [&](const double time, const std::vector<double>& state)
                                            { return evaluate_qmom_continuous_rhs(element, state, moment_count, node_count, time, velocity); });
                                        return result.coefficient_rates;
                                    });
            }));
    }
}

void add_qmom_discontinuous_case(std::vector<BenchmarkCase>& cases, const Element& element, const std::vector<std::size_t>& batch_sizes,
                                 const std::string& case_name, const qmomkit::MomentSource& source, const std::string& callback_label)
{
    const std::size_t moment_count = 4;
    const std::vector<double> moments = qmom_moments_from_terminals({0.4, 0.6}, {-0.5, 1.2}, moment_count);
    const std::vector<double> coefficients = constant_coefficients(element, moments);
    const std::vector<double> local_moments = interpolate_element_state(element, coefficients, moment_count);
    const ElementRhsResult probe = evaluate_qmom_discontinuous_rhs(element, coefficients, moment_count, 0.0, source);
    if (!probe.ok)
    {
        throw std::runtime_error("benchmark setup failed: " + describe_failure(probe));
    }
    const std::vector<double> local_rates = interpolate_element_state(element, probe.coefficient_rates, moment_count);

    const CaseMetadata base{case_name,
                            "",
                            element.basis_size(),
                            element.quadrature_count(),
                            moment_count,
                            1,
                            "2_qmom_nodes",
                            "not_applicable",
                            "not_applicable",
                            callback_label,
                            1,
                            element.quadrature_count(),
                            3 * element.quadrature_count()};
    add_interpolation_projection_rows(cases, base, element, coefficients, local_rates, batch_sizes);

    for (const std::size_t batch_size : batch_sizes)
    {
        CaseMetadata local = base;
        local.category = "qmom_discontinuous_local_rate";
        local.element_batch_size = batch_size;
        cases.push_back(make_case(local, {},
                                  [=]
                                  {
                                      return repeat_batch(batch_size,
                                                          [&]
                                                          {
                                                              std::vector<double> rates(element.quadrature_count() * moment_count, 0.0);
                                                              for (std::size_t q = 0; q < element.quadrature_count(); ++q)
                                                              {
                                                                  const qmomkit::QmomDiscontinuousMomentRateResult result =
                                                                      qmomkit::evaluate_qmom_discontinuous_moment_rates(
                                                                          0.0, local_moment_slice(local_moments, q, moment_count), source);
                                                                  for (std::size_t moment = 0; moment < moment_count; ++moment)
                                                                  {
                                                                      rates[q * moment_count + moment] = result.rates[moment];
                                                                  }
                                                              }
                                                              return rates;
                                                          });
                                  }));

        CaseMetadata rhs = base;
        rhs.category = "complete_element_rhs";
        rhs.element_batch_size = batch_size;
        cases.push_back(make_case(rhs, {},
                                  [=]
                                  {
                                      return repeat_batch(batch_size,
                                                          [&]
                                                          {
                                                              const ElementRhsResult result =
                                                                  evaluate_qmom_discontinuous_rhs(element, coefficients, moment_count, 0.0, source);
                                                              return result.coefficient_rates;
                                                          });
                                  }));

        CaseMetadata rk = base;
        rk.category = "complete_external_ssp_rk3_step";
        rk.element_batch_size = batch_size;
        rk.local_calls_per_element = 3 * element.quadrature_count();
        cases.push_back(make_case(
            rk,
            [=]
            {
                const ElementRhsResult result =
                    ssp_rk3_coefficient_step(coefficients, 0.0, 0.01, [&](const double time, const std::vector<double>& state)
                                             { return evaluate_qmom_discontinuous_rhs(element, state, moment_count, time, source); });
                if (!result.ok)
                {
                    throw std::runtime_error("benchmark setup failed: " + describe_failure(result));
                }
            },
            [=]
            {
                return repeat_batch(batch_size,
                                    [&]
                                    {
                                        const ElementRhsResult result = ssp_rk3_coefficient_step(
                                            coefficients, 0.0, 0.01, [&](const double time, const std::vector<double>& state)
                                            { return evaluate_qmom_discontinuous_rhs(element, state, moment_count, time, source); });
                                        return result.coefficient_rates;
                                    });
            }));
    }
}

[[nodiscard]] qmomkit::CqmomLeanCase make_cqmom_case(const qmomkit::PhaseSpaceHierarchy& hierarchy, const std::vector<double>& weights,
                                                     const std::vector<std::vector<double>>& abscissas,
                                                     const qmomkit::CqmomDiagnosticLevel diagnostic_level, const bool adaptive_branch_counts = false)
{
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::QmomOptions qmom_options;
    if (adaptive_branch_counts)
    {
        qmom_options = {true, 1.0e-10, 1.0e-10, 1.0e-12};
    }
    qmomkit::CqmomLeanCase lean_case{cqmom_moment_data_from_terminals(layout, weights, abscissas),
                                     {},
                                     qmom_options,
                                     adaptive_branch_counts,
                                     diagnostic_level,
                                     qmomkit::VandermondeSolver::automatic};
    return lean_case;
}

void add_cqmom_continuous_case(std::vector<BenchmarkCase>& cases, const Element& element, const std::vector<std::size_t>& batch_sizes,
                               const std::string& case_name, const qmomkit::CqmomLeanCase& lean_case, const std::string& node_configuration,
                               const std::string& conditioning_order, const qmomkit::PhaseSpaceVelocity& velocity, const std::string& callback_label)
{
    const std::size_t moment_count = lean_case.moments.values.size();
    const std::vector<double> coefficients = constant_coefficients(element, lean_case.moments.values);
    const std::vector<double> local_moments = interpolate_element_state(element, coefficients, moment_count);
    const ElementRhsResult probe = evaluate_cqmom_continuous_rhs(element, coefficients, lean_case, 0.0, velocity);
    if (!probe.ok)
    {
        throw std::runtime_error("benchmark setup failed: " + describe_failure(probe));
    }
    const std::vector<double> local_rates = interpolate_element_state(element, probe.coefficient_rates, moment_count);
    const std::string diagnostic_label = lean_case.diagnostic_level == qmomkit::CqmomDiagnosticLevel::none ? "none" : "detailed";

    const CaseMetadata base{case_name,
                            "",
                            element.basis_size(),
                            element.quadrature_count(),
                            moment_count,
                            2,
                            node_configuration,
                            conditioning_order,
                            diagnostic_label,
                            callback_label,
                            1,
                            element.quadrature_count(),
                            3 * element.quadrature_count()};
    add_interpolation_projection_rows(cases, base, element, coefficients, local_rates, batch_sizes);

    for (const std::size_t batch_size : batch_sizes)
    {
        CaseMetadata inversion = base;
        inversion.category = "cqmom_inversion_only";
        inversion.element_batch_size = batch_size;
        cases.push_back(make_case(
            inversion, [=] { validate_cqmom_inversion_fixture(local_moments, lean_case); },
            [=]
            {
                return repeat_batch(batch_size,
                                    [&]
                                    {
                                        std::vector<double> values;
                                        for (std::size_t q = 0; q < element.quadrature_count(); ++q)
                                        {
                                            qmomkit::CqmomLeanCase local_case = lean_case;
                                            local_case.moments.values = local_moment_slice(local_moments, q, moment_count);
                                            const qmomkit::CqmomResult result = qmomkit::compute_cqmom(local_case);
                                            const std::vector<double> flattened = flatten_cqmom_result(result);
                                            values.insert(values.end(), flattened.begin(), flattened.end());
                                        }
                                        return values;
                                    });
            }));

        CaseMetadata local = base;
        local.category = "cqmom_continuous_local_rate";
        local.element_batch_size = batch_size;
        cases.push_back(make_case(
            local,
            [=]
            {
                const ElementRhsResult result = evaluate_cqmom_continuous_rhs(element, coefficients, lean_case, 0.0, velocity);
                if (!result.ok)
                {
                    throw std::runtime_error("benchmark setup failed: " + describe_failure(result));
                }
            },
            [=]
            {
                return repeat_batch(batch_size,
                                    [&]
                                    {
                                        std::vector<double> rates(element.quadrature_count() * moment_count, 0.0);
                                        for (std::size_t q = 0; q < element.quadrature_count(); ++q)
                                        {
                                            qmomkit::CqmomLeanCase local_case = lean_case;
                                            local_case.moments.values = local_moment_slice(local_moments, q, moment_count);
                                            const qmomkit::CqmomContinuousMomentRateResult result =
                                                qmomkit::evaluate_cqmom_continuous_moment_rates(local_case, 0.0, velocity);
                                            for (std::size_t moment = 0; moment < moment_count; ++moment)
                                            {
                                                rates[q * moment_count + moment] = result.rates->values[moment];
                                            }
                                        }
                                        return rates;
                                    });
            }));

        CaseMetadata rhs = base;
        rhs.category = "complete_element_rhs";
        rhs.element_batch_size = batch_size;
        cases.push_back(make_case(
            rhs,
            [=]
            {
                const ElementRhsResult result = evaluate_cqmom_continuous_rhs(element, coefficients, lean_case, 0.0, velocity);
                if (!result.ok)
                {
                    throw std::runtime_error("benchmark setup failed: " + describe_failure(result));
                }
            },
            [=]
            {
                return repeat_batch(batch_size,
                                    [&]
                                    {
                                        const ElementRhsResult result =
                                            evaluate_cqmom_continuous_rhs(element, coefficients, lean_case, 0.0, velocity);
                                        return result.coefficient_rates;
                                    });
            }));
    }
}

void add_cqmom_discontinuous_case(std::vector<BenchmarkCase>& cases, const Element& element, const std::vector<std::size_t>& batch_sizes,
                                  const std::string& case_name, const qmomkit::CqmomLeanCase& lean_case, const std::string& node_configuration,
                                  const std::string& conditioning_order, const qmomkit::CqmomMomentSource& source, const std::string& callback_label)
{
    const std::size_t moment_count = lean_case.moments.values.size();
    const std::vector<double> coefficients = constant_coefficients(element, lean_case.moments.values);
    const std::vector<double> local_moments = interpolate_element_state(element, coefficients, moment_count);
    const ElementRhsResult probe = evaluate_cqmom_discontinuous_rhs(element, coefficients, lean_case.moments, 0.0, source);
    if (!probe.ok)
    {
        throw std::runtime_error("benchmark setup failed: " + describe_failure(probe));
    }
    const std::vector<double> local_rates = interpolate_element_state(element, probe.coefficient_rates, moment_count);

    const CaseMetadata base{case_name,
                            "",
                            element.basis_size(),
                            element.quadrature_count(),
                            moment_count,
                            2,
                            node_configuration,
                            conditioning_order,
                            "not_applicable",
                            callback_label,
                            1,
                            element.quadrature_count(),
                            3 * element.quadrature_count()};
    add_interpolation_projection_rows(cases, base, element, coefficients, local_rates, batch_sizes);

    for (const std::size_t batch_size : batch_sizes)
    {
        CaseMetadata local = base;
        local.category = "cqmom_discontinuous_local_rate";
        local.element_batch_size = batch_size;
        cases.push_back(make_case(local, {},
                                  [=]
                                  {
                                      return repeat_batch(batch_size,
                                                          [&]
                                                          {
                                                              std::vector<double> rates(element.quadrature_count() * moment_count, 0.0);
                                                              for (std::size_t q = 0; q < element.quadrature_count(); ++q)
                                                              {
                                                                  qmomkit::CqmomMomentData local_data = lean_case.moments;
                                                                  local_data.values = local_moment_slice(local_moments, q, moment_count);
                                                                  const qmomkit::CqmomDiscontinuousMomentRateResult result =
                                                                      qmomkit::evaluate_cqmom_discontinuous_moment_rates(0.0, local_data, source);
                                                                  for (std::size_t moment = 0; moment < moment_count; ++moment)
                                                                  {
                                                                      rates[q * moment_count + moment] = result.rates->values[moment];
                                                                  }
                                                              }
                                                              return rates;
                                                          });
                                  }));

        CaseMetadata rhs = base;
        rhs.category = "complete_element_rhs";
        rhs.element_batch_size = batch_size;
        cases.push_back(make_case(rhs, {},
                                  [=]
                                  {
                                      return repeat_batch(batch_size,
                                                          [&]
                                                          {
                                                              const ElementRhsResult result = evaluate_cqmom_discontinuous_rhs(
                                                                  element, coefficients, lean_case.moments, 0.0, source);
                                                              return result.coefficient_rates;
                                                          });
                                  }));

        CaseMetadata rk = base;
        rk.category = "complete_external_ssp_rk3_step";
        rk.element_batch_size = batch_size;
        rk.local_calls_per_element = 3 * element.quadrature_count();
        cases.push_back(make_case(
            rk,
            [=]
            {
                const ElementRhsResult result =
                    ssp_rk3_coefficient_step(coefficients, 0.0, 0.01, [&](const double time, const std::vector<double>& state)
                                             { return evaluate_cqmom_discontinuous_rhs(element, state, lean_case.moments, time, source); });
                if (!result.ok)
                {
                    throw std::runtime_error("benchmark setup failed: " + describe_failure(result));
                }
            },
            [=]
            {
                return repeat_batch(batch_size,
                                    [&]
                                    {
                                        const ElementRhsResult result = ssp_rk3_coefficient_step(
                                            coefficients, 0.0, 0.01, [&](const double time, const std::vector<double>& state)
                                            { return evaluate_cqmom_discontinuous_rhs(element, state, lean_case.moments, time, source); });
                                        return result.coefficient_rates;
                                    });
            }));
    }
}

[[nodiscard]] std::vector<BenchmarkCase> make_benchmark_cases(const BenchmarkOptions& options)
{
    const Element quadratic = make_lagrange_element(2);
    const Element cubic = make_lagrange_element(3);

    const qmomkit::CqmomLeanCase cqmom_identity_none =
        make_cqmom_case(qmomkit::PhaseSpaceHierarchy({2, 2}), {0.2, 0.3, 0.1, 0.4}, {{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}},
                        qmomkit::CqmomDiagnosticLevel::none);
    const qmomkit::CqmomLeanCase cqmom_identity_detailed =
        make_cqmom_case(qmomkit::PhaseSpaceHierarchy({2, 2}), {0.2, 0.3, 0.1, 0.4}, {{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}},
                        qmomkit::CqmomDiagnosticLevel::detailed);
    const qmomkit::CqmomLeanCase cqmom_nonidentity_detailed =
        make_cqmom_case(qmomkit::PhaseSpaceHierarchy({2, 2}, {1, 0}), {0.2, 0.3, 0.1, 0.4}, {{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}},
                        qmomkit::CqmomDiagnosticLevel::detailed);
    const qmomkit::CqmomLeanCase cqmom_larger_fixed =
        make_cqmom_case(qmomkit::PhaseSpaceHierarchy({3, 2}), {0.08, 0.12, 0.12, 0.18, 0.20, 0.30},
                        {{-1.0, -1.0, 0.0, 0.0, 1.0, 1.0}, {-0.5, 1.5, -0.5, 1.5, -0.5, 1.5}}, qmomkit::CqmomDiagnosticLevel::none);
    const qmomkit::CqmomLeanCase cqmom_adaptive = make_cqmom_case(qmomkit::PhaseSpaceHierarchy({2, 2}), {0.1, 0.3, 0.6},
                                                                  {{-1.0, -1.0, 1.0}, {-2.0, 2.0, 3.0}}, qmomkit::CqmomDiagnosticLevel::none, true);

    std::vector<BenchmarkCase> cases;
    add_qmom_continuous_case(cases, quadratic, options.batch_sizes, "qmom_1node_2mom_zero_velocity", qmom_moments_from_terminals({2.0}, {3.0}, 2), 1,
                             zero_velocity(), "zero_velocity");
    add_qmom_continuous_case(cases, quadratic, options.batch_sizes, "qmom_2node_4mom_affine_velocity",
                             qmom_moments_from_terminals({0.4, 0.6}, {-0.5, 1.2}, 4), 2, affine_velocity(), "cheap_affine_velocity");
    add_qmom_continuous_case(cases, cubic, options.batch_sizes, "qmom_2node_6mom_nonlinear_velocity",
                             qmom_moments_from_terminals({0.4, 0.6}, {-0.5, 1.2}, 6), 2, nonlinear_velocity(), "modest_nonlinear_velocity");

    add_qmom_discontinuous_case(cases, quadratic, options.batch_sizes, "qmom_2node_4mom_free_source", free_qmom_source(), "nearly_free_source");
    add_qmom_discontinuous_case(cases, cubic, options.batch_sizes, "qmom_2node_4mom_arithmetic_source", arithmetic_qmom_source(),
                                "modest_arithmetic_source");

    add_cqmom_continuous_case(cases, quadratic, options.batch_sizes, "cqmom_2x2_identity_none_affine_velocity", cqmom_identity_none, "2x2_fixed",
                              "identity", affine_velocity(), "cheap_affine_velocity");
    add_cqmom_continuous_case(cases, quadratic, options.batch_sizes, "cqmom_2x2_identity_detailed_affine_velocity", cqmom_identity_detailed,
                              "2x2_fixed", "identity", affine_velocity(), "cheap_affine_velocity");
    add_cqmom_continuous_case(cases, quadratic, options.batch_sizes, "cqmom_2x2_nonidentity_detailed_nonlinear_velocity", cqmom_nonidentity_detailed,
                              "2x2_fixed", "nonidentity", nonlinear_velocity(), "modest_nonlinear_velocity");
    add_cqmom_continuous_case(cases, cubic, options.batch_sizes, "cqmom_3x2_identity_none_affine_velocity", cqmom_larger_fixed, "3x2_fixed",
                              "identity", affine_velocity(), "cheap_affine_velocity");
    add_cqmom_continuous_case(cases, quadratic, options.batch_sizes, "cqmom_2x2_adaptive_identity_none_zero_velocity", cqmom_adaptive,
                              "2x2_adaptive_branch_counts", "identity", zero_velocity(), "zero_velocity");

    add_cqmom_discontinuous_case(cases, quadratic, options.batch_sizes, "cqmom_2x2_nonidentity_free_source", cqmom_nonidentity_detailed, "2x2_fixed",
                                 "nonidentity", free_cqmom_source(), "nearly_free_source");
    add_cqmom_discontinuous_case(cases, cubic, options.batch_sizes, "cqmom_2x2_nonidentity_arithmetic_source", cqmom_nonidentity_detailed,
                                 "2x2_fixed", "nonidentity", arithmetic_cqmom_source(), "modest_arithmetic_source");
    return cases;
}

void print_usage()
{
    std::cerr << "Usage: moment_rate_element_operator_benchmark [--short] [--csv] [--samples N] [--target-ms MS] "
                 "[--warmup N] [--min-iterations N] [--max-iterations N] [--filter TEXT]\n";
}

[[nodiscard]] BenchmarkOptions parse_options(const int argc, char** argv)
{
    BenchmarkOptions options;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--short")
        {
            options.short_run = true;
            options.samples = 3;
            options.warmup_iterations = 1;
            options.target_sample_seconds = 0.001;
            options.batch_sizes = {1, 10};
        }
        else if (arg == "--csv")
        {
            options.csv = true;
        }
        else if (arg == "--samples" && i + 1 < argc)
        {
            options.samples = static_cast<std::size_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--target-ms" && i + 1 < argc)
        {
            options.target_sample_seconds = std::stod(argv[++i]) * 1.0e-3;
        }
        else if (arg == "--warmup" && i + 1 < argc)
        {
            options.warmup_iterations = static_cast<std::size_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--min-iterations" && i + 1 < argc)
        {
            options.min_iterations = static_cast<std::size_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--max-iterations" && i + 1 < argc)
        {
            options.max_iterations = static_cast<std::size_t>(std::stoull(argv[++i]));
        }
        else if (arg == "--filter" && i + 1 < argc)
        {
            options.filter = argv[++i];
        }
        else
        {
            print_usage();
            throw std::invalid_argument("Unknown benchmark option: " + arg);
        }
    }
    if (options.samples == 0)
    {
        throw std::invalid_argument("--samples must be positive.");
    }
    if (!(options.target_sample_seconds > 0.0))
    {
        throw std::invalid_argument("--target-ms must be positive.");
    }
    if (options.max_iterations < options.min_iterations || options.max_iterations == 0)
    {
        throw std::invalid_argument("--max-iterations must be positive and at least --min-iterations.");
    }
    return options;
}

[[nodiscard]] bool matches_filter(const BenchmarkCase& benchmark_case, const std::string& filter)
{
    if (filter.empty())
    {
        return true;
    }
    const CaseMetadata& metadata = benchmark_case.metadata;
    return metadata.case_name.find(filter) != std::string::npos || metadata.category.find(filter) != std::string::npos ||
           metadata.node_configuration.find(filter) != std::string::npos || metadata.callback_label.find(filter) != std::string::npos ||
           metadata.diagnostic_level.find(filter) != std::string::npos;
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const BenchmarkOptions options = parse_options(argc, argv);
        const std::vector<BenchmarkCase> cases = make_benchmark_cases(options);
        if (options.csv)
        {
            print_csv_header();
        }
        else
        {
            print_human_metadata(options);
        }

        std::size_t measured_cases = 0;
        for (const BenchmarkCase& benchmark_case : cases)
        {
            if (!matches_filter(benchmark_case, options.filter))
            {
                continue;
            }
            const TimingResult timing = measure_case(options, benchmark_case);
            if (options.csv)
            {
                print_csv_row(options, benchmark_case.metadata, timing);
            }
            else
            {
                print_human_row(benchmark_case.metadata, timing);
            }
            ++measured_cases;
        }

        if (measured_cases == 0)
        {
            throw std::runtime_error("No benchmark cases matched the requested filter.");
        }
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
