#include "monte_carlo/comparison.hpp"
#include "monte_carlo/evolution.hpp"
#include "monte_carlo/explicit_distribution.hpp"
#include "monte_carlo/moments.hpp"
#include "monte_carlo/sampling.hpp"
#include "qmomkit/multi_index.hpp"
#include "test_checks.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace
{

using qmomkit::monte_carlo::ComparisonStatus;
using qmomkit::monte_carlo::ExplicitDiscreteDistribution;
using qmomkit::monte_carlo::SamplingOptions;

struct EvolutionComparisonRequest
{
    ExplicitDiscreteDistribution initial_distribution;
    qmomkit::PhaseSpaceVelocity velocity;
    double t0 = 0.0;
    double t1 = 0.0;
    qmomkit::AdaptiveRkOptions adaptive_options;
    SamplingOptions sampling_options;
    std::size_t maximum_order = 0;
};

struct EvolutionComparisonResult
{
    qmomkit::Quadrature1D initial_quadrature;
    qmomkit::QmomNodeAdvectionResult qmomkit_result;
    qmomkit::monte_carlo::ExplicitDiscreteSamplingResult sampling_result;
    qmomkit::monte_carlo::EnsembleEvolutionResult monte_carlo_result;
    qmomkit::monte_carlo::MomentComparisonResult comparison;
    double number_density = 0.0;
    std::size_t sample_count = 0;
    unsigned long long seed = 0;
    double final_time = 0.0;
    qmomkit::AdaptiveRkOptions adaptive_options;
    ComparisonStatus status = ComparisonStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == ComparisonStatus::success;
    }
};

struct MultiDimensionalEvolutionComparisonRequest
{
    ExplicitDiscreteDistribution initial_distribution;
    qmomkit::PhaseSpaceVelocity velocity;
    double t0 = 0.0;
    double t1 = 0.0;
    qmomkit::AdaptiveRkOptions adaptive_options;
    SamplingOptions sampling_options;
    std::vector<std::size_t> maximum_orders;
};

struct MultiDimensionalEvolutionComparisonResult
{
    qmomkit::monte_carlo::ExplicitQuadratureEvolutionResult quadrature_result;
    qmomkit::monte_carlo::ExplicitDiscreteSamplingResult sampling_result;
    qmomkit::monte_carlo::EnsembleEvolutionResult monte_carlo_result;
    qmomkit::monte_carlo::MomentComparisonResult comparison;
    double number_density = 0.0;
    std::size_t sample_count = 0;
    unsigned long long seed = 0;
    double final_time = 0.0;
    qmomkit::AdaptiveRkOptions adaptive_options;
    ComparisonStatus status = ComparisonStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == ComparisonStatus::success;
    }
};

struct SampleCountConvergenceRow
{
    std::size_t sample_count = 0;
    unsigned long long seed = 0;
    qmomkit::MultiIndex multi_index;
    std::size_t flat_index = 0;
    double reference_value = 0.0;
    double monte_carlo_value = 0.0;
    double absolute_error = 0.0;
    double standard_error = 0.0;
    double z_score = 0.0;
};

struct SampleCountConvergenceRequest
{
    ExplicitDiscreteDistribution initial_distribution;
    qmomkit::PhaseSpaceVelocity velocity;
    double t0 = 0.0;
    double t1 = 0.0;
    qmomkit::AdaptiveRkOptions adaptive_options;
    unsigned long long seed = 0;
    std::vector<std::size_t> sample_counts;
    std::vector<std::size_t> maximum_orders;
    std::vector<qmomkit::MultiIndex> selected_moments;
};

struct SampleCountConvergenceResult
{
    std::vector<SampleCountConvergenceRow> rows;
    double number_density = 0.0;
    double final_time = 0.0;
    unsigned long long seed = 0;
    qmomkit::AdaptiveRkOptions adaptive_options;
    ComparisonStatus status = ComparisonStatus::success;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == ComparisonStatus::success;
    }
};

qmomkit::AdaptiveRkOptions accurate_options()
{
    qmomkit::AdaptiveRkOptions options;
    options.relative_tolerance = 1.0e-10;
    options.absolute_tolerance = 1.0e-12;
    options.initial_step = 0.01;
    options.max_step = 0.05;
    return options;
}

void check_close(const double actual, const double expected, const double tolerance)
{
    CHECK_CLOSE(actual, expected, tolerance);
}

std::string row_context(const qmomkit::monte_carlo::MomentComparisonRow& row)
{
    std::ostringstream context;
    context << "flat_index=" << row.order << " reference_value=" << row.reference_value << " computed_value=" << row.computed_value
            << " absolute_error=" << row.absolute_error << " standard_error=" << row.standard_error << " z_score=" << row.z_score;
    return context.str();
}

void check_statistically_close(const qmomkit::monte_carlo::MomentComparisonRow& row)
{
    const double tolerance = std::max(6.0 * row.standard_error, 1.0e-12);
    CHECK_CLOSE_CONTEXT(row.computed_value, row.reference_value, tolerance, 0.0, row_context(row));
    CHECK_CONTEXT(std::abs(row.z_score) <= 6.0 || !std::isfinite(row.z_score), row_context(row));
}

void check_same_samples(const qmomkit::monte_carlo::Ensemble& lhs, const qmomkit::monte_carlo::Ensemble& rhs)
{
    CHECK(lhs.size() == rhs.size());
    for (std::size_t sample = 0; sample < lhs.size(); ++sample)
    {
        CHECK(lhs[sample] == rhs[sample]);
    }
}

const qmomkit::monte_carlo::MomentComparisonRow& find_row(const qmomkit::monte_carlo::MomentComparisonResult& result,
                                                          const qmomkit::MultiIndex& order)
{
    for (const qmomkit::monte_carlo::MomentComparisonRow& row : result.rows)
    {
        if (row.multi_index == order)
        {
            return row;
        }
    }
    CHECK_CONTEXT(false, "requested comparison row was not found");
    return result.rows.front();
}

bool all_finite(const std::vector<double>& values)
{
    return std::all_of(values.begin(), values.end(), [](const double value) { return std::isfinite(value); });
}

double integer_power(const double base, const std::size_t exponent)
{
    double result = 1.0;
    for (std::size_t i = 0; i < exponent; ++i)
    {
        result *= base;
    }
    return result;
}

double explicit_number_density(const ExplicitDiscreteDistribution& distribution)
{
    double number_density = 0.0;
    for (const double weight : distribution.weights)
    {
        number_density += weight;
    }
    return number_density;
}

bool validate_explicit_distribution(const ExplicitDiscreteDistribution& distribution, double& number_density, std::size_t& dimensions,
                                    std::string& message)
{
    if (distribution.nodes.empty() || distribution.nodes.size() != distribution.weights.size())
    {
        message = "Invalid explicit distribution.";
        return false;
    }
    dimensions = distribution.nodes.front().size();
    if (dimensions == 0)
    {
        message = "Invalid explicit distribution dimension.";
        return false;
    }
    for (const std::vector<double>& node : distribution.nodes)
    {
        if (node.size() != dimensions || !all_finite(node))
        {
            message = "Invalid explicit distribution node.";
            return false;
        }
    }
    if (!all_finite(distribution.weights))
    {
        message = "Invalid explicit distribution weights.";
        return false;
    }
    for (const double weight : distribution.weights)
    {
        if (weight < 0.0)
        {
            message = "Invalid explicit distribution weights.";
            return false;
        }
    }
    number_density = explicit_number_density(distribution);
    if (!(number_density > 0.0) || !std::isfinite(number_density))
    {
        message = "Invalid explicit distribution density.";
        return false;
    }
    if (distribution.number_density.has_value())
    {
        const double supplied = *distribution.number_density;
        const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(number_density));
        if (!std::isfinite(supplied) || supplied <= 0.0 || std::abs(supplied - number_density) > tolerance)
        {
            message = "Supplied number density must match weights.";
            return false;
        }
        number_density = supplied;
    }
    return true;
}

std::vector<std::size_t> extents_from_maximum_orders(const std::vector<std::size_t>& maximum_orders)
{
    std::vector<std::size_t> extents(maximum_orders.size(), 0);
    for (std::size_t dimension = 0; dimension < maximum_orders.size(); ++dimension)
    {
        if (maximum_orders[dimension] == std::numeric_limits<std::size_t>::max())
        {
            return {};
        }
        extents[dimension] = maximum_orders[dimension] + 1;
    }
    return extents;
}

double quadrature_moment(const qmomkit::Quadrature1D& quadrature, const std::size_t order)
{
    double value = 0.0;
    for (std::size_t node = 0; node < quadrature.abscissas.size(); ++node)
    {
        value += quadrature.weights[node] * integer_power(quadrature.abscissas[node], order);
    }
    return value;
}

double quadrature_number_density(const qmomkit::Quadrature1D& quadrature)
{
    double number_density = 0.0;
    for (const double weight : quadrature.weights)
    {
        number_density += weight;
    }
    return number_density;
}

double quadrature_moment(const ExplicitDiscreteDistribution& quadrature, const qmomkit::MultiIndex& order)
{
    double value = 0.0;
    for (std::size_t node = 0; node < quadrature.nodes.size(); ++node)
    {
        double monomial = 1.0;
        for (std::size_t dimension = 0; dimension < order.size(); ++dimension)
        {
            monomial *= integer_power(quadrature.nodes[node][dimension], order[dimension]);
        }
        value += quadrature.weights[node] * monomial;
    }
    return value;
}

qmomkit::MomentSet quadrature_moment_set(const qmomkit::Quadrature1D& quadrature, const std::size_t maximum_order)
{
    std::vector<double> values(maximum_order + 1, 0.0);
    for (std::size_t order = 0; order <= maximum_order; ++order)
    {
        values[order] = quadrature_moment(quadrature, order);
    }
    return qmomkit::MomentSet(std::move(values));
}

qmomkit::RawMomentTensor quadrature_moment_tensor(const ExplicitDiscreteDistribution& quadrature, const std::vector<std::size_t>& maximum_orders)
{
    const std::vector<std::size_t> extents = extents_from_maximum_orders(maximum_orders);
    if (extents.empty())
    {
        return {};
    }
    const std::vector<qmomkit::MultiIndex> orders = qmomkit::enumerate_multi_indices(extents);
    std::vector<double> values(orders.size(), 0.0);
    for (std::size_t flat = 0; flat < orders.size(); ++flat)
    {
        values[flat] = quadrature_moment(quadrature, orders[flat]);
    }
    return qmomkit::RawMomentTensor{extents, std::move(values)};
}

qmomkit::Quadrature1D make_one_dimensional_quadrature(const ExplicitDiscreteDistribution& distribution)
{
    std::vector<double> abscissas;
    abscissas.reserve(distribution.nodes.size());
    for (const std::vector<double>& node : distribution.nodes)
    {
        abscissas.push_back(node.front());
    }
    return qmomkit::Quadrature1D{std::move(abscissas), distribution.weights, qmomkit::QuadratureStatus::success, ""};
}

qmomkit::monte_carlo::MomentComparisonResult compare_quadrature_to_ensemble(const qmomkit::Quadrature1D& quadrature,
                                                                            const qmomkit::monte_carlo::Ensemble& ensemble,
                                                                            const double number_density, const std::size_t maximum_order)
{
    if (!quadrature.ok() || quadrature.abscissas.empty() || quadrature.abscissas.size() != quadrature.weights.size() ||
        !all_finite(quadrature.abscissas) || !all_finite(quadrature.weights))
    {
        return qmomkit::monte_carlo::MomentComparisonResult{{}, ComparisonStatus::invalid_input, "Invalid quadrature."};
    }
    for (const double weight : quadrature.weights)
    {
        if (weight < 0.0)
        {
            return qmomkit::monte_carlo::MomentComparisonResult{{}, ComparisonStatus::invalid_input, "Invalid quadrature weights."};
        }
    }
    const double quadrature_density = quadrature_number_density(quadrature);
    const double density_tolerance = 64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(quadrature_density));
    if (!std::isfinite(number_density) || number_density <= 0.0 || std::abs(quadrature_density - number_density) > density_tolerance)
    {
        return qmomkit::monte_carlo::MomentComparisonResult{{}, ComparisonStatus::invalid_input, "Invalid number density."};
    }
    const qmomkit::monte_carlo::OneDimensionalMomentEstimate estimate =
        qmomkit::monte_carlo::compute_ensemble_moments(ensemble, number_density, maximum_order);
    if (!estimate.ok())
    {
        return qmomkit::monte_carlo::MomentComparisonResult{{}, ComparisonStatus::moment_failed, estimate.message};
    }
    return qmomkit::monte_carlo::compare_moments(quadrature_moment_set(quadrature, maximum_order), estimate.moments, estimate.standard_errors);
}

qmomkit::monte_carlo::MomentComparisonResult compare_quadrature_to_ensemble(const ExplicitDiscreteDistribution& quadrature,
                                                                            const qmomkit::monte_carlo::Ensemble& ensemble,
                                                                            const double number_density,
                                                                            const std::vector<std::size_t>& maximum_orders)
{
    double validated_density = 0.0;
    std::size_t dimensions = 0;
    std::string message;
    if (!validate_explicit_distribution(quadrature, validated_density, dimensions, message) || maximum_orders.size() != dimensions)
    {
        return qmomkit::monte_carlo::MomentComparisonResult{{}, ComparisonStatus::invalid_input, std::move(message)};
    }
    const double density_tolerance = 64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(validated_density));
    if (!std::isfinite(number_density) || number_density <= 0.0 || std::abs(validated_density - number_density) > density_tolerance)
    {
        return qmomkit::monte_carlo::MomentComparisonResult{{}, ComparisonStatus::invalid_input, "Invalid number density."};
    }
    const qmomkit::monte_carlo::TensorMomentEstimate estimate =
        qmomkit::monte_carlo::compute_ensemble_moments(ensemble, number_density, maximum_orders);
    if (!estimate.ok())
    {
        return qmomkit::monte_carlo::MomentComparisonResult{{}, ComparisonStatus::moment_failed, estimate.message};
    }
    return qmomkit::monte_carlo::compare_moments(quadrature_moment_tensor(quadrature, maximum_orders), estimate.moments, estimate.standard_errors);
}

double local_z_score(const double error, const double standard_error)
{
    if (standard_error == 0.0)
    {
        return error == 0.0 ? 0.0 : std::copysign(std::numeric_limits<double>::infinity(), error);
    }
    return error / standard_error;
}

EvolutionComparisonResult failed_evolution_result(const EvolutionComparisonRequest& request, const ComparisonStatus status, std::string message)
{
    EvolutionComparisonResult result;
    result.status = status;
    result.message = std::move(message);
    result.sample_count = request.sampling_options.sample_count;
    result.seed = request.sampling_options.seed;
    result.final_time = request.t0;
    result.adaptive_options = request.adaptive_options;
    return result;
}

MultiDimensionalEvolutionComparisonResult failed_multidimensional_result(const MultiDimensionalEvolutionComparisonRequest& request,
                                                                         const ComparisonStatus status, std::string message)
{
    MultiDimensionalEvolutionComparisonResult result;
    result.status = status;
    result.message = std::move(message);
    result.sample_count = request.sampling_options.sample_count;
    result.seed = request.sampling_options.seed;
    result.final_time = request.t0;
    result.adaptive_options = request.adaptive_options;
    return result;
}

EvolutionComparisonResult compare_evolved_quadrature_to_ensemble(const EvolutionComparisonRequest& request)
{
    if (!request.velocity)
    {
        return failed_evolution_result(request, ComparisonStatus::invalid_input, "Phase-space velocity callback must be provided.");
    }
    if (request.sampling_options.sample_count == 0)
    {
        return failed_evolution_result(request, ComparisonStatus::invalid_input, "Sample count must be positive.");
    }
    if (!std::isfinite(request.t0) || !std::isfinite(request.t1) || request.t0 == request.t1)
    {
        return failed_evolution_result(request, ComparisonStatus::invalid_input, "Evolution time interval is invalid.");
    }

    EvolutionComparisonResult result;
    result.sample_count = request.sampling_options.sample_count;
    result.seed = request.sampling_options.seed;
    result.final_time = request.t1;
    result.adaptive_options = request.adaptive_options;

    double number_density = 0.0;
    std::size_t dimensions = 0;
    std::string message;
    if (!validate_explicit_distribution(request.initial_distribution, number_density, dimensions, message))
    {
        return failed_evolution_result(request, ComparisonStatus::invalid_input, std::move(message));
    }
    if (dimensions != 1)
    {
        return failed_evolution_result(request, ComparisonStatus::invalid_input,
                                       "One-dimensional comparison requires one-dimensional explicit nodes.");
    }
    result.number_density = number_density;

    result.initial_quadrature = make_one_dimensional_quadrature(request.initial_distribution);
    result.qmomkit_result =
        qmomkit::monte_carlo::evolve_quadrature_nodes(result.initial_quadrature, request.velocity, request.t0, request.t1, request.adaptive_options);
    if (!result.qmomkit_result.ok())
    {
        result.status = ComparisonStatus::qmomkit_failed;
        result.message = result.qmomkit_result.message;
        result.final_time = result.qmomkit_result.diagnostics.final_time;
        return result;
    }

    result.sampling_result = qmomkit::monte_carlo::sample_explicit_discrete_distribution(request.initial_distribution, request.sampling_options);
    if (!result.sampling_result.ok())
    {
        result.status = ComparisonStatus::sampling_failed;
        result.message = result.sampling_result.message;
        return result;
    }

    result.monte_carlo_result =
        qmomkit::monte_carlo::evolve_ensemble(result.sampling_result.ensemble, request.velocity, request.t0, request.t1, request.adaptive_options);
    if (!result.monte_carlo_result.ok())
    {
        result.status = ComparisonStatus::monte_carlo_failed;
        result.message = result.monte_carlo_result.message;
        result.final_time = result.monte_carlo_result.final_time;
        return result;
    }

    result.comparison =
        compare_quadrature_to_ensemble(result.qmomkit_result.quadrature, result.monte_carlo_result.ensemble, number_density, request.maximum_order);
    result.status = result.comparison.status;
    result.message = result.comparison.message;
    return result;
}

MultiDimensionalEvolutionComparisonResult compare_evolved_quadrature_to_ensemble(const MultiDimensionalEvolutionComparisonRequest& request)
{
    if (!request.velocity)
    {
        return failed_multidimensional_result(request, ComparisonStatus::invalid_input, "Phase-space velocity callback must be provided.");
    }
    if (request.sampling_options.sample_count == 0)
    {
        return failed_multidimensional_result(request, ComparisonStatus::invalid_input, "Sample count must be positive.");
    }
    if (!std::isfinite(request.t0) || !std::isfinite(request.t1) || request.t0 == request.t1)
    {
        return failed_multidimensional_result(request, ComparisonStatus::invalid_input, "Evolution time interval is invalid.");
    }

    MultiDimensionalEvolutionComparisonResult result;
    result.sample_count = request.sampling_options.sample_count;
    result.seed = request.sampling_options.seed;
    result.final_time = request.t1;
    result.adaptive_options = request.adaptive_options;

    double number_density = 0.0;
    std::size_t dimensions = 0;
    std::string message;
    if (!validate_explicit_distribution(request.initial_distribution, number_density, dimensions, message))
    {
        return failed_multidimensional_result(request, ComparisonStatus::invalid_input, std::move(message));
    }
    if (request.maximum_orders.size() != dimensions)
    {
        return failed_multidimensional_result(request, ComparisonStatus::invalid_input,
                                              "Maximum-order rank must match explicit quadrature dimension.");
    }
    result.number_density = number_density;

    result.quadrature_result = qmomkit::monte_carlo::evolve_explicit_quadrature(request.initial_distribution, request.velocity, request.t0,
                                                                                request.t1, request.adaptive_options);
    if (!result.quadrature_result.ok())
    {
        result.status = result.quadrature_result.status;
        result.message = result.quadrature_result.message;
        result.final_time = result.quadrature_result.final_time;
        return result;
    }

    result.sampling_result = qmomkit::monte_carlo::sample_explicit_discrete_distribution(request.initial_distribution, request.sampling_options);
    if (!result.sampling_result.ok())
    {
        result.status = ComparisonStatus::sampling_failed;
        result.message = result.sampling_result.message;
        return result;
    }

    result.monte_carlo_result =
        qmomkit::monte_carlo::evolve_ensemble(result.sampling_result.ensemble, request.velocity, request.t0, request.t1, request.adaptive_options);
    if (!result.monte_carlo_result.ok())
    {
        result.status = ComparisonStatus::monte_carlo_failed;
        result.message = result.monte_carlo_result.message;
        result.final_time = result.monte_carlo_result.final_time;
        return result;
    }

    result.comparison = compare_quadrature_to_ensemble(result.quadrature_result.evolved_distribution, result.monte_carlo_result.ensemble,
                                                       number_density, request.maximum_orders);
    result.status = result.comparison.status;
    result.message = result.comparison.message;
    return result;
}

SampleCountConvergenceResult run_sample_count_convergence(const SampleCountConvergenceRequest& request)
{
    SampleCountConvergenceResult result;
    result.final_time = request.t0;
    result.seed = request.seed;
    result.adaptive_options = request.adaptive_options;

    double number_density = 0.0;
    std::size_t dimensions = 0;
    std::string message;
    if (request.sample_counts.empty())
    {
        result.status = ComparisonStatus::invalid_input;
        result.message = "At least one sample count is required.";
        return result;
    }
    if (!validate_explicit_distribution(request.initial_distribution, number_density, dimensions, message))
    {
        result.status = ComparisonStatus::invalid_input;
        result.message = std::move(message);
        return result;
    }
    result.number_density = number_density;

    const std::vector<std::size_t> extents = extents_from_maximum_orders(request.maximum_orders);
    std::vector<qmomkit::MultiIndex> selected = request.selected_moments;
    if (selected.empty())
    {
        selected = qmomkit::enumerate_multi_indices(extents);
    }

    const qmomkit::monte_carlo::ExplicitQuadratureEvolutionResult reference_evolution = qmomkit::monte_carlo::evolve_explicit_quadrature(
        request.initial_distribution, request.velocity, request.t0, request.t1, request.adaptive_options);
    if (!reference_evolution.ok())
    {
        result.status = reference_evolution.status;
        result.message = reference_evolution.message;
        result.final_time = reference_evolution.final_time;
        return result;
    }
    result.final_time = request.t1;
    const qmomkit::RawMomentTensor reference = quadrature_moment_tensor(reference_evolution.evolved_distribution, request.maximum_orders);

    for (const std::size_t sample_count : request.sample_counts)
    {
        if (sample_count == 0)
        {
            result.status = ComparisonStatus::invalid_input;
            result.message = "Sample counts must be positive.";
            return result;
        }
        const SamplingOptions sampling_options{sample_count, request.seed};
        const qmomkit::monte_carlo::ExplicitDiscreteSamplingResult sampling_result =
            qmomkit::monte_carlo::sample_explicit_discrete_distribution(request.initial_distribution, sampling_options);
        if (!sampling_result.ok())
        {
            result.status = ComparisonStatus::sampling_failed;
            result.message = sampling_result.message;
            return result;
        }
        const qmomkit::monte_carlo::EnsembleEvolutionResult monte_carlo_result =
            qmomkit::monte_carlo::evolve_ensemble(sampling_result.ensemble, request.velocity, request.t0, request.t1, request.adaptive_options);
        if (!monte_carlo_result.ok())
        {
            result.status = ComparisonStatus::monte_carlo_failed;
            result.message = monte_carlo_result.message;
            result.final_time = monte_carlo_result.final_time;
            return result;
        }
        const qmomkit::monte_carlo::TensorMomentEstimate estimate =
            qmomkit::monte_carlo::compute_ensemble_moments(monte_carlo_result.ensemble, number_density, request.maximum_orders);
        if (!estimate.ok())
        {
            result.status = ComparisonStatus::moment_failed;
            result.message = estimate.message;
            return result;
        }

        for (const qmomkit::MultiIndex& order : selected)
        {
            const std::size_t flat = qmomkit::flatten_multi_index(order, extents);
            const double reference_value = reference.values[flat];
            const double monte_carlo_value = estimate.moments.values[flat];
            const double error = monte_carlo_value - reference_value;
            result.rows.push_back(SampleCountConvergenceRow{sample_count, request.seed, order, flat, reference_value, monte_carlo_value,
                                                            std::abs(error), estimate.standard_errors[flat],
                                                            local_z_score(error, estimate.standard_errors[flat])});
        }
    }

    result.status = ComparisonStatus::success;
    result.message.clear();
    return result;
}

void write_comparison_text(std::ostream& output, const EvolutionComparisonResult& result)
{
    qmomkit::monte_carlo::write_comparison_text(output, result.comparison,
                                                qmomkit::monte_carlo::ComparisonOutputMetadata{result.final_time, result.sample_count, result.seed,
                                                                                               result.number_density, result.adaptive_options});
}

void write_comparison_text(std::ostream& output, const MultiDimensionalEvolutionComparisonResult& result)
{
    qmomkit::monte_carlo::write_comparison_text(output, result.comparison,
                                                qmomkit::monte_carlo::ComparisonOutputMetadata{result.final_time, result.sample_count, result.seed,
                                                                                               result.number_density, result.adaptive_options});
}

void test_explicit_quadrature_construction()
{
    const ExplicitDiscreteDistribution distribution{{{-1.0}, {0.5}, {2.0}}, {0.4, 1.1, 0.7}, 2.2};
    const qmomkit::Quadrature1D quadrature{{-1.0, 0.5, 2.0}, distribution.weights, qmomkit::QuadratureStatus::success, ""};

    const qmomkit::monte_carlo::MomentComparisonResult direct = compare_quadrature_to_ensemble(quadrature, {{-1.0}, {0.5}, {2.0}}, 2.2, 0);
    CHECK(direct.ok());
    check_close(direct.rows.front().reference_value, 2.2, 1.0e-15);

    const qmomkit::monte_carlo::MomentComparisonResult mismatched_comparison =
        compare_quadrature_to_ensemble(quadrature, {{-1.0}, {0.5}, {2.0}}, 3.0, 0);
    CHECK(!mismatched_comparison.ok());

    const qmomkit::Quadrature1D negative_weight{{0.0}, {-1.0}, qmomkit::QuadratureStatus::success, ""};
    const qmomkit::monte_carlo::MomentComparisonResult negative_comparison = compare_quadrature_to_ensemble(negative_weight, {{0.0}}, 1.0, 0);
    CHECK(!negative_comparison.ok());

    const ExplicitDiscreteDistribution mismatched_density{{{0.0}}, {1.0}, 2.0};
    const qmomkit::monte_carlo::ExplicitDiscreteSamplingResult invalid =
        qmomkit::monte_carlo::sample_explicit_discrete_distribution(mismatched_density, SamplingOptions{10, 1});
    CHECK(!invalid.ok());
}

void test_direct_moment_comparison_utilities()
{
    const qmomkit::MomentSet reference_1d{2.0, 0.0, 8.0};
    const qmomkit::MomentSet same_1d{2.0, 0.0, 8.0};
    const qmomkit::monte_carlo::MomentComparisonResult identical = qmomkit::monte_carlo::compare_moments(reference_1d, same_1d, {0.0, 0.0, 0.0});
    CHECK(identical.ok());
    for (const qmomkit::monte_carlo::MomentComparisonRow& row : identical.rows)
    {
        check_close(row.signed_error, 0.0, 0.0);
        check_close(row.absolute_error, 0.0, 0.0);
        check_close(row.z_score, 0.0, 0.0);
    }

    const qmomkit::MomentSet shifted_1d{2.0, 0.2, 7.8};
    const qmomkit::monte_carlo::MomentComparisonResult shifted = qmomkit::monte_carlo::compare_moments(reference_1d, shifted_1d, {0.0, 0.1, 0.2});
    CHECK(shifted.ok());
    check_close(shifted.rows[1].signed_error, 0.2, 1.0e-15);
    check_close(shifted.rows[1].absolute_error, 0.2, 1.0e-15);
    check_close(shifted.rows[2].signed_error, -0.2, 1.0e-14);
    check_close(shifted.rows[2].absolute_error, 0.2, 1.0e-14);
    check_close(shifted.rows[1].z_score, 2.0, 1.0e-15);
    check_close(shifted.rows[2].z_score, -1.0, 1.0e-14);
    CHECK(std::isinf(shifted.rows[0].z_score) == false);

    const qmomkit::RawMomentTensor reference_tensor{{3, 4, 2}, std::vector<double>(24, 0.0)};
    qmomkit::RawMomentTensor monte_carlo_tensor = reference_tensor;
    monte_carlo_tensor.values[qmomkit::flatten_multi_index({2, 3, 1}, reference_tensor.extents)] = 0.25;
    std::vector<double> standard_errors(24, 0.0);
    standard_errors[qmomkit::flatten_multi_index({2, 3, 1}, reference_tensor.extents)] = 0.5;
    const qmomkit::monte_carlo::MomentComparisonResult tensor_result =
        qmomkit::monte_carlo::compare_moments(reference_tensor, monte_carlo_tensor, standard_errors);
    CHECK(tensor_result.ok());
    CHECK(tensor_result.rows.size() == 24);
    const std::size_t flat = qmomkit::flatten_multi_index({2, 3, 1}, reference_tensor.extents);
    CHECK(tensor_result.rows[flat].multi_index == qmomkit::MultiIndex({2, 3, 1}));
    check_close(tensor_result.rows[flat].signed_error, 0.25, 1.0e-15);
    check_close(tensor_result.rows[flat].absolute_error, 0.25, 1.0e-15);
    check_close(tensor_result.rows[flat].z_score, 0.5, 1.0e-15);

    const qmomkit::MomentSet zero_reference{0.0, 0.0};
    const qmomkit::MomentSet nonzero_with_zero_error{0.0, 1.0};
    const qmomkit::monte_carlo::MomentComparisonResult infinite_z =
        qmomkit::monte_carlo::compare_moments(zero_reference, nonzero_with_zero_error, {0.0, 0.0});
    CHECK(infinite_z.ok());
    check_close(infinite_z.rows[0].z_score, 0.0, 0.0);
    CHECK(std::isinf(infinite_z.rows[1].z_score));
    CHECK(infinite_z.rows[1].z_score > 0.0);
}

void test_multidimensional_sampling_and_dimension_validation()
{
    const ExplicitDiscreteDistribution distribution{{{-0.2, 0.1}, {0.4, -0.1}, {0.1, 0.3}}, {0.8, 0.9, 0.7}, {}};
    const SamplingOptions options{4000, 5101};
    const qmomkit::monte_carlo::ExplicitDiscreteSamplingResult first =
        qmomkit::monte_carlo::sample_explicit_discrete_distribution(distribution, options);
    const qmomkit::monte_carlo::ExplicitDiscreteSamplingResult second =
        qmomkit::monte_carlo::sample_explicit_discrete_distribution(distribution, options);
    CHECK(first.ok());
    CHECK(second.ok());
    check_close(first.number_density, 2.4, 1.0e-15);
    check_same_samples(first.ensemble, second.ensemble);
    for (const qmomkit::monte_carlo::Sample& sample : first.ensemble)
    {
        CHECK(sample.size() == 2);
    }

    const ExplicitDiscreteDistribution inconsistent{{{0.0, 1.0}, {2.0}}, {1.0, 1.0}, {}};
    const qmomkit::monte_carlo::ExplicitDiscreteSamplingResult invalid =
        qmomkit::monte_carlo::sample_explicit_discrete_distribution(inconsistent, options);
    CHECK(!invalid.ok());
}

void test_deterministic_sampling_and_number_density()
{
    const ExplicitDiscreteDistribution distribution{{{-1.0}, {0.5}, {2.0}}, {0.4, 1.1, 0.7}, {}};
    const SamplingOptions options{50000, 3001};
    const qmomkit::monte_carlo::ExplicitDiscreteSamplingResult first =
        qmomkit::monte_carlo::sample_explicit_discrete_distribution(distribution, options);
    const qmomkit::monte_carlo::ExplicitDiscreteSamplingResult second =
        qmomkit::monte_carlo::sample_explicit_discrete_distribution(distribution, options);

    CHECK(first.ok());
    CHECK(second.ok());
    check_close(first.number_density, 2.2, 1.0e-15);
    check_same_samples(first.ensemble, second.ensemble);

    const qmomkit::monte_carlo::OneDimensionalMomentEstimate moments =
        qmomkit::monte_carlo::compute_ensemble_moments(first.ensemble, first.number_density, 2);
    CHECK(moments.ok());
    check_close(moments.moments[0], first.number_density, 1.0e-15);
    check_close(moments.standard_errors[0], 0.0, 0.0);
}

void test_two_dimensional_constant_velocity_comparison()
{
    const ExplicitDiscreteDistribution distribution{{{-0.2, 0.1}, {0.25, -0.15}, {0.4, 0.3}}, {0.8, 1.0, 0.9}, {}};
    const double t0 = 0.0;
    const double t1 = 0.25;
    const double cx = 0.12;
    const double cy = -0.08;
    const qmomkit::PhaseSpaceVelocity callback = [cx, cy](double, const std::vector<double>& coordinates)
    {
        CHECK(coordinates.size() == 2);
        return std::vector<double>{cx, cy};
    };

    MultiDimensionalEvolutionComparisonRequest request;
    request.initial_distribution = distribution;
    request.velocity = callback;
    request.t0 = t0;
    request.t1 = t1;
    request.adaptive_options = accurate_options();
    request.sampling_options = SamplingOptions{7000, 5201};
    request.maximum_orders = {2, 2};

    const MultiDimensionalEvolutionComparisonResult result = compare_evolved_quadrature_to_ensemble(request);
    CHECK(result.ok());
    check_close(result.number_density, 2.7, 1.0e-15);
    CHECK(result.monte_carlo_result.ensemble.size() == request.sampling_options.sample_count);
    check_close(result.monte_carlo_result.final_time, t1, 0.0);

    for (std::size_t node = 0; node < distribution.nodes.size(); ++node)
    {
        check_close(result.quadrature_result.evolved_distribution.nodes[node][0], distribution.nodes[node][0] + cx * (t1 - t0), 1.0e-12);
        check_close(result.quadrature_result.evolved_distribution.nodes[node][1], distribution.nodes[node][1] + cy * (t1 - t0), 1.0e-12);
    }
    check_close(find_row(result.comparison, {0, 0}).reference_value, result.number_density, 1.0e-15);
    for (const qmomkit::monte_carlo::MomentComparisonRow& row : result.comparison.rows)
    {
        check_statistically_close(row);
    }
}

void test_two_dimensional_affine_comparison()
{
    const ExplicitDiscreteDistribution distribution{{{-0.15, 0.2}, {0.25, -0.1}, {0.45, 0.35}}, {0.7, 1.1, 0.8}, {}};
    const double t0 = -0.1;
    const double t1 = 0.2;
    const double ax = -0.2;
    const double bx = 0.35;
    const double ay = 0.15;
    const double by = -0.2;
    const qmomkit::PhaseSpaceVelocity callback = [ax, bx, ay, by](double, const std::vector<double>& coordinates)
    {
        CHECK(coordinates.size() == 2);
        return std::vector<double>{ax * coordinates[0] + bx, ay * coordinates[1] + by};
    };

    MultiDimensionalEvolutionComparisonRequest request;
    request.initial_distribution = distribution;
    request.velocity = callback;
    request.t0 = t0;
    request.t1 = t1;
    request.adaptive_options = accurate_options();
    request.sampling_options = SamplingOptions{8000, 5301};
    request.maximum_orders = {2, 2};

    const MultiDimensionalEvolutionComparisonResult result = compare_evolved_quadrature_to_ensemble(request);
    CHECK(result.ok());

    const double elapsed = t1 - t0;
    for (std::size_t node = 0; node < distribution.nodes.size(); ++node)
    {
        const double expected_x = (distribution.nodes[node][0] + bx / ax) * std::exp(ax * elapsed) - bx / ax;
        const double expected_y = (distribution.nodes[node][1] + by / ay) * std::exp(ay * elapsed) - by / ay;
        check_close(result.quadrature_result.evolved_distribution.nodes[node][0], expected_x, 1.0e-9);
        check_close(result.quadrature_result.evolved_distribution.nodes[node][1], expected_y, 1.0e-9);
    }
    check_statistically_close(find_row(result.comparison, {1, 1}));
    check_statistically_close(find_row(result.comparison, {2, 1}));
}

void test_two_dimensional_coupled_linear_comparison()
{
    const ExplicitDiscreteDistribution distribution{{{-0.25, 0.1}, {0.15, -0.2}, {0.35, 0.25}}, {0.9, 0.8, 1.0}, {}};
    const double t0 = 0.0;
    const double t1 = 0.2;
    const double a = -0.3;
    const double b = 0.4;
    const double c = -0.1;
    const qmomkit::PhaseSpaceVelocity callback = [a, b, c](double, const std::vector<double>& coordinates)
    {
        CHECK(coordinates.size() == 2);
        return std::vector<double>{a * coordinates[0], b * coordinates[0] + c * coordinates[1]};
    };

    MultiDimensionalEvolutionComparisonRequest request;
    request.initial_distribution = distribution;
    request.velocity = callback;
    request.t0 = t0;
    request.t1 = t1;
    request.adaptive_options = accurate_options();
    request.sampling_options = SamplingOptions{9000, 5401};
    request.maximum_orders = {2, 2};

    const MultiDimensionalEvolutionComparisonResult result = compare_evolved_quadrature_to_ensemble(request);
    CHECK(result.ok());

    const double elapsed = t1 - t0;
    for (std::size_t node = 0; node < distribution.nodes.size(); ++node)
    {
        const double x0 = distribution.nodes[node][0];
        const double y0 = distribution.nodes[node][1];
        const double expected_x = x0 * std::exp(a * elapsed);
        const double expected_y = y0 * std::exp(c * elapsed) + b * x0 * (std::exp(a * elapsed) - std::exp(c * elapsed)) / (a - c);
        check_close(result.quadrature_result.evolved_distribution.nodes[node][0], expected_x, 1.0e-9);
        check_close(result.quadrature_result.evolved_distribution.nodes[node][1], expected_y, 1.0e-9);
    }
    check_statistically_close(find_row(result.comparison, {1, 1}));
    check_statistically_close(find_row(result.comparison, {1, 2}));
}

void test_asymmetric_tensor_indexing()
{
    const ExplicitDiscreteDistribution quadrature{{{0.2, -0.1, 0.3}, {-0.3, 0.25, -0.2}, {0.1, 0.4, 0.15}}, {1.0, 1.0, 1.0}, {}};
    const qmomkit::monte_carlo::Ensemble ensemble{{0.2, -0.1, 0.3}, {-0.3, 0.25, -0.2}, {0.1, 0.4, 0.15}};
    const std::vector<std::size_t> maximum_orders{2, 3, 1};
    const qmomkit::monte_carlo::MomentComparisonResult result = compare_quadrature_to_ensemble(quadrature, ensemble, 3.0, maximum_orders);
    CHECK(result.ok());
    CHECK(result.rows.size() == 24);

    const std::vector<std::size_t> extents{3, 4, 2};
    for (std::size_t flat = 0; flat < result.rows.size(); ++flat)
    {
        const qmomkit::MultiIndex expected = qmomkit::unflatten_multi_index(flat, extents);
        CHECK(result.rows[flat].multi_index == expected);
        CHECK(result.rows[flat].order == flat);
        check_close(result.rows[flat].reference_value, result.rows[flat].computed_value, 1.0e-14);
    }
}

void test_constant_velocity_comparison()
{
    const ExplicitDiscreteDistribution distribution{{{-0.2}, {0.1}, {0.35}}, {0.8, 1.1, 0.9}, {}};
    const double t0 = 0.1;
    const double t1 = 0.6;
    const double velocity = -0.35;
    const qmomkit::PhaseSpaceVelocity callback = [velocity](double, const std::vector<double>& coordinates)
    {
        CHECK(coordinates.size() == 1);
        return std::vector<double>{velocity};
    };

    EvolutionComparisonRequest request;
    request.initial_distribution = distribution;
    request.velocity = callback;
    request.t0 = t0;
    request.t1 = t1;
    request.adaptive_options = accurate_options();
    request.sampling_options = SamplingOptions{6000, 4101};
    request.maximum_order = 2;

    const EvolutionComparisonResult result = compare_evolved_quadrature_to_ensemble(request);
    CHECK(result.ok());
    CHECK(result.qmomkit_result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk23);
    check_close(result.number_density, 2.8, 1.0e-15);
    check_close(result.comparison.rows[0].reference_value, result.number_density, 1.0e-15);
    check_close(result.comparison.rows[0].computed_value, result.number_density, 1.0e-15);
    CHECK(result.monte_carlo_result.ensemble.size() == request.sampling_options.sample_count);
    check_close(result.monte_carlo_result.final_time, t1, 0.0);

    for (std::size_t node = 0; node < distribution.nodes.size(); ++node)
    {
        const double expected = distribution.nodes[node][0] + velocity * (t1 - t0);
        check_close(result.qmomkit_result.quadrature.abscissas[node], expected, 1.0e-12);
    }
    for (const qmomkit::monte_carlo::MomentComparisonRow& row : result.comparison.rows)
    {
        check_statistically_close(row);
    }
}

void test_affine_evolution_comparison()
{
    const ExplicitDiscreteDistribution distribution{{{-0.3}, {0.2}, {0.55}}, {0.7, 1.0, 0.8}, {}};
    const double t0 = -0.2;
    const double t1 = 0.4;
    const double rate = -0.25;
    const double shift = 0.8;
    const qmomkit::PhaseSpaceVelocity callback = [rate, shift](double, const std::vector<double>& coordinates)
    {
        CHECK(coordinates.size() == 1);
        return std::vector<double>{rate * coordinates[0] + shift};
    };

    EvolutionComparisonRequest request;
    request.initial_distribution = distribution;
    request.velocity = callback;
    request.t0 = t0;
    request.t1 = t1;
    request.adaptive_options = accurate_options();
    request.sampling_options = SamplingOptions{8000, 4201};
    request.maximum_order = 2;

    const EvolutionComparisonResult result = compare_evolved_quadrature_to_ensemble(request);
    CHECK(result.ok());
    CHECK(result.qmomkit_result.diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk23);

    const double elapsed = t1 - t0;
    for (std::size_t node = 0; node < distribution.nodes.size(); ++node)
    {
        const double expected = (distribution.nodes[node][0] + shift / rate) * std::exp(rate * elapsed) - shift / rate;
        check_close(result.qmomkit_result.quadrature.abscissas[node], expected, 1.0e-9);
    }
    for (const qmomkit::monte_carlo::MomentComparisonRow& row : result.comparison.rows)
    {
        check_statistically_close(row);
    }
}

void test_text_output()
{
    const ExplicitDiscreteDistribution distribution{{{-1.0}, {2.0}}, {0.75, 1.25}, {}};
    EvolutionComparisonRequest request;
    request.initial_distribution = distribution;
    request.velocity = [](double, const std::vector<double>&) { return std::vector<double>{0.1}; };
    request.t0 = 0.0;
    request.t1 = 0.2;
    request.adaptive_options = accurate_options();
    request.sampling_options = SamplingOptions{20000, 4301};
    request.maximum_order = 1;

    const EvolutionComparisonResult result = compare_evolved_quadrature_to_ensemble(request);
    CHECK(result.ok());

    std::ostringstream output;
    write_comparison_text(output, result);
    const std::string text = output.str();
    CHECK(text.find("# time=0.2") != std::string::npos);
    CHECK(text.find("# sample_count=20000") != std::string::npos);
    CHECK(text.find("# seed=4301") != std::string::npos);
    CHECK(text.find("# number_density=2") != std::string::npos);
    CHECK(text.find("# rk23_relative_tolerance=1e-10") != std::string::npos);
    CHECK(text.find("# time sample_count seed number_density alpha_0 reference_value computed_value signed_error absolute_error relative_error "
                    "standard_error z_score") != std::string::npos);
    CHECK(text.find(',') == std::string::npos);

    MultiDimensionalEvolutionComparisonRequest multi_request;
    multi_request.initial_distribution = ExplicitDiscreteDistribution{{{0.0, 0.1}, {0.2, -0.1}}, {1.0, 1.0}, {}};
    multi_request.velocity = [](double, const std::vector<double>& coordinates) { return std::vector<double>{0.05, -0.02 * coordinates[1]}; };
    multi_request.t0 = 0.0;
    multi_request.t1 = 0.1;
    multi_request.adaptive_options = accurate_options();
    multi_request.sampling_options = SamplingOptions{3000, 5501};
    multi_request.maximum_orders = {1, 1};

    const MultiDimensionalEvolutionComparisonResult multi_result = compare_evolved_quadrature_to_ensemble(multi_request);
    CHECK(multi_result.ok());

    std::ostringstream multi_output;
    write_comparison_text(multi_output, multi_result);
    const std::string multi_text = multi_output.str();
    CHECK(multi_text.find("# time sample_count seed number_density alpha_0 alpha_1 reference_value computed_value signed_error absolute_error "
                          "relative_error standard_error z_score") != std::string::npos);
}

void test_sample_count_convergence()
{
    SampleCountConvergenceRequest request;
    request.initial_distribution = ExplicitDiscreteDistribution{{{-1.0}, {0.0}, {1.0}}, {1.0, 1.0, 1.0}, {}};
    request.velocity = [](double, const std::vector<double>&) { return std::vector<double>{0.0}; };
    request.t0 = 0.0;
    request.t1 = 0.2;
    request.adaptive_options = accurate_options();
    request.seed = 6101;
    request.sample_counts = {500, 2000, 8000};
    request.maximum_orders = {2};
    request.selected_moments = {{2}};

    const SampleCountConvergenceResult first = run_sample_count_convergence(request);
    const SampleCountConvergenceResult second = run_sample_count_convergence(request);
    CHECK(first.ok());
    CHECK(second.ok());
    CHECK(first.rows.size() == 3);
    CHECK(first.rows.size() == second.rows.size());
    for (std::size_t row = 0; row < first.rows.size(); ++row)
    {
        CHECK(first.rows[row].sample_count == request.sample_counts[row]);
        CHECK(first.rows[row].multi_index == qmomkit::MultiIndex({2}));
        check_close(first.rows[row].reference_value, 2.0, 1.0e-12);
        check_close(first.rows[row].monte_carlo_value, second.rows[row].monte_carlo_value, 0.0);
        check_close(first.rows[row].standard_error, second.rows[row].standard_error, 0.0);
        CHECK(first.rows[row].standard_error > 0.0);
    }
    CHECK(first.rows[1].standard_error < first.rows[0].standard_error);
    CHECK(first.rows[2].standard_error < first.rows[1].standard_error);

    std::vector<qmomkit::monte_carlo::MomentComparisonRow> output_rows;
    for (const SampleCountConvergenceRow& row : first.rows)
    {
        const double signed_error = row.monte_carlo_value - row.reference_value;
        output_rows.push_back(qmomkit::monte_carlo::MomentComparisonRow{row.multi_index, row.flat_index, row.reference_value, row.monte_carlo_value,
                                                                        signed_error, row.absolute_error, 0.0, row.standard_error, row.z_score});
    }
    std::ostringstream output;
    qmomkit::monte_carlo::write_comparison_text(
        output, qmomkit::monte_carlo::MomentComparisonResult{std::move(output_rows), ComparisonStatus::success, ""},
        qmomkit::monte_carlo::ComparisonOutputMetadata{first.final_time, first.rows.front().sample_count, first.seed, first.number_density,
                                                       first.adaptive_options});
    const std::string text = output.str();
    CHECK(text.find("# time sample_count seed number_density alpha_0 reference_value computed_value signed_error absolute_error relative_error "
                    "standard_error z_score") != std::string::npos);
    CHECK(text.find("# rk23_relative_tolerance=1e-10") != std::string::npos);
}

void test_failure_propagation()
{
    EvolutionComparisonRequest request;
    request.initial_distribution = ExplicitDiscreteDistribution{{{0.0}}, {1.0}, {}};
    request.t0 = 0.0;
    request.t1 = 1.0;
    request.sampling_options = SamplingOptions{10, 1};
    request.maximum_order = 1;

    const EvolutionComparisonResult missing_callback = compare_evolved_quadrature_to_ensemble(request);
    CHECK(!missing_callback.ok());
    CHECK(missing_callback.status == ComparisonStatus::invalid_input);

    request.velocity = [](double, const std::vector<double>&) { return std::vector<double>{0.0}; };
    request.initial_distribution = ExplicitDiscreteDistribution{{{0.0}}, {1.0}, 2.0};
    const EvolutionComparisonResult bad_density = compare_evolved_quadrature_to_ensemble(request);
    CHECK(!bad_density.ok());
    CHECK(bad_density.status == ComparisonStatus::invalid_input);

    request.initial_distribution = ExplicitDiscreteDistribution{{{0.0}}, {1.0}, {}};
    request.sampling_options = SamplingOptions{0, 1};
    const EvolutionComparisonResult zero_samples = compare_evolved_quadrature_to_ensemble(request);
    CHECK(!zero_samples.ok());
    CHECK(zero_samples.status == ComparisonStatus::invalid_input);

    request.sampling_options = SamplingOptions{10, 1};
    request.adaptive_options = accurate_options();
    request.adaptive_options.max_steps = 1;
    request.adaptive_options.max_step = 0.01;
    request.velocity = [](double, const std::vector<double>& coordinates) { return std::vector<double>{coordinates[0] + 1.0}; };
    const EvolutionComparisonResult rk_failure = compare_evolved_quadrature_to_ensemble(request);
    CHECK(!rk_failure.ok());
    CHECK(rk_failure.status == ComparisonStatus::qmomkit_failed || rk_failure.status == ComparisonStatus::monte_carlo_failed);

    MultiDimensionalEvolutionComparisonRequest multidimensional;
    multidimensional.initial_distribution = ExplicitDiscreteDistribution{{{0.0, 1.0}, {2.0}}, {1.0, 1.0}, {}};
    multidimensional.velocity = [](double, const std::vector<double>& coordinates) { return std::vector<double>(coordinates.size(), 0.0); };
    multidimensional.t0 = 0.0;
    multidimensional.t1 = 1.0;
    multidimensional.sampling_options = SamplingOptions{10, 1};
    multidimensional.maximum_orders = {1, 1};
    const MultiDimensionalEvolutionComparisonResult bad_dimension = compare_evolved_quadrature_to_ensemble(multidimensional);
    CHECK(!bad_dimension.ok());

    multidimensional.initial_distribution = ExplicitDiscreteDistribution{{{0.0, 1.0}}, {1.0}, {}};
    multidimensional.velocity = [](double, const std::vector<double>&) { return std::vector<double>{0.0}; };
    const MultiDimensionalEvolutionComparisonResult bad_callback = compare_evolved_quadrature_to_ensemble(multidimensional);
    CHECK(!bad_callback.ok());
}

} // namespace

int main()
{
    test_explicit_quadrature_construction();
    test_direct_moment_comparison_utilities();
    test_multidimensional_sampling_and_dimension_validation();
    test_deterministic_sampling_and_number_density();
    test_two_dimensional_constant_velocity_comparison();
    test_two_dimensional_affine_comparison();
    test_two_dimensional_coupled_linear_comparison();
    test_asymmetric_tensor_indexing();
    test_constant_velocity_comparison();
    test_affine_evolution_comparison();
    test_text_output();
    test_sample_count_convergence();
    test_failure_propagation();
    return 0;
}
