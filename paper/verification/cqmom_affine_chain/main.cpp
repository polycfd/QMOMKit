#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{

constexpr double verification_tolerance = 5.0e-8;
constexpr std::size_t nodes_per_rule = 3;
constexpr long double affine_intercept = 1.0L / 16.0L;
constexpr long double affine_parent_coefficient = 3.0L / 16.0L;
constexpr long double affine_offset_coefficient = 13.0L / 16.0L;
constexpr long double offsets[nodes_per_rule]{-5.0L / 8.0L, 1.0L / 8.0L, 3.0L / 4.0L};
constexpr long double probabilities[nodes_per_rule]{1.0L / 4.0L, 1.0L / 2.0L, 1.0L / 4.0L};

struct VectorHash
{
    std::size_t operator()(const std::vector<std::size_t>& values) const noexcept
    {
        std::size_t hash = values.size();
        for (const std::size_t value : values)
        {
            hash ^= value + 0x9e3779b9U + (hash << 6U) + (hash >> 2U);
        }
        return hash;
    }
};

struct ExactTerminal
{
    std::vector<long double> coordinates;
    long double weight = 0.0L;
    std::vector<std::size_t> path;
};

struct RecoveredTerminal
{
    qmomkit::ConditioningPath path;
    std::vector<double> coordinates;
    double weight = 0.0;
};

struct ErrorMetrics
{
    double absolute_node = 0.0;
    double scaled_node = 0.0;
    double absolute_weight = 0.0;
    double scaled_weight = 0.0;
    double absolute_moment = 0.0;
    double scaled_moment = 0.0;
};

struct RunResult
{
    std::size_t dimension = 0;
    std::size_t leaf_count = 0;
    std::size_t required_moment_count = 0;
    std::string status = "not_run";
    ErrorMetrics errors;
    bool passed = false;
};

class AffineMomentEvaluator
{
public:
    explicit AffineMomentEvaluator(const std::size_t dimension) : maximum_exponent_(5 + 2 * (dimension - 1))
    {
        conditional_coefficients_.resize(maximum_exponent_ + 1);
        for (std::size_t exponent = 0; exponent <= maximum_exponent_; ++exponent)
        {
            conditional_coefficients_[exponent].assign(exponent + 1, 0.0L);
            for (std::size_t power = 0; power <= exponent; ++power)
            {
                long double binomial = 1.0L;
                for (std::size_t factor = 1; factor <= power; ++factor)
                {
                    binomial *= static_cast<long double>(exponent + 1 - factor);
                    binomial /= static_cast<long double>(factor);
                }
                long double averaged_offset_power = 0.0L;
                for (std::size_t node = 0; node < nodes_per_rule; ++node)
                {
                    const long double local_offset = affine_intercept + affine_offset_coefficient * offsets[node];
                    averaged_offset_power += probabilities[node] * integer_power(local_offset, exponent - power);
                }
                conditional_coefficients_[exponent][power] = binomial * integer_power(affine_parent_coefficient, power) * averaged_offset_power;
            }
        }
    }

    long double moment(const std::vector<std::size_t>& exponents)
    {
        if (exponents.empty())
        {
            return 1.0L;
        }
        const auto found = cache_.find(exponents);
        if (found != cache_.end())
        {
            return found->second;
        }

        const std::size_t exponent = exponents.back();
        if (exponent > maximum_exponent_)
        {
            throw std::runtime_error("Affine-chain moment evaluator exceeded its precomputed exponent range.");
        }

        long double value = 0.0L;
        if (exponents.size() == 1)
        {
            for (std::size_t node = 0; node < nodes_per_rule; ++node)
            {
                value += probabilities[node] * integer_power(offsets[node], exponent);
            }
        }
        else
        {
            std::vector<std::size_t> reduced = exponents;
            reduced.pop_back();
            const std::size_t parent_exponent = reduced.back();
            for (std::size_t power = 0; power <= exponent; ++power)
            {
                reduced.back() = parent_exponent + power;
                value += conditional_coefficients_[exponent][power] * moment(reduced);
            }
        }
        cache_.emplace(exponents, value);
        return value;
    }

private:
    static long double integer_power(const long double base, const std::size_t exponent)
    {
        long double value = 1.0L;
        for (std::size_t power = 0; power < exponent; ++power)
        {
            value *= base;
        }
        return value;
    }

    std::size_t maximum_exponent_ = 0;
    std::vector<std::vector<long double>> conditional_coefficients_;
    std::unordered_map<std::vector<std::size_t>, long double, VectorHash> cache_;
};

const char* status_name(const qmomkit::CqmomStatus status)
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

bool all_finite(const std::vector<double>& values)
{
    return std::all_of(values.begin(), values.end(), [](const double value) { return std::isfinite(value); });
}

void update_error(double& maximum_absolute, double& maximum_scaled, const double actual, const double expected)
{
    const double absolute = std::abs(actual - expected);
    maximum_absolute = std::max(maximum_absolute, absolute);
    maximum_scaled = std::max(maximum_scaled, absolute / std::max(1.0, std::abs(expected)));
}

std::vector<ExactTerminal> exact_terminals(const std::size_t dimension)
{
    std::vector<ExactTerminal> terminals{{std::vector<long double>(dimension, 0.0L), 1.0L, {}}};
    for (std::size_t level = 0; level < dimension; ++level)
    {
        std::vector<ExactTerminal> children;
        children.reserve(terminals.size() * nodes_per_rule);
        for (const ExactTerminal& parent : terminals)
        {
            for (std::size_t node = 0; node < nodes_per_rule; ++node)
            {
                ExactTerminal child = parent;
                child.path.push_back(node);
                child.weight *= probabilities[node];
                child.coordinates[level] = level == 0 ? offsets[node]
                                                      : affine_intercept + affine_parent_coefficient * parent.coordinates[level - 1] +
                                                            affine_offset_coefficient * offsets[node];
                children.push_back(std::move(child));
            }
        }
        terminals = std::move(children);
    }
    return terminals;
}

std::vector<double> exact_required_moments(const qmomkit::CqmomMomentLayout& layout)
{
    AffineMomentEvaluator evaluator(layout.hierarchy.number_of_dimensions());
    std::vector<double> values;
    values.reserve(layout.required_indices.size());
    for (const qmomkit::MultiIndex& index : layout.required_indices)
    {
        const long double value = evaluator.moment(index);
        const double converted = static_cast<double>(value);
        if (!std::isfinite(converted))
        {
            throw std::runtime_error("The affine-chain exact moment is not finite.");
        }
        values.push_back(converted);
    }
    return values;
}

std::vector<RecoveredTerminal> recovered_terminal_coordinates(const qmomkit::CqmomQuadrature& quadrature)
{
    if (!quadrature.has_rectangular_storage())
    {
        throw std::runtime_error("The fixed affine-chain reconstruction returned adaptive terminal storage.");
    }

    std::vector<RecoveredTerminal> terminals;
    terminals.reserve(quadrature.total_nodes());
    for (const qmomkit::MultiIndex& node_index : qmomkit::enumerate_quadrature_node_multi_indices(quadrature.hierarchy()))
    {
        std::vector<double> coordinates(quadrature.number_of_dimensions(), 0.0);
        for (std::size_t dimension = 0; dimension < coordinates.size(); ++dimension)
        {
            coordinates[dimension] = quadrature.abscissa(dimension, node_index);
        }
        terminals.push_back({qmomkit::ConditioningPath(node_index), std::move(coordinates), quadrature.weight(node_index)});
    }
    return terminals;
}

RunResult run_dimension(const std::size_t dimension)
{
    const qmomkit::PhaseSpaceHierarchy hierarchy(std::vector<std::size_t>(dimension, nodes_per_rule));
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::vector<ExactTerminal> expected = exact_terminals(dimension);
    RunResult result;
    result.dimension = dimension;
    result.leaf_count = expected.size();
    result.required_moment_count = layout.required_indices.size();
    const std::vector<double> input_values = exact_required_moments(layout);

    qmomkit::CqmomLeanCase lean_case{
        qmomkit::CqmomMomentData{layout, input_values}, {}, qmomkit::QmomOptions{}, false, qmomkit::CqmomDiagnosticLevel::basic,
        qmomkit::VandermondeSolver::automatic,          {}};
    const qmomkit::CqmomResult reconstructed = qmomkit::compute_cqmom(lean_case);
    result.status = status_name(reconstructed.status);
    if (!reconstructed.ok() || reconstructed.quadrature.total_nodes() != expected.size() || !all_finite(reconstructed.quadrature.weights()))
    {
        return result;
    }
    for (const std::vector<double>& coordinate : reconstructed.quadrature.abscissas())
    {
        if (!all_finite(coordinate))
        {
            return result;
        }
    }
    if (reconstructed.diagnostics.size() != (expected.size() - 1) / (nodes_per_rule - 1))
    {
        return result;
    }
    for (const qmomkit::CqmomDiagnostics& diagnostic : reconstructed.diagnostics)
    {
        if (diagnostic.requested_node_count != nodes_per_rule || diagnostic.actual_node_count != nodes_per_rule ||
            diagnostic.adaptive_node_count_reduced || diagnostic.qmom_status != qmomkit::QmomStatus::success)
        {
            return result;
        }
    }

    std::vector<RecoveredTerminal> recovered = recovered_terminal_coordinates(reconstructed.quadrature);
    std::vector<std::pair<qmomkit::ConditioningPath, ExactTerminal>> expected_by_path;
    expected_by_path.reserve(expected.size());
    for (const ExactTerminal& terminal : expected)
    {
        expected_by_path.emplace_back(qmomkit::ConditioningPath(terminal.path), terminal);
    }
    std::sort(recovered.begin(), recovered.end(), [](const auto& lhs, const auto& rhs) { return lhs.path < rhs.path; });
    std::sort(expected_by_path.begin(), expected_by_path.end(), [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

    for (std::size_t terminal = 0; terminal < expected_by_path.size(); ++terminal)
    {
        if (recovered[terminal].path != expected_by_path[terminal].first)
        {
            return result;
        }
        for (std::size_t axis = 0; axis < dimension; ++axis)
        {
            update_error(result.errors.absolute_node, result.errors.scaled_node, recovered[terminal].coordinates[axis],
                         static_cast<double>(expected_by_path[terminal].second.coordinates[axis]));
        }
    }

    for (std::size_t terminal = 0; terminal < expected_by_path.size(); ++terminal)
    {
        update_error(result.errors.absolute_weight, result.errors.scaled_weight, recovered[terminal].weight,
                     static_cast<double>(expected_by_path[terminal].second.weight));
    }

    const qmomkit::CqmomMomentData recovered_moments = qmomkit::reconstruct_cqmom_moments(reconstructed.quadrature, layout);
    if (recovered_moments.values.size() != input_values.size() || !all_finite(recovered_moments.values))
    {
        return result;
    }
    for (std::size_t moment = 0; moment < input_values.size(); ++moment)
    {
        update_error(result.errors.absolute_moment, result.errors.scaled_moment, recovered_moments.values[moment], input_values[moment]);
    }

    result.passed = result.errors.scaled_node <= verification_tolerance && result.errors.scaled_weight <= verification_tolerance &&
                    result.errors.scaled_moment <= verification_tolerance;
    return result;
}

void write_csv(const std::vector<RunResult>& results)
{
    std::ofstream output("cqmom_affine_chain.csv");
    if (!output)
    {
        throw std::runtime_error("Could not open cqmom_affine_chain.csv for writing.");
    }
    output << "dimension,nodes_per_conditional_rule,leaf_count,required_moment_count,status,max_absolute_node_error,max_scaled_node_error,"
              "max_absolute_weight_error,max_scaled_weight_error,max_absolute_moment_error,max_scaled_moment_error\n";
    output << std::scientific << std::setprecision(17);
    for (const RunResult& result : results)
    {
        output << result.dimension << ',' << nodes_per_rule << ',' << result.leaf_count << ',' << result.required_moment_count << ',' << result.status
               << ',' << result.errors.absolute_node << ',' << result.errors.scaled_node << ',' << result.errors.absolute_weight << ','
               << result.errors.scaled_weight << ',' << result.errors.absolute_moment << ',' << result.errors.scaled_moment << '\n';
    }
}

void print_table(const std::vector<RunResult>& results)
{
    std::cout << "Recursive CQMOM affine-chain reconstruction\n";
    std::cout << std::right << std::setw(10) << "dimension" << std::setw(12) << "nodes/rule" << std::setw(12) << "leaves" << std::setw(20)
              << "required moments" << std::setw(22) << "scaled moment error" << "  status\n";
    std::cout << std::scientific << std::setprecision(6);
    for (const RunResult& result : results)
    {
        std::cout << std::right << std::setw(10) << result.dimension << std::setw(12) << nodes_per_rule << std::setw(12) << result.leaf_count
                  << std::setw(20) << result.required_moment_count << std::setw(22) << result.errors.scaled_moment << "  " << result.status << '\n';
    }
}

} // namespace

int main()
{
    try
    {
        std::vector<RunResult> results;
        for (const std::size_t dimension : {std::size_t{2}, std::size_t{3}, std::size_t{5}, std::size_t{7}, std::size_t{10}})
        {
            results.push_back(run_dimension(dimension));
        }
        write_csv(results);
        print_table(results);

        const bool passed = std::all_of(results.begin(), results.end(), [](const RunResult& result) { return result.passed; });
        if (!passed)
        {
            std::cerr << "One or more affine-chain reconstructions did not meet the verification criteria.\n";
            return 1;
        }
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
