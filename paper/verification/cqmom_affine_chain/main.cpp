#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{

// The reconstruction is exact, so the errors are round-off amplified by the
// conditioning of the hierarchy, which grows with the dimension.
constexpr double verification_tolerance = 1.0e-10;
constexpr std::size_t nodes_per_rule = 3;
constexpr long double affine_intercept = 1.0L / 16.0L;
constexpr long double affine_parent_coefficient = 3.0L / 16.0L;
constexpr long double affine_offset_coefficient = 13.0L / 16.0L;
constexpr long double offsets[nodes_per_rule]{-5.0L / 8.0L, 1.0L / 8.0L, 3.0L / 4.0L};
constexpr long double probabilities[nodes_per_rule]{1.0L / 4.0L, 1.0L / 2.0L, 1.0L / 4.0L};

struct ExactTerminal
{
    std::vector<long double> coordinates;
    long double weight = 0.0L;
    qmomkit::ConditioningPath path;
};

struct RunResult
{
    std::size_t dimension = 0;
    std::size_t leaf_count = 0;
    std::size_t required_moment_count = 0;
    double maximum_scaled_node_error = 0.0;
    double maximum_scaled_weight_error = 0.0;
    double maximum_scaled_moment_error = 0.0;
};

long double integer_power(const long double base, const std::size_t exponent)
{
    long double value = 1.0L;
    for (std::size_t power = 0; power < exponent; ++power)
    {
        value *= base;
    }
    return value;
}

// The affine-chain reference is independent of CQMOM.  Conditioning the last
// coordinate on its parent gives a polynomial recurrence, memoized here only
// to make the fixed 10D reference practical.
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

    long double moment(const qmomkit::MultiIndex& exponents)
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
            qmomkit::MultiIndex reduced = exponents;
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
    std::size_t maximum_exponent_ = 0;
    std::vector<std::vector<long double>> conditional_coefficients_;
    std::map<qmomkit::MultiIndex, long double> cache_;
};

// Each error is scaled by the size of the exact rule it belongs to, as in the
// one-dimensional quadrature study, so the metric does not change when a
// coordinate or the number density is rescaled. A non-finite reconstructed
// value counts as an infinite error.
void update_scaled_error(double& maximum_scaled, const double actual, const double expected, const long double scale)
{
    const double scaled = std::isfinite(actual) ? std::abs(actual - expected) / static_cast<double>(scale) : std::numeric_limits<double>::infinity();
    maximum_scaled = std::max(maximum_scaled, scaled);
}

std::vector<ExactTerminal> exact_terminals(const std::size_t dimension)
{
    std::vector<ExactTerminal> terminals{{std::vector<long double>(dimension, 0.0L), 1.0L, qmomkit::ConditioningPath{}}};
    for (std::size_t level = 0; level < dimension; ++level)
    {
        std::vector<ExactTerminal> children;
        children.reserve(terminals.size() * nodes_per_rule);
        for (const ExactTerminal& parent : terminals)
        {
            for (std::size_t node = 0; node < nodes_per_rule; ++node)
            {
                ExactTerminal child = parent;
                child.path = parent.path.with_child(node);
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

// The sum over the exact terminal nodes of w |x_1|^k_1 ... |x_d|^k_d, the
// magnitude the moment of order k is formed from. The terminal nodes are ordered
// by conditioning path, so the subtree below a node is a contiguous block of
// them, and the sum is accumulated along the tree.
long double absolute_moment_sum(const std::vector<ExactTerminal>& terminals, const qmomkit::MultiIndex& order, const std::size_t level,
                                const std::size_t first, const std::size_t count)
{
    if (level == order.size())
    {
        return 1.0L;
    }
    const std::size_t block = count / nodes_per_rule;
    long double sum = 0.0L;
    for (std::size_t node = 0; node < nodes_per_rule; ++node)
    {
        const std::size_t block_first = first + node * block;
        sum += probabilities[node] * integer_power(std::abs(terminals[block_first].coordinates[level]), order[level]) *
               absolute_moment_sum(terminals, order, level + 1, block_first, block);
    }
    return sum;
}

std::vector<double> exact_required_moments(const qmomkit::CqmomMomentLayout& layout)
{
    AffineMomentEvaluator evaluator(layout.hierarchy.number_of_dimensions());
    std::vector<double> values;
    values.reserve(layout.required_indices.size());
    for (const qmomkit::MultiIndexView index : layout.required_indices)
    {
        values.push_back(static_cast<double>(evaluator.moment(index.to_multi_index())));
    }
    return values;
}

// Compares every terminal node with the exact weighted point of the same
// conditioning path, visiting the nodes in the public node enumeration, and
// records the largest coordinate and weight errors. Coordinates are scaled by
// the extent max_i |x_i,d| of their own dimension in the exact rule, weights by
// its total weight m_0.
void compare_terminal_quadrature(const qmomkit::CqmomQuadrature& quadrature, const std::vector<ExactTerminal>& expected, RunResult& result)
{
    if (!quadrature.has_rectangular_storage() || quadrature.total_nodes() != expected.size())
    {
        throw std::runtime_error("The affine-chain reconstruction returned an unexpected terminal quadrature layout.");
    }

    std::vector<long double> extents(expected.front().coordinates.size(), 0.0L);
    long double total_weight = 0.0L;
    for (const ExactTerminal& terminal : expected)
    {
        for (std::size_t dimension = 0; dimension < extents.size(); ++dimension)
        {
            extents[dimension] = std::max(extents[dimension], std::abs(terminal.coordinates[dimension]));
        }
        total_weight += terminal.weight;
    }

    const std::vector<qmomkit::MultiIndex> node_indices = qmomkit::enumerate_quadrature_node_multi_indices(quadrature.hierarchy());
    for (std::size_t terminal = 0; terminal < node_indices.size(); ++terminal)
    {
        const qmomkit::MultiIndex& node_index = node_indices[terminal];
        if (qmomkit::ConditioningPath(node_index) != expected[terminal].path)
        {
            throw std::runtime_error("The affine-chain terminal path does not match the exact reference.");
        }
        update_scaled_error(result.maximum_scaled_weight_error, quadrature.weight(node_index), static_cast<double>(expected[terminal].weight),
                            total_weight);
        for (std::size_t dimension = 0; dimension < extents.size(); ++dimension)
        {
            update_scaled_error(result.maximum_scaled_node_error, quadrature.abscissa(dimension, node_index),
                                static_cast<double>(expected[terminal].coordinates[dimension]), extents[dimension]);
        }
    }
}

RunResult run_dimension(const std::size_t dimension)
{
    const qmomkit::PhaseSpaceHierarchy hierarchy(std::vector<std::size_t>(dimension, nodes_per_rule));
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::vector<double> exact_moments = exact_required_moments(layout);
    qmomkit::CqmomCase cqmom_case{qmomkit::CqmomMomentData{layout, exact_moments}};
    cqmom_case.qmom_options.adaptive = false;
    const qmomkit::CqmomResult reconstructed = qmomkit::compute_cqmom(cqmom_case);
    if (!reconstructed.ok())
    {
        throw std::runtime_error("Affine-chain CQMOM reconstruction failed: " + reconstructed.message);
    }

    const std::vector<ExactTerminal> expected = exact_terminals(dimension);
    RunResult result{dimension, expected.size(), layout.required_indices.size()};
    compare_terminal_quadrature(reconstructed.quadrature, expected, result);

    const qmomkit::CqmomMomentData recovered_moments = qmomkit::reconstruct_cqmom_moments(reconstructed.quadrature, layout);
    // Moments are scaled by sum_i w_i |x_i,1|^k_1 ... |x_i,d|^k_d of the exact rule.
    for (std::size_t moment = 0; moment < exact_moments.size(); ++moment)
    {
        update_scaled_error(result.maximum_scaled_moment_error, recovered_moments.values[moment], exact_moments[moment],
                            absolute_moment_sum(expected, layout.required_indices[moment].to_multi_index(), 0, 0, expected.size()));
    }
    if (result.maximum_scaled_node_error > verification_tolerance || result.maximum_scaled_weight_error > verification_tolerance ||
        result.maximum_scaled_moment_error > verification_tolerance)
    {
        throw std::runtime_error("The affine-chain terminal quadrature or its moments exceed the verification tolerance.");
    }
    return result;
}

void write_csv(const std::vector<RunResult>& results)
{
    std::ofstream output("cqmom_affine_chain.csv");
    if (!output)
    {
        throw std::runtime_error("Could not open cqmom_affine_chain.csv for writing.");
    }
    output << "dimension,leaves,required_moments,max_scaled_node_error,max_scaled_weight_error,max_scaled_moment_error\n";
    output << std::scientific << std::setprecision(17);
    for (const RunResult& result : results)
    {
        output << result.dimension << ',' << result.leaf_count << ',' << result.required_moment_count << ',' << result.maximum_scaled_node_error
               << ',' << result.maximum_scaled_weight_error << ',' << result.maximum_scaled_moment_error << '\n';
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
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
