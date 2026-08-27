#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/quadrature.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

constexpr long double reference_tolerance = 2.0e-12L;
constexpr double node_weight_tolerance = 2.0e-8;
// Raw odd Hermite moments near order 19 are cancellation-sensitive in double
// precision even when the reconstructed nodes and weights remain accurate.
constexpr double moment_tolerance = 2.0e-5;

enum class Family
{
    legendre,
    laguerre,
    normal
};

struct ReferenceQuadrature
{
    std::vector<long double> nodes;
    std::vector<long double> weights;
};

struct ErrorMetrics
{
    double scaled_node = 0.0;
    double scaled_weight = 0.0;
    double scaled_moment = 0.0;
};

struct RunResult
{
    std::string family;
    std::size_t node_count = 0;
    std::string status;
    std::size_t active_node_count = 0;
    ErrorMetrics errors;
    bool passed = false;
};

const char* family_name(const Family family)
{
    switch (family)
    {
    case Family::legendre:
        return "legendre_unit_weight_minus1_to_1";
    case Family::laguerre:
        return "laguerre_alpha0_exp_minus_x";
    case Family::normal:
        return "hermite_probabilists_normalized_standard_normal";
    }
    return "unknown";
}

const char* status_name(const qmomkit::QmomStatus status)
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

bool finite(const double value)
{
    return std::isfinite(value);
}

bool finite(const std::vector<double>& values)
{
    return std::all_of(values.begin(), values.end(), [](const double value) { return finite(value); });
}

long double integer_power(const long double base, const std::size_t exponent)
{
    long double value = 1.0L;
    for (std::size_t power = 0; power < exponent; ++power)
    {
        value *= base;
    }
    return value;
}

std::vector<long double> analytical_moments(const Family family, const std::size_t moment_count)
{
    std::vector<long double> moments(moment_count, 0.0L);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        if (family == Family::legendre)
        {
            moments[order] = order % 2 == 0 ? 2.0L / static_cast<long double>(order + 1) : 0.0L;
        }
        else if (family == Family::laguerre)
        {
            long double factorial = 1.0L;
            for (std::size_t factor = 2; factor <= order; ++factor)
            {
                factorial *= static_cast<long double>(factor);
            }
            moments[order] = factorial;
        }
        else if (order % 2 == 0)
        {
            long double double_factorial = 1.0L;
            for (std::size_t factor = 1; factor < order; factor += 2)
            {
                double_factorial *= static_cast<long double>(factor);
            }
            moments[order] = double_factorial;
        }
    }
    return moments;
}

ReferenceQuadrature golub_welsch_reference(const Family family, const std::size_t node_count)
{
    // This case-local Golub--Welsch construction is the independent classical
    // reference; QMOMKit is used only after the input moments are formed.
    std::vector<std::vector<long double>> matrix(node_count, std::vector<long double>(node_count, 0.0L));
    long double mass = 1.0L;

    for (std::size_t row = 0; row < node_count; ++row)
    {
        if (family == Family::laguerre)
        {
            matrix[row][row] = 2.0L * static_cast<long double>(row) + 1.0L;
        }
    }
    if (family == Family::legendre)
    {
        mass = 2.0L;
    }

    for (std::size_t index = 1; index < node_count; ++index)
    {
        const long double i = static_cast<long double>(index);
        long double off_diagonal = 0.0L;
        if (family == Family::legendre)
        {
            off_diagonal = i / std::sqrt(4.0L * i * i - 1.0L);
        }
        else if (family == Family::laguerre)
        {
            off_diagonal = i;
        }
        else
        {
            off_diagonal = std::sqrt(i);
        }
        matrix[index - 1][index] = off_diagonal;
        matrix[index][index - 1] = off_diagonal;
    }

    std::vector<std::vector<long double>> eigenvectors(node_count, std::vector<long double>(node_count, 0.0L));
    for (std::size_t index = 0; index < node_count; ++index)
    {
        eigenvectors[index][index] = 1.0L;
    }

    const std::size_t maximum_iterations = 1000 * node_count * node_count;
    for (std::size_t iteration = 0; iteration < maximum_iterations; ++iteration)
    {
        std::size_t pivot_row = 0;
        std::size_t pivot_column = 0;
        long double maximum_off_diagonal = 0.0L;
        for (std::size_t row = 0; row < node_count; ++row)
        {
            for (std::size_t column = row + 1; column < node_count; ++column)
            {
                const long double candidate = std::abs(matrix[row][column]);
                if (candidate > maximum_off_diagonal)
                {
                    maximum_off_diagonal = candidate;
                    pivot_row = row;
                    pivot_column = column;
                }
            }
        }
        if (maximum_off_diagonal <= 1.0e-28L)
        {
            break;
        }
        if (iteration + 1 == maximum_iterations)
        {
            throw std::runtime_error("The case-local Golub--Welsch Jacobi eigensolver did not converge.");
        }

        const long double app = matrix[pivot_row][pivot_row];
        const long double aqq = matrix[pivot_column][pivot_column];
        const long double apq = matrix[pivot_row][pivot_column];
        const long double angle = 0.5L * std::atan2(2.0L * apq, aqq - app);
        const long double cosine = std::cos(angle);
        const long double sine = std::sin(angle);

        for (std::size_t row = 0; row < node_count; ++row)
        {
            if (row == pivot_row || row == pivot_column)
            {
                continue;
            }
            const long double arp = matrix[row][pivot_row];
            const long double arq = matrix[row][pivot_column];
            matrix[row][pivot_row] = cosine * arp - sine * arq;
            matrix[pivot_row][row] = matrix[row][pivot_row];
            matrix[row][pivot_column] = sine * arp + cosine * arq;
            matrix[pivot_column][row] = matrix[row][pivot_column];
        }
        matrix[pivot_row][pivot_row] = cosine * cosine * app - 2.0L * sine * cosine * apq + sine * sine * aqq;
        matrix[pivot_column][pivot_column] = sine * sine * app + 2.0L * sine * cosine * apq + cosine * cosine * aqq;
        matrix[pivot_row][pivot_column] = 0.0L;
        matrix[pivot_column][pivot_row] = 0.0L;

        for (std::size_t row = 0; row < node_count; ++row)
        {
            const long double vrp = eigenvectors[row][pivot_row];
            const long double vrq = eigenvectors[row][pivot_column];
            eigenvectors[row][pivot_row] = cosine * vrp - sine * vrq;
            eigenvectors[row][pivot_column] = sine * vrp + cosine * vrq;
        }
    }

    std::vector<std::pair<long double, long double>> rule;
    rule.reserve(node_count);
    for (std::size_t node = 0; node < node_count; ++node)
    {
        rule.emplace_back(matrix[node][node], mass * eigenvectors[0][node] * eigenvectors[0][node]);
    }
    std::sort(rule.begin(), rule.end(), [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

    ReferenceQuadrature reference;
    reference.nodes.reserve(node_count);
    reference.weights.reserve(node_count);
    for (const auto& [node, weight] : rule)
    {
        reference.nodes.push_back(node);
        reference.weights.push_back(weight);
    }
    return reference;
}

std::vector<long double> moments_from_reference(const ReferenceQuadrature& reference, const std::size_t moment_count)
{
    std::vector<long double> moments(moment_count, 0.0L);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        for (std::size_t node = 0; node < reference.nodes.size(); ++node)
        {
            moments[order] += reference.weights[node] * integer_power(reference.nodes[node], order);
        }
    }
    return moments;
}

void verify_reference(const Family family, const ReferenceQuadrature& reference, const std::size_t moment_count)
{
    // Exact family moments validate the independently constructed rule before
    // it supplies the raw moments to the reconstruction under study.
    const std::vector<long double> actual = moments_from_reference(reference, moment_count);
    const std::vector<long double> expected = analytical_moments(family, moment_count);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        long double absolute_sum = 0.0L;
        for (std::size_t node = 0; node < reference.nodes.size(); ++node)
        {
            absolute_sum += std::abs(reference.weights[node] * integer_power(reference.nodes[node], order));
        }
        const long double scale = std::max({1.0L, std::abs(expected[order]), absolute_sum});
        const long double scaled_error = std::abs(actual[order] - expected[order]) / scale;
        if (!std::isfinite(static_cast<double>(actual[order])) || scaled_error > reference_tolerance)
        {
            std::ostringstream message;
            message << "The case-local classical reference quadrature failed its analytical moment check for " << family_name(family) << " at order "
                    << order << ", scaled error " << std::scientific << std::setprecision(17) << static_cast<double>(scaled_error) << '.';
            throw std::runtime_error(message.str());
        }
    }
}

void update_scaled_error(double& maximum_scaled, const double actual, const double expected)
{
    const double absolute = std::abs(actual - expected);
    maximum_scaled = std::max(maximum_scaled, absolute / std::max(1.0, std::abs(expected)));
}

RunResult run_case(const Family family, const std::size_t node_count)
{
    const std::size_t moment_count = 2 * node_count;
    const ReferenceQuadrature reference = golub_welsch_reference(family, node_count);
    verify_reference(family, reference, moment_count);
    const std::vector<long double> reference_moments = moments_from_reference(reference, moment_count);

    std::vector<double> input_moments;
    input_moments.reserve(moment_count);
    for (const long double moment : reference_moments)
    {
        input_moments.push_back(static_cast<double>(moment));
    }

    // Native QMOMKit inversion and raw-moment reconstruction begin here.
    const qmomkit::QmomResult reconstructed = qmomkit::compute_qmom(qmomkit::MomentSet{input_moments}, node_count, {});
    RunResult result;
    result.family = family_name(family);
    result.node_count = node_count;
    result.status = status_name(reconstructed.status);
    result.active_node_count = reconstructed.weights.size();
    if (!reconstructed.ok() || reconstructed.weights.size() != node_count || reconstructed.abscissas.size() != node_count ||
        !finite(reconstructed.abscissas) || !finite(reconstructed.weights))
    {
        return result;
    }

    std::vector<std::pair<double, double>> recovered_rule;
    recovered_rule.reserve(node_count);
    for (std::size_t node = 0; node < node_count; ++node)
    {
        recovered_rule.emplace_back(reconstructed.abscissas[node], reconstructed.weights[node]);
    }
    std::sort(recovered_rule.begin(), recovered_rule.end(), [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

    for (std::size_t node = 0; node < node_count; ++node)
    {
        update_scaled_error(result.errors.scaled_node, recovered_rule[node].first, static_cast<double>(reference.nodes[node]));
        update_scaled_error(result.errors.scaled_weight, recovered_rule[node].second, static_cast<double>(reference.weights[node]));
    }

    const qmomkit::Quadrature1D quadrature{reconstructed.abscissas, reconstructed.weights, qmomkit::QuadratureStatus::success, ""};
    const std::vector<double> reconstructed_moments = qmomkit::reconstruct_raw_moments(quadrature, moment_count);
    if (!finite(reconstructed_moments))
    {
        return result;
    }
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        update_scaled_error(result.errors.scaled_moment, reconstructed_moments[order], static_cast<double>(reference_moments[order]));
    }

    result.passed = result.errors.scaled_node <= node_weight_tolerance && result.errors.scaled_weight <= node_weight_tolerance &&
                    result.errors.scaled_moment <= moment_tolerance;
    return result;
}

void write_csv(const std::vector<RunResult>& results)
{
    std::ofstream output("quadrature_1d.csv");
    if (!output)
    {
        throw std::runtime_error("Could not open quadrature_1d.csv for writing.");
    }
    output << "family,node_count,status,active_node_count,max_scaled_node_error,max_scaled_weight_error,max_scaled_moment_error\n";
    output << std::scientific << std::setprecision(17);
    for (const RunResult& result : results)
    {
        output << result.family << ',' << result.node_count << ',' << result.status << ',' << result.active_node_count << ','
               << result.errors.scaled_node << ',' << result.errors.scaled_weight << ',' << result.errors.scaled_moment << '\n';
    }
}

void print_table(const std::vector<RunResult>& results)
{
    std::cout << "Classical quadrature reconstruction\n";
    std::cout << "Hermite convention: normalized standard normal, not physicists' exp(-x^2).\n";
    std::cout << std::left << std::setw(46) << "Family" << std::right << std::setw(4) << "N" << std::setw(18) << "node scaled" << std::setw(18)
              << "weight scaled" << std::setw(18) << "moment scaled" << "  status\n";
    std::cout << std::scientific << std::setprecision(6);
    for (const RunResult& result : results)
    {
        std::cout << std::left << std::setw(46) << result.family << std::right << std::setw(4) << result.node_count << std::setw(18)
                  << result.errors.scaled_node << std::setw(18) << result.errors.scaled_weight << std::setw(18) << result.errors.scaled_moment << "  "
                  << result.status << '\n';
    }
}

} // namespace

int main()
{
    try
    {
        std::vector<RunResult> results;
        for (const Family family : {Family::legendre, Family::laguerre, Family::normal})
        {
            for (std::size_t node_count = 1; node_count <= 10; ++node_count)
            {
                results.push_back(run_case(family, node_count));
            }
        }
        write_csv(results);
        print_table(results);

        const bool passed = std::all_of(results.begin(), results.end(), [](const RunResult& result) { return result.passed; });
        if (!passed)
        {
            std::cerr << "One or more classical quadrature reconstructions did not meet the verification criteria.\n";
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
