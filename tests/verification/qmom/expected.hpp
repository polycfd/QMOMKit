#pragma once

#include "qmomkit/qmom.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace qmomkit::verification
{

constexpr double qmom_analytical_tolerance = 1.0e-12;

struct QmomAnalyticalCase
{
    std::string name;
    std::vector<double> weights;
    std::vector<double> abscissas;
    enum class Family
    {
        discrete,
        uniform,
        exponential,
        gaussian
    };
    Family family = Family::discrete;
};

struct ErrorMetrics
{
    double maximum_absolute = 0.0;
    double maximum_relative = 0.0;
};

struct QmomVerificationResult
{
    std::string name;
    QmomStatus status = QmomStatus::invalid_input;
    std::string message;
    ErrorMetrics weight_errors;
    ErrorMetrics abscissa_errors;
    ErrorMetrics moment_errors;
    ErrorMetrics mass_errors;
    ErrorMetrics family_errors;
    bool property_checks_passed = true;

    [[nodiscard]] double maximum_absolute_error() const noexcept
    {
        return std::max({weight_errors.maximum_absolute, abscissa_errors.maximum_absolute, moment_errors.maximum_absolute,
                         mass_errors.maximum_absolute, family_errors.maximum_absolute});
    }

    [[nodiscard]] double maximum_relative_error() const noexcept
    {
        return std::max({weight_errors.maximum_relative, abscissa_errors.maximum_relative, moment_errors.maximum_relative,
                         mass_errors.maximum_relative, family_errors.maximum_relative});
    }

    [[nodiscard]] bool passed() const noexcept
    {
        return status == QmomStatus::success && property_checks_passed && maximum_absolute_error() <= qmom_analytical_tolerance &&
               maximum_relative_error() <= qmom_analytical_tolerance;
    }
};

inline const QmomAnalyticalCase* find_qmom_analytical_case(const std::string& identifier)
{
    // These cases are independent references for one-dimensional QMOM:
    // manufactured discrete rules plus classical Gaussian quadrature rules for
    // uniform, exponential, and standard-normal moment families. The configured
    // .qmom files supply raw moments from these references; compute_qmom is
    // expected to recover the same weights, nodes, mass, and reconstructed
    // moments without normalizing m0 away.
    static const std::vector<QmomAnalyticalCase> cases = {
        {"single_dirac", {2.5}, {-0.7}},
        {"two_node_discrete", {1.0 / 3.0, 2.0 / 3.0}, {1.0, 3.0}},
        {"three_node_discrete", {0.2, 0.3, 0.5}, {-1.5, 0.4, 2.2}},
        {"uniform_n2", {0.5, 0.5}, {0.5 * (1.0 - 1.0 / std::sqrt(3.0)), 0.5 * (1.0 + 1.0 / std::sqrt(3.0))}, QmomAnalyticalCase::Family::uniform},
        {"uniform_n3",
         {5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0},
         {0.5 * (1.0 - std::sqrt(3.0 / 5.0)), 0.5, 0.5 * (1.0 + std::sqrt(3.0 / 5.0))},
         QmomAnalyticalCase::Family::uniform},
        {"exponential_n2",
         {0.25 * (2.0 + std::sqrt(2.0)), 0.25 * (2.0 - std::sqrt(2.0))},
         {2.0 - std::sqrt(2.0), 2.0 + std::sqrt(2.0)},
         QmomAnalyticalCase::Family::exponential},
        {"exponential_n3",
         {0.71109300992917301545, 0.27851773356924084880, 0.01038925650158613575},
         {0.41577455678347908331, 2.2942803602790417198, 6.2899450829374791969},
         QmomAnalyticalCase::Family::exponential},
        {"gaussian_n2", {0.5, 0.5}, {-1.0, 1.0}, QmomAnalyticalCase::Family::gaussian},
        {"gaussian_n3", {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}, {-std::sqrt(3.0), 0.0, std::sqrt(3.0)}, QmomAnalyticalCase::Family::gaussian}};
    for (const QmomAnalyticalCase& test_case : cases)
    {
        if (test_case.name == identifier)
        {
            return &test_case;
        }
    }
    return nullptr;
}

inline void update_errors(ErrorMetrics& errors, const double actual, const double expected)
{
    const double absolute = std::abs(actual - expected);
    const double relative = absolute / std::max(std::abs(expected), 1.0);
    errors.maximum_absolute = std::max(errors.maximum_absolute, absolute);
    errors.maximum_relative = std::max(errors.maximum_relative, relative);
}

inline void append_failure(QmomVerificationResult& verification, const std::string& message)
{
    verification.property_checks_passed = false;
    if (!verification.message.empty())
    {
        verification.message += '\n';
    }
    verification.message += message;
}

inline bool near_zero(const double value)
{
    return std::abs(value) <= qmom_analytical_tolerance;
}

inline void verify_distribution_properties(const QmomAnalyticalCase& test_case, QmomVerificationResult& verification,
                                           const std::vector<std::pair<double, double>>& recovered)
{
    // Family-specific checks protect mathematical structure that is easy to
    // lose even when moment reconstruction is close: support for uniform and
    // exponential cases, symmetry for uniform and Gaussian rules, and strictly
    // positive quadrature weights.
    const std::size_t node_count = recovered.size();
    for (const auto& node : recovered)
    {
        if (!(node.second > 0.0))
        {
            append_failure(verification, "Recovered quadrature contains a non-positive weight.");
        }
    }

    switch (test_case.family)
    {
    case QmomAnalyticalCase::Family::discrete:
        break;
    case QmomAnalyticalCase::Family::uniform:
        for (const auto& node : recovered)
        {
            if (node.first < -qmom_analytical_tolerance || node.first > 1.0 + qmom_analytical_tolerance)
            {
                append_failure(verification, "Uniform quadrature node lies outside [0, 1].");
            }
        }
        for (std::size_t left = 0; left < node_count; ++left)
        {
            const std::size_t right = node_count - 1 - left;
            update_errors(verification.family_errors, recovered[left].first + recovered[right].first, 1.0);
            update_errors(verification.family_errors, recovered[left].second, recovered[right].second);
        }
        break;
    case QmomAnalyticalCase::Family::exponential:
        for (const auto& node : recovered)
        {
            if (!(node.first > 0.0))
            {
                append_failure(verification, "Exponential quadrature node is not positive.");
            }
        }
        break;
    case QmomAnalyticalCase::Family::gaussian:
        for (std::size_t left = 0; left < node_count; ++left)
        {
            const std::size_t right = node_count - 1 - left;
            update_errors(verification.family_errors, recovered[left].first + recovered[right].first, 0.0);
            update_errors(verification.family_errors, recovered[left].second, recovered[right].second);
        }
        break;
    }
}

inline QmomVerificationResult verify_qmom_case(const QmomAnalyticalCase& test_case, const QmomCase& qmom_case)
{
    // Verify a complete public QMOM analytical workflow. The reference is the
    // independently listed quadrature rule above; the test compares recovered
    // nodes and weights, then reconstructs every configured raw moment from the
    // recovered quadrature as the user-visible mathematical guarantee.
    QmomVerificationResult verification;
    verification.name = test_case.name;
    if (qmom_case.node_count != test_case.weights.size() || qmom_case.moments.empty())
    {
        verification.message = "Configured node count or moment sequence does not match the analytical case.";
        return verification;
    }

    long double expected_mass = 0.0L;
    for (const double weight : test_case.weights)
    {
        expected_mass += static_cast<long double>(weight);
    }
    update_errors(verification.mass_errors, qmom_case.moments[0], static_cast<double>(expected_mass));

    const QmomResult quadrature = compute_qmom(qmom_case);
    verification.status = quadrature.status;
    verification.message = quadrature.message;
    if (!quadrature.ok() || quadrature.weights.size() != test_case.weights.size() || quadrature.abscissas.size() != test_case.abscissas.size())
    {
        return verification;
    }

    std::vector<std::pair<double, double>> recovered;
    std::vector<std::pair<double, double>> expected;
    long double recovered_mass = 0.0L;
    for (std::size_t node = 0; node < test_case.weights.size(); ++node)
    {
        recovered.emplace_back(quadrature.abscissas[node], quadrature.weights[node]);
        expected.emplace_back(test_case.abscissas[node], test_case.weights[node]);
        recovered_mass += static_cast<long double>(quadrature.weights[node]);
    }
    update_errors(verification.mass_errors, static_cast<double>(recovered_mass), qmom_case.moments[0]);
    std::sort(recovered.begin(), recovered.end());
    std::sort(expected.begin(), expected.end());
    for (std::size_t node = 0; node < expected.size(); ++node)
    {
        update_errors(verification.abscissa_errors, recovered[node].first, expected[node].first);
        update_errors(verification.weight_errors, recovered[node].second, expected[node].second);
    }

    for (std::size_t order = 0; order < qmom_case.moments.size(); ++order)
    {
        long double reconstructed = 0.0L;
        for (const auto& node : recovered)
        {
            long double power = 1.0L;
            for (std::size_t exponent = 0; exponent < order; ++exponent)
            {
                power *= static_cast<long double>(node.first);
            }
            reconstructed += static_cast<long double>(node.second) * power;
        }
        update_errors(verification.moment_errors, static_cast<double>(reconstructed), qmom_case.moments[order]);
        if (test_case.family == QmomAnalyticalCase::Family::gaussian && order % 2 == 1 && !near_zero(static_cast<double>(reconstructed)))
        {
            append_failure(verification, "Gaussian odd reconstructed moment does not vanish within tolerance.");
        }
    }
    verify_distribution_properties(test_case, verification, recovered);
    return verification;
}

} // namespace qmomkit::verification
