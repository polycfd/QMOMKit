#include "qmomkit/wheeler.hpp"

#include "internal/linear_algebra.hpp"
#include "internal/numeric_utils.hpp"
#include "internal/performance.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace qmomkit
{

MomentSet::MomentSet(std::vector<double> moments) : moments_(std::move(moments)) {}

MomentSet::MomentSet(std::initializer_list<double> moments) : moments_(moments) {}

std::size_t MomentSet::size() const noexcept
{
    return moments_.size();
}

bool MomentSet::empty() const noexcept
{
    return moments_.empty();
}

const std::vector<double>& MomentSet::values() const noexcept
{
    return moments_;
}

double MomentSet::operator[](const std::size_t moment_order) const
{
    if (moment_order >= moments_.size())
    {
        throw std::out_of_range("Moment order is outside the stored sequence.");
    }
    return moments_[moment_order];
}

namespace
{

using internal::DenseMatrix;
using internal::EigenDecomposition;

struct RecurrenceCoefficients
{
    std::vector<double> a;
    std::vector<double> b;
    DenseMatrix sigma;
};

[[nodiscard]] Quadrature1D one_node_result(const MomentSet& moments, const QuadratureStatus status, std::string message)
{
    // Non-success payloads are diagnostic only; callers must check status before
    // interpreting the node and weight as an accepted quadrature.
    return Quadrature1D{{moments[1] / moments[0]}, {moments[0]}, status, std::move(message)};
}

[[nodiscard]] Quadrature1D negative_mass_result(const MomentSet& moments)
{
    // Preserve a finite diagnostic payload for failure reporting, not as a
    // physically valid quadrature rule.
    return Quadrature1D{
        {moments[1] / moments[0]}, {std::abs(moments[0])}, QuadratureStatus::not_realizable, "Moment sequence is not realizable: m0 < 0."};
}

[[nodiscard]] RecurrenceCoefficients compute_recurrence(const MomentSet& moments, const std::size_t n)
{
    RecurrenceCoefficients coefficients{std::vector<double>(n, 0.0), std::vector<double>(n, 0.0), DenseMatrix(2 * n + 1, 2 * n + 1)};

    // Wheeler's sigma table converts raw moments into the three-term
    // recurrence coefficients of the orthogonal polynomial sequence.
    for (std::size_t i = 1; i <= 2 * n; ++i)
    {
        coefficients.sigma(1, i) = moments[i - 1];
    }

    coefficients.a[0] = moments[1] / moments[0];
    coefficients.b[0] = 0.0;

    for (std::size_t k = 2; k <= n; ++k)
    {
        for (std::size_t l = k; l <= 2 * n - k + 1; ++l)
        {
            coefficients.sigma(k, l) = coefficients.sigma(k - 1, l + 1) - coefficients.a[k - 2] * coefficients.sigma(k - 1, l) -
                                       coefficients.b[k - 2] * coefficients.sigma(k - 2, l);
        }

        coefficients.a[k - 1] =
            coefficients.sigma(k, k + 1) / coefficients.sigma(k, k) - coefficients.sigma(k - 1, k) / coefficients.sigma(k - 1, k - 1);
        coefficients.b[k - 1] = coefficients.sigma(k, k) / coefficients.sigma(k - 1, k - 1);
    }

    return coefficients;
}

[[nodiscard]] DenseMatrix make_jacobi_matrix(const std::vector<double>& a, const std::vector<double>& b, const std::size_t n)
{
    // The Jacobi matrix encodes the recurrence; its eigenvalues are
    // quadrature nodes and its first eigenvector row determines weights.
    DenseMatrix jacobi(n, n);

    for (std::size_t i = 0; i < n; ++i)
    {
        jacobi(i, i) = a[i];
    }

    for (std::size_t i = 1; i < n; ++i)
    {
        const double off_diagonal = std::sqrt(b[i]);
        jacobi(i - 1, i) = off_diagonal;
        jacobi(i, i - 1) = off_diagonal;
    }

    return jacobi;
}

[[nodiscard]] EigenDecomposition symmetric_jacobi_eigen(DenseMatrix matrix)
{
    const std::size_t n = matrix.rows();
    constexpr std::size_t max_sweeps_per_node = 50;
    constexpr double tolerance = 1.0e-14;
    const std::size_t max_sweeps = max_sweeps_per_node * std::max<std::size_t>(n, 1);
    return internal::symmetric_jacobi_eigen(matrix, {tolerance, max_sweeps, 1});
}

[[nodiscard]] bool recurrence_is_finite(const RecurrenceCoefficients& coefficients)
{
    return internal::all_finite(coefficients.a) && internal::all_finite(coefficients.b);
}

[[nodiscard]] bool has_negative_recurrence_b(const std::vector<double>& b)
{
    return *std::min_element(b.begin(), b.end()) < 0.0;
}

[[nodiscard]] double adaptive_sigma_tolerance(const MomentSet& moments, const WheelerOptions& options)
{
    constexpr double roundoff_multiplier = 100.0;
    double moment_scale = 1.0;
    for (const double moment : moments.values())
    {
        moment_scale = std::max(moment_scale, std::abs(moment));
    }
    return std::max(options.cutoff, roundoff_multiplier * std::numeric_limits<double>::epsilon() * moment_scale);
}

[[nodiscard]] std::vector<double> compute_weights(const double zeroth_moment, const DenseMatrix& eigenvectors)
{
    // Gaussian quadrature weights are m0 times the squared first component of
    // each normalized eigenvector of the Jacobi matrix.
    std::vector<double> weights(eigenvectors.cols(), 0.0);
    for (std::size_t i = 0; i < eigenvectors.cols(); ++i)
    {
        weights[i] = zeroth_moment * eigenvectors(0, i) * eigenvectors(0, i);
    }
    return weights;
}

[[nodiscard]] bool adaptive_accepts(const std::vector<double>& abscissas, const std::vector<double>& weights, const WheelerOptions& options)
{
    // Adaptive Wheeler keeps a reduced quadrature only when node weights remain
    // significant and abscissas remain separated at the requested tolerances.
    const auto minmax_weight = std::minmax_element(weights.begin(), weights.end());
    if (*minmax_weight.first / *minmax_weight.second <= options.rmin)
    {
        return false;
    }

    std::vector<double> nearest_distances(abscissas.size(), 0.0);
    std::vector<double> farthest_distances(abscissas.size(), 0.0);

    for (std::size_t i = abscissas.size() - 1; i > 0; --i)
    {
        double nearest = std::numeric_limits<double>::infinity();
        double farthest = 0.0;
        for (std::size_t j = 0; j < i; ++j)
        {
            const double distance = std::abs(abscissas[i] - abscissas[j]);
            nearest = std::min(nearest, distance);
            farthest = std::max(farthest, distance);
        }
        nearest_distances[i] = nearest;
        farthest_distances[i] = farthest;
    }

    const double min_nearest = *std::min_element(nearest_distances.begin() + 1, nearest_distances.end());
    double max_farthest = *std::max_element(farthest_distances.begin() + 1, farthest_distances.end());
    if (abscissas.size() == 2)
    {
        max_farthest = 1.0;
    }

    return min_nearest / max_farthest > options.eabs;
}

} // namespace

Quadrature1D compute_quadrature(const MomentSet& moments, const std::size_t max_nodes, const WheelerOptions& options)
{
    internal::performance::record_wheeler_call();
    // Wheeler inversion path:
    // raw moments -> recurrence coefficients -> Jacobi matrix -> eigenproblem
    // -> quadrature abscissas and weights.
    if (!std::isfinite(options.rmin) || !std::isfinite(options.eabs) || !std::isfinite(options.cutoff) || options.rmin < 0.0 || options.eabs < 0.0 ||
        options.cutoff <= 0.0)
    {
        return Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "Invalid Wheeler numerical option."};
    }
    if (max_nodes == 0)
    {
        return Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "max_nodes must be positive."};
    }
    if (moments.size() < 2 * max_nodes)
    {
        return Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "Moment sequence has fewer than 2*n entries."};
    }
    if (!internal::all_finite(moments.values()))
    {
        return Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "Moment sequence contains non-finite values."};
    }

    if (moments[0] < 0.0)
    {
        return negative_mass_result(moments);
    }

    if (moments[0] < options.cutoff)
    {
        return Quadrature1D{{0.0}, {1.0}, QuadratureStatus::not_realizable, "Moment sequence has m0 < cutoff."};
    }

    if (max_nodes == 1 || (options.adaptive && moments[0] <= options.rmin))
    {
        if (max_nodes > 1 && options.adaptive)
        {
            internal::performance::record_wheeler_adaptive_reduction();
        }
        return one_node_result(moments, QuadratureStatus::success, "");
    }

    std::size_t n = max_nodes;
    RecurrenceCoefficients coefficients = compute_recurrence(moments, n);

    if (options.adaptive)
    {
        const double sigma_tolerance = adaptive_sigma_tolerance(moments, options);
        for (std::size_t k = max_nodes; k > 1; --k)
        {
            if (coefficients.sigma(k, k) < -sigma_tolerance)
            {
                std::ostringstream message;
                message << std::scientific << std::setprecision(17) << "Adaptive Wheeler recurrence has a negative sigma diagonal: "
                        << "sigma(" << k << ", " << k << ") = " << coefficients.sigma(k, k) << ", tolerance = " << sigma_tolerance << '.';

                return one_node_result(moments, QuadratureStatus::not_realizable, message.str());
            }
            if (std::abs(coefficients.sigma(k, k)) <= sigma_tolerance)
            {
                n = k - 1;
                if (n == 1)
                {
                    internal::performance::record_wheeler_adaptive_reduction();
                    return one_node_result(moments, QuadratureStatus::success, "");
                }
            }
        }

        coefficients = compute_recurrence(moments, n);
        if (!recurrence_is_finite(coefficients))
        {
            return one_node_result(moments, QuadratureStatus::not_realizable, "Reduced Wheeler recurrence produced non-finite values.");
        }
    }
    else if (!recurrence_is_finite(coefficients))
    {
        return one_node_result(moments, QuadratureStatus::not_realizable, "Wheeler recurrence produced non-finite values.");
    }

    if (has_negative_recurrence_b(coefficients.b))
    {
        const double minimum_b = *std::min_element(coefficients.b.begin(), coefficients.b.end());

        std::ostringstream message;
        message << std::scientific << std::setprecision(17) << "Moment sequence has negative recurrence b = " << minimum_b << '.';

        return one_node_result(moments, QuadratureStatus::not_realizable, message.str());
    }

    for (std::size_t n1 = n; n1 > 0; --n1)
    {
        if (n1 == 1)
        {
            return one_node_result(moments, QuadratureStatus::success, "");
        }

        const DenseMatrix jacobi = make_jacobi_matrix(coefficients.a, coefficients.b, n1);
        EigenDecomposition eigen = symmetric_jacobi_eigen(jacobi);
        std::vector<double> weights = compute_weights(moments[0], eigen.eigenvectors);

        if (!internal::all_finite(eigen.eigenvalues) || !internal::all_finite(weights))
        {
            return one_node_result(moments, QuadratureStatus::not_realizable, "Jacobi matrix diagonalization produced non-finite values.");
        }

        if (options.adaptive)
        {
            if (adaptive_accepts(eigen.eigenvalues, weights, options))
            {
                if (n1 < max_nodes)
                {
                    internal::performance::record_wheeler_adaptive_reduction();
                }
                return Quadrature1D{std::move(eigen.eigenvalues), std::move(weights), QuadratureStatus::success, ""};
            }
        }
        else
        {
            return Quadrature1D{std::move(eigen.eigenvalues), std::move(weights), QuadratureStatus::success, ""};
        }
    }

    return one_node_result(moments, QuadratureStatus::no_adaptive_quadrature, "No adaptive quadrature passed thresholds.");
}

} // namespace qmomkit
