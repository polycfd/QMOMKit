#include "qmomkit/wheeler.hpp"

#include "internal/bounded_quadrature.hpp"
#include "internal/checked_size.hpp"
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

struct WheelerAttemptResult
{
    Quadrature1D quadrature;
    bool negative_sigma_diagonal = false;
};

[[nodiscard]] Quadrature1D one_node_result(const MomentSet& moments, const QuadratureStatus status, std::string message)
{
    const double abscissa = moments[1] / moments[0];
    if (!std::isfinite(abscissa))
    {
        return Quadrature1D{{}, {}, QuadratureStatus::not_realizable, "One-node Wheeler quadrature produced non-finite values."};
    }

    // Non-success payloads are diagnostic only; callers must check status before
    // interpreting the node and weight as an accepted quadrature.
    return Quadrature1D{{abscissa}, {moments[0]}, status, std::move(message)};
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
    return internal::symmetric_jacobi_eigen(matrix, {tolerance, max_sweeps});
}

[[nodiscard]] bool recurrence_is_finite(const RecurrenceCoefficients& coefficients)
{
    return internal::all_finite(coefficients.a) && internal::all_finite(coefficients.b);
}

[[nodiscard]] bool has_negative_recurrence_b(const std::vector<double>& b)
{
    return *std::min_element(b.begin(), b.end()) < 0.0;
}

[[nodiscard]] double adaptive_sigma_tolerance(const MomentSet& moments, const WheelerOptions& options, const std::size_t moment_count)
{
    constexpr double roundoff_multiplier = 100.0;
    double moment_scale = 1.0;
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        moment_scale = std::max(moment_scale, std::abs(moments[order]));
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
    const double max_farthest = *std::max_element(farthest_distances.begin() + 1, farthest_distances.end());
    if (!std::isfinite(min_nearest) || !std::isfinite(max_farthest) || max_farthest <= 0.0)
    {
        return false;
    }

    return min_nearest / max_farthest > options.eabs;
}

[[nodiscard]] bool required_moments_are_finite(const MomentSet& moments, const std::size_t moment_count)
{
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        if (!std::isfinite(moments[order]))
        {
            return false;
        }
    }
    return true;
}

[[nodiscard]] WheelerAttemptResult attempt_quadrature_count(const MomentSet& moments, const std::size_t max_nodes, const WheelerOptions& options)
{
    std::size_t n = max_nodes;
    RecurrenceCoefficients coefficients = compute_recurrence(moments, n);

    if (options.adaptive)
    {
        const double sigma_tolerance = adaptive_sigma_tolerance(moments, options, 2 * max_nodes);
        for (std::size_t k = max_nodes; k > 1; --k)
        {
            if (coefficients.sigma(k, k) < -sigma_tolerance)
            {
                std::ostringstream message;
                message << std::scientific << std::setprecision(17) << "Adaptive Wheeler recurrence has a negative sigma diagonal: "
                        << "sigma(" << k << ", " << k << ") = " << coefficients.sigma(k, k) << ", tolerance = " << sigma_tolerance << '.';

                return {one_node_result(moments, QuadratureStatus::not_realizable, message.str()), true};
            }
            if (std::abs(coefficients.sigma(k, k)) <= sigma_tolerance)
            {
                n = k - 1;
                if (n == 1)
                {
                    internal::performance::record_wheeler_adaptive_reduction();
                    return {one_node_result(moments, QuadratureStatus::success, ""), false};
                }
            }
        }

        coefficients = compute_recurrence(moments, n);
        if (!recurrence_is_finite(coefficients))
        {
            return {one_node_result(moments, QuadratureStatus::not_realizable, "Reduced Wheeler recurrence produced non-finite values."), false};
        }
    }
    else if (!recurrence_is_finite(coefficients))
    {
        return {one_node_result(moments, QuadratureStatus::not_realizable, "Wheeler recurrence produced non-finite values."), false};
    }

    if (has_negative_recurrence_b(coefficients.b))
    {
        const double minimum_b = *std::min_element(coefficients.b.begin(), coefficients.b.end());

        std::ostringstream message;
        message << std::scientific << std::setprecision(17) << "Moment sequence has negative recurrence b = " << minimum_b << '.';

        return {one_node_result(moments, QuadratureStatus::not_realizable, message.str()), false};
    }

    for (std::size_t n1 = n; n1 > 1; --n1)
    {
        const DenseMatrix jacobi = make_jacobi_matrix(coefficients.a, coefficients.b, n1);
        EigenDecomposition eigen = symmetric_jacobi_eigen(jacobi);
        std::vector<double> weights = compute_weights(moments[0], eigen.eigenvectors);

        if (!internal::all_finite(eigen.eigenvalues) || !internal::all_finite(weights))
        {
            return {one_node_result(moments, QuadratureStatus::not_realizable, "Jacobi matrix diagonalization produced non-finite values."), false};
        }

        if (options.adaptive)
        {
            if (adaptive_accepts(eigen.eigenvalues, weights, options))
            {
                if (n1 < max_nodes)
                {
                    internal::performance::record_wheeler_adaptive_reduction();
                }
                return {Quadrature1D{std::move(eigen.eigenvalues), std::move(weights), QuadratureStatus::success, ""}, false};
            }
        }
        else
        {
            return {Quadrature1D{std::move(eigen.eigenvalues), std::move(weights), QuadratureStatus::success, ""}, false};
        }
    }

    return {one_node_result(moments, QuadratureStatus::success, ""), false};
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
    const std::optional<std::size_t> required_moment_count = internal::checked_size_multiply(2, max_nodes);
    if (!required_moment_count.has_value())
    {
        return Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "max_nodes is too large to derive the required moment count."};
    }
    if (moments.size() < *required_moment_count)
    {
        return Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "Moment sequence has fewer than 2*n entries."};
    }
    if (!required_moments_are_finite(moments, *required_moment_count))
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

    if (max_nodes == 1)
    {
        return one_node_result(moments, QuadratureStatus::success, "");
    }

    if (options.adaptive)
    {
        for (std::size_t candidate_count = max_nodes; candidate_count > 1; --candidate_count)
        {
            const WheelerAttemptResult attempt = attempt_quadrature_count(moments, candidate_count, options);
            if (attempt.quadrature.ok())
            {
                if (candidate_count < max_nodes && attempt.quadrature.weights.size() == candidate_count)
                {
                    internal::performance::record_wheeler_adaptive_reduction();
                }
                return attempt.quadrature;
            }
            if (!attempt.negative_sigma_diagonal)
            {
                return attempt.quadrature;
            }
        }

        internal::performance::record_wheeler_adaptive_reduction();
        return one_node_result(moments, QuadratureStatus::success, "");
    }

    return attempt_quadrature_count(moments, max_nodes, options).quadrature;
}

namespace internal
{

BoundedQuadratureSelection select_bounded_quadrature(const std::vector<double>& moment_sequence, const BoundedQuadratureSelectionRequest& request,
                                                     const AbscissaBounds interval, const WheelerOptions& wheeler_options)
{
    const auto evaluate = [&](const std::size_t node_count, std::string attempt_type)
    {
        BoundedQuadratureCandidate candidate;
        candidate.node_count = node_count;
        candidate.attempt_type = std::move(attempt_type);
        const std::optional<std::size_t> required_moment_count = internal::checked_size_multiply(2, node_count);
        if (node_count == 0 || !required_moment_count.has_value() || moment_sequence.size() < *required_moment_count)
        {
            candidate.quadrature.status = QuadratureStatus::invalid_input;
            candidate.message = "Bounded-support candidate does not have enough moments.";
            return candidate;
        }
        if (!std::isfinite(interval.lower) || !std::isfinite(interval.upper) || !(interval.lower < interval.upper))
        {
            candidate.quadrature.status = QuadratureStatus::invalid_input;
            candidate.message = "Bounded-support interval must have finite lower < upper values.";
            return candidate;
        }

        candidate.inversion_attempted = true;
        candidate.quadrature = compute_quadrature(
            MomentSet(std::vector<double>(moment_sequence.begin(), moment_sequence.begin() + static_cast<std::ptrdiff_t>(*required_moment_count))),
            node_count, wheeler_options);
        if (!candidate.quadrature.ok())
        {
            candidate.message = "Bounded-support Wheeler inversion failed: " + candidate.quadrature.message;
            return candidate;
        }
        if (candidate.quadrature.weights.size() != node_count || candidate.quadrature.abscissas.size() != node_count)
        {
            candidate.quadrature.status = QuadratureStatus::not_realizable;
            candidate.message = "Bounded-support Wheeler inversion returned an unexpected node count.";
            return candidate;
        }

        constexpr double coordinate_support_tolerance = 1.0e-12;
        const bool nodes_inside_interval = std::all_of(candidate.quadrature.abscissas.begin(), candidate.quadrature.abscissas.end(),
                                                       [interval](const double node)
                                                       {
                                                           return std::isfinite(node) && node >= interval.lower - coordinate_support_tolerance &&
                                                                  node <= interval.upper + coordinate_support_tolerance;
                                                       });
        if (!nodes_inside_interval)
        {
            candidate.quadrature.status = QuadratureStatus::not_realizable;
            candidate.message = "Bounded-support Wheeler inversion returned nodes outside the configured interval.";
            return candidate;
        }

        candidate.accepted = true;
        candidate.message = "Bounded-support candidate accepted.";
        return candidate;
    };

    BoundedQuadratureSelection selection;
    if (!request.adaptive)
    {
        BoundedQuadratureCandidate fixed = evaluate(request.configured_node_count, "fixed");
        selection.quadrature = fixed.quadrature;
        selection.selected_node_count_ceiling = request.configured_node_count;
        selection.success = fixed.accepted;
        selection.diagnostic_candidate = std::move(fixed);
        selection.accepted_fallback_node_count = selection.success ? request.configured_node_count : 0;
        if (!selection.success)
        {
            selection.quadrature.message = selection.diagnostic_candidate.message;
        }
        return selection;
    }

    const std::size_t first_count = request.initial_node_count;
    std::optional<BoundedQuadratureCandidate> rejected_candidate;
    for (std::size_t candidate_count = first_count; candidate_count >= 1; --candidate_count)
    {
        const char* attempt_type =
            candidate_count == first_count ? (request.previous_active_node_count.has_value() ? "current" : "configured") : "fallback";
        BoundedQuadratureCandidate candidate = evaluate(candidate_count, attempt_type);
        if (candidate.accepted)
        {
            const bool reduced = candidate_count < first_count;
            selection.diagnostic_candidate = reduced ? *rejected_candidate : candidate;
            selection.quadrature = std::move(candidate.quadrature);
            selection.selected_node_count_ceiling = candidate_count;
            selection.reduced = reduced;
            selection.success = true;
            selection.accepted_fallback_node_count = candidate_count;
            break;
        }
        if (!rejected_candidate.has_value())
        {
            rejected_candidate = candidate;
        }
        if (candidate_count == 1)
        {
            selection.quadrature = std::move(candidate.quadrature);
            selection.selected_node_count_ceiling = first_count;
            selection.diagnostic_candidate = std::move(*rejected_candidate);
            selection.quadrature.message = candidate.message;
            return selection;
        }
    }

    if (request.allow_growth && request.adaptive && wheeler_options.adaptive && request.previous_active_node_count.has_value() &&
        selection.quadrature.ok() && selection.quadrature.weights.size() == *request.previous_active_node_count &&
        *request.previous_active_node_count < request.configured_node_count)
    {
        const std::size_t growth_count = *request.previous_active_node_count + 1;
        BoundedQuadratureCandidate growth = evaluate(growth_count, "growth");
        selection.diagnostic_candidate = growth;
        if (growth.accepted)
        {
            selection.quadrature = std::move(growth.quadrature);
            selection.selected_node_count_ceiling = growth_count;
            selection.reduced = false;
            selection.accepted_fallback_node_count = growth_count;
        }
        else
        {
            selection.accepted_fallback_node_count = selection.quadrature.weights.size();
        }
    }
    return selection;
}

} // namespace internal

} // namespace qmomkit
