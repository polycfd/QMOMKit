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
#include <string>
#include <utility>
#include <vector>

namespace qmomkit
{

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

[[nodiscard]] bool recurrence_is_finite(const RecurrenceCoefficients& coefficients)
{
    return internal::all_finite(coefficients.a) && internal::all_finite(coefficients.b);
}

[[nodiscard]] bool has_negative_recurrence_b(const std::vector<double>& b)
{
    return *std::min_element(b.begin(), b.end()) < 0.0;
}

// The tolerance below which adaptive Wheeler treats sigma(k, k) as zero, for
// each k: a running bound on its roundoff, the sigma recurrence summed in
// absolute values. sigma(k, k) and the bound both scale as m0 times the
// coordinate to the power 2(k - 1), so the test depends neither on the unit of
// the coordinate nor on m0, and the bound grows with the cancellation of a
// distribution far from the origin. Conditional moments are not data: the
// solve that produced them leaves rounding relative to the scale of their
// coordinate, not to their own magnitudes, so a positive coordinate_scale s
// adds the floor m0 s^(2(k - 1)).
[[nodiscard]] std::vector<double> adaptive_sigma_tolerances(const MomentSet& moments, const RecurrenceCoefficients& coefficients, const std::size_t n,
                                                            const double coordinate_scale)
{
    constexpr double roundoff_multiplier = 100.0;
    DenseMatrix magnitude(2 * n + 1, 2 * n + 1);
    for (std::size_t i = 1; i <= 2 * n; ++i)
    {
        magnitude(1, i) = std::abs(moments[i - 1]);
    }
    for (std::size_t k = 2; k <= n; ++k)
    {
        for (std::size_t l = k; l <= 2 * n - k + 1; ++l)
        {
            magnitude(k, l) = magnitude(k - 1, l + 1) + std::abs(coefficients.a[k - 2]) * magnitude(k - 1, l) +
                              std::abs(coefficients.b[k - 2]) * magnitude(k - 2, l);
        }
    }

    std::vector<double> tolerances(n + 1, 0.0);
    double noise_floor = std::abs(moments[0]);
    for (std::size_t k = 1; k <= n; ++k)
    {
        tolerances[k] = roundoff_multiplier * std::numeric_limits<double>::epsilon() * std::max(magnitude(k, k), noise_floor);
        noise_floor *= coordinate_scale * coordinate_scale;
    }
    return tolerances;
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

    if (abscissas.size() == 2)
    {
        // A two-node candidate has only one pairwise distance, so normalizing it
        // by the maximum pairwise distance (as done below for three or more
        // nodes) is always exactly 1.0 and cannot discriminate separation.
        // Normalize by the larger abscissa magnitude instead, per the documented
        // two-node contract.
        const double distance = std::abs(abscissas[0] - abscissas[1]);
        const double magnitude = std::max(std::abs(abscissas[0]), std::abs(abscissas[1]));
        if (!std::isfinite(distance) || !std::isfinite(magnitude) || magnitude <= 0.0)
        {
            return false;
        }
        return distance / magnitude > options.eabs;
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

[[nodiscard]] WheelerAttemptResult attempt_quadrature_count(const MomentSet& moments, const std::size_t max_nodes, const WheelerOptions& options,
                                                            const double coordinate_scale)
{
    std::size_t n = max_nodes;
    RecurrenceCoefficients coefficients = compute_recurrence(moments, n);

    if (options.adaptive)
    {
        const std::vector<double> sigma_tolerances = adaptive_sigma_tolerances(moments, coefficients, max_nodes, coordinate_scale);
        for (std::size_t k = max_nodes; k > 1; --k)
        {
            const double sigma_tolerance = sigma_tolerances[k];
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

        // The sigma recursion reads only the moments, so the coefficients for a
        // reduced count are an exact prefix of the ones already computed at
        // max_nodes. Recompute only when the scan above actually reduced n.
        if (n != max_nodes)
        {
            coefficients = compute_recurrence(moments, n);
        }
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
        EigenDecomposition eigen = internal::symmetric_jacobi_eigen(make_jacobi_matrix(coefficients.a, coefficients.b, n1), {});
        if (!eigen.converged)
        {
            return {one_node_result(moments, QuadratureStatus::not_realizable, "Jacobi matrix eigen solve did not converge."), false};
        }
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

[[nodiscard]] Quadrature1D compute_quadrature_at_scale(const MomentSet& moments, const std::size_t max_nodes, const WheelerOptions& options,
                                                       const double coordinate_scale)
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
            WheelerAttemptResult attempt = attempt_quadrature_count(moments, candidate_count, options, coordinate_scale);
            if (attempt.quadrature.ok())
            {
                if (candidate_count < max_nodes && attempt.quadrature.weights.size() == candidate_count)
                {
                    internal::performance::record_wheeler_adaptive_reduction();
                }
                return std::move(attempt.quadrature);
            }
            if (!attempt.negative_sigma_diagonal)
            {
                return std::move(attempt.quadrature);
            }
        }

        internal::performance::record_wheeler_adaptive_reduction();
        return one_node_result(moments, QuadratureStatus::success, "");
    }

    return attempt_quadrature_count(moments, max_nodes, options, coordinate_scale).quadrature;
}

} // namespace

Quadrature1D compute_quadrature(const MomentSet& moments, const std::size_t max_nodes, const WheelerOptions& options)
{
    return compute_quadrature_at_scale(moments, max_nodes, options, 0.0);
}

namespace internal
{

Quadrature1D compute_conditional_quadrature(const MomentSet& moments, const std::size_t max_nodes, const WheelerOptions& options,
                                            const double coordinate_scale)
{
    return compute_quadrature_at_scale(moments, max_nodes, options, coordinate_scale);
}

BoundedQuadratureSelection select_bounded_quadrature(const std::vector<double>& moment_sequence, const BoundedQuadratureSelectionRequest& request,
                                                     const AbscissaBounds interval, const WheelerOptions& wheeler_options,
                                                     const double coordinate_scale)
{
    // Both selection branches below assume at least one requested node. The
    // adaptive branch counts candidates down to one, so a zero initial count
    // skips its loop entirely and would leave the default-constructed selection,
    // whose quadrature still carries the default success status, as the result.
    // Rejecting the malformed request here keeps the bounded path's status for a
    // zero node count the same as the unbounded Wheeler path's.
    if (request.initial_node_count == 0 || request.configured_node_count == 0)
    {
        BoundedQuadratureSelection rejected;
        rejected.quadrature.status = QuadratureStatus::invalid_input;
        rejected.quadrature.message = "Bounded-support quadrature requires a positive node count.";
        return rejected;
    }
    // No candidate count exceeds the larger of the initial and configured
    // counts, so one check covers every candidate. Too few moments are rejected
    // exactly as the unbounded Wheeler path rejects them, not answered with the
    // rule the available moments happen to support.
    const std::optional<std::size_t> required_moment_count =
        internal::checked_size_multiply(2, std::max(request.initial_node_count, request.configured_node_count));
    if (!required_moment_count.has_value() || moment_sequence.size() < *required_moment_count)
    {
        BoundedQuadratureSelection rejected;
        rejected.quadrature.status = QuadratureStatus::invalid_input;
        rejected.quadrature.message = required_moment_count.has_value() ? "Moment sequence has fewer than 2*n entries."
                                                                        : "max_nodes is too large to derive the required moment count.";
        return rejected;
    }

    const auto evaluate = [&](const std::size_t node_count, std::string attempt_type)
    {
        BoundedQuadratureCandidate candidate;
        candidate.node_count = node_count;
        candidate.attempt_type = std::move(attempt_type);
        if (!std::isfinite(interval.lower) || !std::isfinite(interval.upper) || !(interval.lower < interval.upper))
        {
            candidate.quadrature.status = QuadratureStatus::invalid_input;
            candidate.message = "Bounded-support interval must have finite lower < upper values.";
            return candidate;
        }

        candidate.inversion_attempted = true;
        candidate.quadrature = compute_quadrature_at_scale(
            MomentSet(std::vector<double>(moment_sequence.begin(), moment_sequence.begin() + static_cast<std::ptrdiff_t>(2 * node_count))),
            node_count, wheeler_options, coordinate_scale);
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

        // A node on a bound comes back off it by the rounding of the moments,
        // amplified by the inversion, which scales with the nodes. Every node
        // inside the interval is bounded by max(|lower|, |upper|), so a tolerance
        // relative to that magnitude gives the same decision in any unit.
        const double coordinate_support_tolerance = 1.0e-12 * std::max(std::abs(interval.lower), std::abs(interval.upper));
        const bool nodes_inside_interval = std::all_of(candidate.quadrature.abscissas.begin(), candidate.quadrature.abscissas.end(),
                                                       [interval, coordinate_support_tolerance](const double node)
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
        const bool accepted = fixed.accepted;
        selection.quadrature = fixed.quadrature;
        selection.selected_node_count_ceiling = request.configured_node_count;
        selection.diagnostic_candidate = std::move(fixed);
        selection.accepted_fallback_node_count = accepted ? request.configured_node_count : 0;
        if (!accepted)
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
