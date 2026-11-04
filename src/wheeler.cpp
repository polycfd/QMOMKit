#include "qmomkit/wheeler.hpp"

#include "internal/bounded_quadrature.hpp"
#include "internal/checked_size.hpp"
#include "internal/hyperbolic_quadrature.hpp"
#include "internal/linear_algebra.hpp"
#include "internal/numeric_utils.hpp"
#include "internal/performance.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <numeric>
#include <optional>
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

struct SigmaDiagonal
{
    double sigma = 0.0;
    double tolerance = 0.0;
};

// Wheeler's recurrence for the moments m_0, ..., m_{2N-1}, indexed as in
// docs/developer/qmom.md: a_0, ..., a_{N-1}, b_0 = 0, b_1, ..., b_{N-1}, and for
// k = 1, ..., N - 1 the diagonal sigma(k, k) with its rounding bound tau_k.
// Entry 0 of the diagonal is unused. The moments m_0, ..., m_{2N} give one
// entry more: b_N and sigma(N, N) with tau_N, and a_N, which they do not
// determine, left zero.
struct Recurrence
{
    std::vector<double> a;
    std::vector<double> b;
    std::vector<SigmaDiagonal> diagonal;
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

// Row k of the sigma table of M = moment_count moments, sigma(k, l) for
// l = k, ..., M - k - 1, depends only on rows k - 1 and k - 2, so the rows are
// computed one at a time from sigma(-1, l) = 0 and sigma(0, l) = m_l. The same
// pass runs the recurrence with every term in absolute value, whose diagonal
// bounds the rounding of sigma(k, k). sigma(k, k) and the bound both scale as
// m0 times the coordinate to the power 2k, so the tests on them depend neither
// on the unit of the coordinate nor on m0, and the bound grows with the
// cancellation of a distribution far from the origin. Conditional moments are
// not data: the solve that produced them leaves rounding relative to the scale
// of their coordinate, not to their own magnitudes, so a positive
// coordinate_scale s adds the floor m0 s^(2k).
//
// a_k needs the entry l = k + 1 of row k. From M = 2N moments that gives rows
// and coefficients up to k = N - 1; from M = 2N + 1 also row N, whose only
// entry is sigma(N, N), with b_N and tau_N. The entries for k < N do not
// depend on m_{2N}.
[[nodiscard]] Recurrence compute_recurrence(const MomentSet& moments, const std::size_t moment_count, const double coordinate_scale)
{
    constexpr double roundoff_multiplier = 100.0;
    const std::size_t size = (moment_count + 1) / 2;
    Recurrence recurrence{std::vector<double>(size, 0.0), std::vector<double>(size, 0.0), std::vector<SigmaDiagonal>(size)};

    // Rows k - 2, k - 1 and k of sigma and of its absolute-value bound,
    // indexed by l, in one workspace.
    const std::size_t width = moment_count;
    std::vector<double> rows(6 * width, 0.0);
    double* sigma_second_previous = rows.data();
    double* sigma_previous = sigma_second_previous + width;
    double* sigma_current = sigma_previous + width;
    double* bound_second_previous = sigma_current + width;
    double* bound_previous = bound_second_previous + width;
    double* bound_current = bound_previous + width;
    for (std::size_t l = 0; l < width; ++l)
    {
        sigma_previous[l] = moments[l];
        bound_previous[l] = std::abs(moments[l]);
    }

    recurrence.a[0] = moments[1] / moments[0];
    recurrence.b[0] = 0.0;
    double noise_floor = std::abs(moments[0]);

    for (std::size_t k = 1; k < size; ++k)
    {
        const double a = recurrence.a[k - 1];
        const double b = recurrence.b[k - 1];
        for (std::size_t l = k; l < width - k; ++l)
        {
            sigma_current[l] = sigma_previous[l + 1] - a * sigma_previous[l] - b * sigma_second_previous[l];
            bound_current[l] = bound_previous[l + 1] + std::abs(a) * bound_previous[l] + std::abs(b) * bound_second_previous[l];
        }

        if (k + 1 < width - k)
        {
            recurrence.a[k] = sigma_current[k + 1] / sigma_current[k] - sigma_previous[k] / sigma_previous[k - 1];
        }
        recurrence.b[k] = sigma_current[k] / sigma_previous[k - 1];
        noise_floor *= coordinate_scale * coordinate_scale;
        recurrence.diagonal[k] = {sigma_current[k],
                                  roundoff_multiplier * std::numeric_limits<double>::epsilon() * std::max(bound_current[k], noise_floor)};

        std::swap(sigma_second_previous, sigma_previous);
        std::swap(sigma_previous, sigma_current);
        std::swap(bound_second_previous, bound_previous);
        std::swap(bound_previous, bound_current);
    }

    return recurrence;
}

[[nodiscard]] bool recurrence_is_finite(const Recurrence& recurrence, const std::size_t n)
{
    const auto finite = [](const double value) { return std::isfinite(value); };
    const auto end = static_cast<std::ptrdiff_t>(n);
    return std::all_of(recurrence.a.begin(), recurrence.a.begin() + end, finite) &&
           std::all_of(recurrence.b.begin(), recurrence.b.begin() + end, finite);
}

[[nodiscard]] double minimum_b(const Recurrence& recurrence, const std::size_t n)
{
    return *std::min_element(recurrence.b.begin(), recurrence.b.begin() + static_cast<std::ptrdiff_t>(n));
}

[[nodiscard]] Quadrature1D negative_b_result(const MomentSet& moments, const double minimum)
{
    std::ostringstream message;
    message << std::scientific << std::setprecision(17) << "Moment sequence has negative recurrence b = " << minimum << '.';
    return one_node_result(moments, QuadratureStatus::not_realizable, message.str());
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

// The n-node Gauss rule of the leading n recurrence coefficients, with its
// abscissas in ascending order.
[[nodiscard]] Quadrature1D gauss_rule(const MomentSet& moments, const Recurrence& recurrence, const std::size_t n)
{
    EigenDecomposition eigen = internal::symmetric_jacobi_eigen(make_jacobi_matrix(recurrence.a, recurrence.b, n), {});
    if (!eigen.converged)
    {
        return one_node_result(moments, QuadratureStatus::not_realizable, "Jacobi matrix eigen solve did not converge.");
    }
    std::vector<double> weights = compute_weights(moments[0], eigen.eigenvectors);
    if (!internal::all_finite(eigen.eigenvalues) || !internal::all_finite(weights))
    {
        return one_node_result(moments, QuadratureStatus::not_realizable, "Jacobi matrix diagonalization produced non-finite values.");
    }
    return Quadrature1D{std::move(eigen.eigenvalues), std::move(weights), QuadratureStatus::success, ""};
}

// Adaptive inversion keeps a rule only when its weights remain significant and
// its abscissas remain separated at the requested tolerances.
[[nodiscard]] bool adaptive_accepts(const Quadrature1D& rule, const WheelerOptions& options)
{
    const std::vector<double>& abscissas = rule.abscissas;
    const std::vector<double>& weights = rule.weights;
    const auto minmax_weight = std::minmax_element(weights.begin(), weights.end());
    if (*minmax_weight.first / *minmax_weight.second <= options.rmin)
    {
        return false;
    }

    if (abscissas.size() == 2)
    {
        // A single distance has ratio one to the span, so a two-node rule
        // normalizes it by the larger abscissa magnitude instead.
        const double distance = std::abs(abscissas[0] - abscissas[1]);
        const double magnitude = std::max(std::abs(abscissas[0]), std::abs(abscissas[1]));
        if (!std::isfinite(distance) || !std::isfinite(magnitude) || magnitude <= 0.0)
        {
            return false;
        }
        return distance / magnitude > options.eabs;
    }

    // The abscissas are sorted and rounded subtraction is monotone, so the
    // smallest computed distance between two abscissas is one between
    // neighbors and the largest is the span.
    double smallest_gap = std::numeric_limits<double>::infinity();
    for (std::size_t i = 1; i < abscissas.size(); ++i)
    {
        smallest_gap = std::min(smallest_gap, abscissas[i] - abscissas[i - 1]);
    }
    const double span = abscissas.back() - abscissas.front();
    if (!std::isfinite(smallest_gap) || !std::isfinite(span) || span <= 0.0)
    {
        return false;
    }
    return smallest_gap / span > options.eabs;
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

[[nodiscard]] std::optional<Quadrature1D> option_failure(const WheelerOptions& options)
{
    if (!std::isfinite(options.rmin) || !std::isfinite(options.eabs) || !std::isfinite(options.cutoff) || options.rmin < 0.0 || options.eabs < 0.0 ||
        options.cutoff <= 0.0)
    {
        return Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "Invalid Wheeler numerical option."};
    }
    return std::nullopt;
}

// The checks on the first moment_count moments, after those on the node count:
// all finite, m0 >= 0 and m0 >= cutoff.
[[nodiscard]] std::optional<Quadrature1D> moment_failure(const MomentSet& moments, const std::size_t moment_count, const WheelerOptions& options)
{
    if (!required_moments_are_finite(moments, moment_count))
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
    return std::nullopt;
}

// The checks before a Gauss rule of node_count nodes: a non-finite recurrence
// or a negative b is not_realizable.
[[nodiscard]] std::optional<Quadrature1D> recurrence_failure(const MomentSet& moments, const Recurrence& recurrence, const std::size_t node_count)
{
    if (!recurrence_is_finite(recurrence, node_count))
    {
        return one_node_result(moments, QuadratureStatus::not_realizable, "Wheeler recurrence produced non-finite values.");
    }
    const double smallest_b = minimum_b(recurrence, node_count);
    if (smallest_b < 0.0)
    {
        return negative_b_result(moments, smallest_b);
    }
    return std::nullopt;
}

// sigma(k, k) = m0 b_1 ... b_k is positive while m_0, ..., m_{2k} are the
// moments of at least k + 1 points, zero when k points carry them, and
// negative when no non-negative density has them. The first k with sigma(k, k)
// at or below its rounding bound therefore caps the usable count at k; without
// one, every entry of the recurrence is usable.
[[nodiscard]] std::size_t supported_count(const Recurrence& recurrence)
{
    for (std::size_t k = 1; k < recurrence.diagonal.size(); ++k)
    {
        if (recurrence.diagonal[k].sigma <= recurrence.diagonal[k].tolerance)
        {
            return k;
        }
    }
    return recurrence.diagonal.size();
}

// From usable_count down, the first Gauss rule that passes the acceptance tests
// is returned, and the one-node rule if none does.
[[nodiscard]] Quadrature1D adaptive_selection(const MomentSet& moments, const Recurrence& recurrence, const std::size_t usable_count,
                                              const WheelerOptions& options)
{
    if (usable_count == 1)
    {
        return one_node_result(moments, QuadratureStatus::success, "");
    }
    if (const std::optional<Quadrature1D> failure = recurrence_failure(moments, recurrence, usable_count))
    {
        return *failure;
    }

    for (std::size_t node_count = usable_count; node_count > 1; --node_count)
    {
        Quadrature1D rule = gauss_rule(moments, recurrence, node_count);
        if (!rule.ok() || adaptive_accepts(rule, options))
        {
            return rule;
        }
    }
    return one_node_result(moments, QuadratureStatus::success, "");
}

// A fixed count, one node per entry of the recurrence, must be supported by
// the moments. A sigma(k, k) that is zero within its rounding bound puts the
// sequence on the boundary at k nodes, and any larger rule would carry nodes
// made of rounding. The first such k is reached while the coefficients before
// it are finite; a recurrence that overflows instead is reported as non-finite.
[[nodiscard]] Quadrature1D fixed_quadrature(const MomentSet& moments, const Recurrence& recurrence)
{
    const std::size_t node_count = recurrence.diagonal.size();
    for (std::size_t k = 1; k < node_count; ++k)
    {
        const double sigma = recurrence.diagonal[k].sigma;
        const double tolerance = recurrence.diagonal[k].tolerance;
        if (std::isfinite(tolerance) && std::abs(sigma) <= tolerance)
        {
            std::ostringstream message;
            message << std::scientific << std::setprecision(17) << "Moment sequence supports at most " << k << (k == 1 ? " node" : " nodes")
                    << ": sigma(" << k << ", " << k << ") = " << sigma << " is zero within the tolerance " << tolerance << '.';
            return one_node_result(moments, QuadratureStatus::not_realizable, message.str());
        }
    }

    if (const std::optional<Quadrature1D> failure = recurrence_failure(moments, recurrence, node_count))
    {
        return *failure;
    }
    return gauss_rule(moments, recurrence, node_count);
}

// The HyQMOM closure for a_n, the mean of a_0, ..., a_{n-1} (Fox & Laurent
// 2022, eq. 4.16, eq. 4.13 in standardized moments).
[[nodiscard]] double hyperbolic_closure(const std::vector<double>& a, const std::size_t n)
{
    const auto end = static_cast<std::ptrdiff_t>(n);
    return std::accumulate(a.begin(), a.begin() + end, 0.0) / static_cast<double>(n);
}

[[nodiscard]] double dot(const std::vector<double>& lhs, const std::vector<double>& rhs)
{
    return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), 0.0);
}

// The recurrence a_0, ..., a_{N-1}, b_1, ..., b_{N-1} of N weighted points,
// that of their moments m_0, ..., m_{2N-1}, by Lanczos on diag(x) from the unit
// vector q_0 with components sqrt(w_i / m_0). Its vector q_k holds
// p_k(x_i) sqrt(w_i / m_0) for the orthonormal polynomials p_k of the points,
// so a_k = q_k^T diag(x) q_k and b_{k+1} is the squared norm of diag(x) q_k
// after its projections onto q_0, ..., q_k are removed. Removing them twice
// keeps the vectors orthogonal in floating point; once is not enough where the
// projections cancel most of diag(x) q_k, as for points of small spread far
// from the origin. Orthogonal transformations avoid the cancellation of the
// moment table near the boundary of moment space. A b_{k+1} that is zero, as
// for points all at zero or coordinates that underflow, or not finite ends the
// recurrence, with the later entries zero; points that coincide elsewhere
// leave a b_{k+1} of rounding size. The diagonal of the Wheeler table is left
// empty.
[[nodiscard]] Recurrence weighted_point_recurrence(const std::vector<double>& abscissas, const std::vector<double>& weights)
{
    const std::size_t count = abscissas.size();
    const double m0 = std::accumulate(weights.begin(), weights.end(), 0.0);
    Recurrence recurrence{std::vector<double>(count, 0.0), std::vector<double>(count, 0.0), {}};
    std::vector<std::vector<double>> basis;
    basis.reserve(count);
    std::vector<double> current(count);
    for (std::size_t i = 0; i < count; ++i)
    {
        current[i] = std::sqrt(weights[i] / m0);
    }
    for (std::size_t k = 0; k < count; ++k)
    {
        basis.push_back(current);
        std::vector<double> next(count);
        for (std::size_t i = 0; i < count; ++i)
        {
            next[i] = abscissas[i] * current[i];
        }
        recurrence.a[k] = dot(current, next);
        if (k + 1 == count)
        {
            break;
        }
        for (int pass = 0; pass < 2; ++pass)
        {
            for (const std::vector<double>& earlier : basis)
            {
                const double projection = dot(earlier, next);
                for (std::size_t i = 0; i < count; ++i)
                {
                    next[i] -= projection * earlier[i];
                }
            }
        }
        const double norm_squared = dot(next, next);
        recurrence.b[k + 1] = norm_squared;
        if (norm_squared == 0.0 || !std::isfinite(norm_squared))
        {
            break;
        }
        const double norm = std::sqrt(norm_squared);
        for (std::size_t i = 0; i < count; ++i)
        {
            current[i] = next[i] / norm;
        }
    }
    return recurrence;
}

[[nodiscard]] Quadrature1D compute_quadrature_at_scale(const MomentSet& moments, const std::size_t max_nodes, const WheelerOptions& options,
                                                       const double coordinate_scale)
{
    internal::performance::record_wheeler_call();
    // Wheeler inversion path:
    // raw moments -> recurrence coefficients -> Jacobi matrix -> eigenproblem
    // -> quadrature abscissas and weights.
    if (const std::optional<Quadrature1D> failure = option_failure(options))
    {
        return *failure;
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
    if (const std::optional<Quadrature1D> failure = moment_failure(moments, *required_moment_count, options))
    {
        return *failure;
    }

    if (max_nodes == 1)
    {
        return one_node_result(moments, QuadratureStatus::success, "");
    }

    const Recurrence recurrence = compute_recurrence(moments, *required_moment_count, coordinate_scale);
    return options.adaptive ? adaptive_selection(moments, recurrence, supported_count(recurrence), options) : fixed_quadrature(moments, recurrence);
}

} // namespace

Quadrature1D compute_quadrature(const MomentSet& moments, const std::size_t max_nodes, const WheelerOptions& options)
{
    return compute_quadrature_at_scale(moments, max_nodes, options, 0.0);
}

WheelerRecurrence compute_wheeler_recurrence(const MomentSet& moments, const std::size_t node_count)
{
    const auto failure = [](const QuadratureStatus status, std::string message)
    {
        WheelerRecurrence result;
        result.status = status;
        result.message = std::move(message);
        return result;
    };

    // The input checks of compute_quadrature_at_scale, in its order and with
    // its messages, except that these name node_count; a change to one set of
    // checks belongs in the other. The option and cutoff checks have no
    // counterpart here, and since the recurrence divides by m0, m0 = 0 is
    // refused instead.
    if (node_count == 0)
    {
        return failure(QuadratureStatus::invalid_input, "node_count must be positive.");
    }
    const std::optional<std::size_t> required_moment_count = internal::checked_size_multiply(2, node_count);
    if (!required_moment_count.has_value())
    {
        return failure(QuadratureStatus::invalid_input, "node_count is too large to derive the required moment count.");
    }
    if (moments.size() < *required_moment_count)
    {
        return failure(QuadratureStatus::invalid_input, "Moment sequence has fewer than 2*n entries.");
    }
    if (!required_moments_are_finite(moments, *required_moment_count))
    {
        return failure(QuadratureStatus::invalid_input, "Moment sequence contains non-finite values.");
    }
    if (moments[0] < 0.0)
    {
        return failure(QuadratureStatus::not_realizable, "Moment sequence is not realizable: m0 < 0.");
    }
    if (moments[0] == 0.0)
    {
        return failure(QuadratureStatus::not_realizable, "Moment sequence has m0 = 0.");
    }

    // The recurrence of the QMOM path, which passes coordinate scale zero.
    // compute_recurrence leaves diagonal entry 0 unused, since no test is made
    // at k = 0; sigma(0, 0) is m0. No record_wheeler_call(): the counter
    // counts inversions, and this function inverts nothing.
    const Recurrence recurrence = compute_recurrence(moments, *required_moment_count, 0.0);
    WheelerRecurrence result;
    result.a = recurrence.a;
    result.b = recurrence.b;
    result.sigma.assign(node_count, 0.0);
    result.tolerance.assign(node_count, 0.0);
    result.sigma[0] = moments[0];
    for (std::size_t k = 1; k < node_count; ++k)
    {
        result.sigma[k] = recurrence.diagonal[k].sigma;
        result.tolerance[k] = recurrence.diagonal[k].tolerance;
    }
    return result;
}

namespace internal
{

Quadrature1D compute_conditional_quadrature(const MomentSet& moments, const std::size_t max_nodes, const WheelerOptions& options,
                                            const double coordinate_scale)
{
    return compute_quadrature_at_scale(moments, max_nodes, options, coordinate_scale);
}

Quadrature1D compute_hyperbolic_quadrature(const MomentSet& moments, const std::size_t node_count, const WheelerOptions& options)
{
    performance::record_wheeler_call();
    // The input checks of compute_quadrature_at_scale for m_0, ..., m_{2n},
    // with messages naming node_count and no one-node shortcut: n = 1 gives
    // the two nodes of J_2. A change to these checks belongs in
    // compute_quadrature_at_scale as well.
    if (const std::optional<Quadrature1D> failure = option_failure(options))
    {
        return *failure;
    }
    if (node_count == 0)
    {
        return Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "node_count must be positive."};
    }
    const std::optional<std::size_t> even_moment_count = checked_size_multiply(2, node_count);
    if (!even_moment_count.has_value())
    {
        return Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "node_count is too large to derive the required moment count."};
    }
    // 2n is even, so 2n + 1 cannot overflow once 2n does not.
    const std::size_t moment_count = *even_moment_count + 1;
    if (moments.size() < moment_count)
    {
        return Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "Moment sequence has fewer than 2*n+1 entries."};
    }
    if (const std::optional<Quadrature1D> failure = moment_failure(moments, moment_count, options))
    {
        return *failure;
    }

    // The recurrence of m_0, ..., m_{2n} at coordinate scale zero, with the
    // closure for a_n. Its n + 1 entries give the Jacobi matrix J_{n+1}.
    Recurrence recurrence = compute_recurrence(moments, moment_count, 0.0);
    recurrence.a[node_count] = hyperbolic_closure(recurrence.a, node_count);
    if (!options.adaptive)
    {
        return fixed_quadrature(moments, recurrence);
    }

    // Moments on or below the boundary take the classical selection from the
    // first k with sigma(k, k) <= tau_k, which is HyQMOM's own limit there.
    // Interior moments give the rule of J_{n+1} without the acceptance tests:
    // the rounding test already excludes nodes made of rounding, and near the
    // boundary the extra node tends to a_n, with a weight that vanishes on the
    // boundary.
    const std::size_t usable_count = supported_count(recurrence);
    if (usable_count <= node_count)
    {
        return adaptive_selection(moments, recurrence, usable_count, options);
    }
    if (const std::optional<Quadrature1D> failure = recurrence_failure(moments, recurrence, node_count + 1))
    {
        return *failure;
    }
    return gauss_rule(moments, recurrence, node_count + 1);
}

HyperbolicCharacteristicSpeeds compute_hyperbolic_characteristic_speeds(const std::vector<double>& abscissas, const std::vector<double>& weights)
{
    const auto failure = [](std::string message) { return HyperbolicCharacteristicSpeeds{{}, QuadratureStatus::not_realizable, std::move(message)}; };
    const std::size_t node_count = abscissas.size() - 1;
    Recurrence recurrence = weighted_point_recurrence(abscissas, weights);
    recurrence.a[node_count] = hyperbolic_closure(recurrence.a, node_count);
    if (!recurrence_is_finite(recurrence, node_count + 1))
    {
        return failure("Characteristic speed recurrence produced non-finite values.");
    }
    const auto end = static_cast<std::ptrdiff_t>(node_count + 1);
    if (std::any_of(recurrence.b.begin() + 1, recurrence.b.begin() + end, [](const double b) { return b <= 0.0; }))
    {
        return failure("Characteristic speed recurrence of the nodes vanished.");
    }

    // J_n uses b_1, ..., b_{n-1}. K_{n+1} adds b_n, replaced by beta_n =
    // (2n + 1) b_n / n (Fox & Laurent 2022, eqs. 4.13, 4.16 and 4.18).
    const EigenDecomposition j_eigen = symmetric_jacobi_eigen(make_jacobi_matrix(recurrence.a, recurrence.b, node_count), {});
    recurrence.b[node_count] = static_cast<double>((2 * node_count) + 1) * recurrence.b[node_count] / static_cast<double>(node_count);
    const EigenDecomposition k_eigen = symmetric_jacobi_eigen(make_jacobi_matrix(recurrence.a, recurrence.b, node_count + 1), {});
    for (const EigenDecomposition* eigen : {&j_eigen, &k_eigen})
    {
        if (!eigen->converged)
        {
            return failure("Jacobi matrix eigen solve did not converge.");
        }
        if (!all_finite(eigen->eigenvalues))
        {
            return failure("Jacobi matrix diagonalization produced non-finite values.");
        }
    }

    // Both eigenvalue sequences are ascending.
    std::vector<double> speeds((2 * node_count) + 1);
    std::merge(j_eigen.eigenvalues.begin(), j_eigen.eigenvalues.end(), k_eigen.eigenvalues.begin(), k_eigen.eigenvalues.end(), speeds.begin());
    return HyperbolicCharacteristicSpeeds{std::move(speeds), QuadratureStatus::success, ""};
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
    // counts, so one check covers every candidate. Too few moments, or a
    // non-finite one among those required, are rejected with the status and
    // message of the unbounded Wheeler path, not answered with a rule of fewer
    // nodes that never reads them.
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
    if (!std::all_of(moment_sequence.begin(), moment_sequence.begin() + static_cast<std::ptrdiff_t>(*required_moment_count),
                     [](const double value) { return std::isfinite(value); }))
    {
        BoundedQuadratureSelection rejected;
        rejected.quadrature.status = QuadratureStatus::invalid_input;
        rejected.quadrature.message = "Moment sequence contains non-finite values.";
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
