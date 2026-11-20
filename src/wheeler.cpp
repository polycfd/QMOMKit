#include "qmomkit/wheeler.hpp"

#include "internal/bounded_quadrature.hpp"
#include "internal/checked_size.hpp"
#include "internal/generalized_quadrature.hpp"
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
// determine, left zero. A continued recurrence, such as normal_continuation
// returns, holds only a and b, with an empty diagonal.
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

// Whether every abscissa lies in the closed interval, up to a margin for
// rounding. A margin relative to the largest abscissa magnitude gives the
// same decision in any unit, whatever the bounds. The abscissas are not empty
// and in ascending order, so the largest magnitude is that of the first or
// the last abscissa.
[[nodiscard]] bool abscissas_inside_bounds(const std::vector<double>& abscissas, const AbscissaBounds bounds)
{
    const double margin = 1.0e-12 * std::max(std::abs(abscissas.front()), std::abs(abscissas.back()));
    return std::all_of(abscissas.begin(), abscissas.end(), [bounds, margin](const double node)
                       { return std::isfinite(node) && node >= bounds.lower - margin && node <= bounds.upper + margin; });
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

// A fixed count of node_count nodes, one per leading entry of the recurrence,
// must be supported by the moments. A sigma(k, k) that is zero within its
// rounding bound puts the sequence on the boundary at k nodes, and any larger
// quadrature would carry nodes made of rounding. The first such k is reached
// while the coefficients before it are finite; a recurrence that overflows
// instead is reported as non-finite. Every fixed-mode reconstruction makes
// these checks before its eigensolve.
[[nodiscard]] std::optional<Quadrature1D> fixed_count_failure(const MomentSet& moments, const Recurrence& recurrence, const std::size_t node_count)
{
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
    return recurrence_failure(moments, recurrence, node_count);
}

[[nodiscard]] Quadrature1D fixed_quadrature(const MomentSet& moments, const Recurrence& recurrence)
{
    const std::size_t node_count = recurrence.diagonal.size();
    if (const std::optional<Quadrature1D> failure = fixed_count_failure(moments, recurrence, node_count))
    {
        return *failure;
    }
    return gauss_rule(moments, recurrence, node_count);
}

// The input checks of compute_quadrature_at_scale for a reconstruction from
// m_0, ..., m_{2N}, N = node_count, with messages naming node_count and no
// one-node shortcut. Once they pass, 2N + 1 fits in std::size_t, since 2N is
// even and does. A change to these checks belongs in
// compute_quadrature_at_scale as well.
[[nodiscard]] std::optional<Quadrature1D> input_failure_through_m2n(const MomentSet& moments, const std::size_t node_count,
                                                                    const WheelerOptions& options)
{
    if (const std::optional<Quadrature1D> failure = option_failure(options))
    {
        return failure;
    }
    if (node_count == 0)
    {
        return Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "node_count must be positive."};
    }
    const std::optional<std::size_t> even_moment_count = internal::checked_size_multiply(2, node_count);
    if (!even_moment_count.has_value())
    {
        return Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "node_count is too large to derive the required moment count."};
    }
    const std::size_t moment_count = *even_moment_count + 1;
    if (moments.size() < moment_count)
    {
        return Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "Moment sequence has fewer than 2*n+1 entries."};
    }
    return moment_failure(moments, moment_count, options);
}

// The normal continuation of GQMOM (Fox, Laurent & Passalacqua 2023,
// Definition 1, eqs. 21 and 22).
[[nodiscard]] Recurrence normal_continuation(const Recurrence& recurrence, const std::size_t k, const std::size_t continued_count, const double nu)
{
    Recurrence continued{std::vector<double>(continued_count, 0.0), std::vector<double>(continued_count, 0.0), {}};
    const auto known = static_cast<std::ptrdiff_t>(k);
    std::copy(recurrence.a.begin(), recurrence.a.begin() + known, continued.a.begin());
    std::copy(recurrence.b.begin(), recurrence.b.begin() + known + 1, continued.b.begin());
    const double mean = std::accumulate(recurrence.a.begin(), recurrence.a.begin() + known, 0.0) / static_cast<double>(k);
    for (std::size_t i = k; i < continued_count; ++i)
    {
        continued.a[i] = mean;
    }
    for (std::size_t i = k + 1; i < continued_count; ++i)
    {
        continued.b[i] = std::pow(static_cast<double>(i) / static_cast<double>(k), nu) * recurrence.b[k];
    }
    return continued;
}

// The zeta form of the recurrence of k nodes, shifted by lower (Fox, Laurent &
// Passalacqua 2023, eq. 23): zeta[i] holds zeta_i for i = 1, ..., 2k, in a
// vector of the size the continuation to N_G nodes fills. zeta[0] is unused.
[[nodiscard]] std::vector<double> known_zetas(const Recurrence& recurrence, const std::size_t k, const std::size_t continued_count,
                                              const double lower)
{
    std::vector<double> zeta((2 * continued_count) + 1, 0.0);
    zeta[1] = recurrence.a[0] - lower;
    for (std::size_t j = 1; j <= k; ++j)
    {
        zeta[2 * j] = recurrence.b[j] / zeta[(2 * j) - 1];
        if (j < k)
        {
            zeta[(2 * j) + 1] = (recurrence.a[j] - lower) - zeta[2 * j];
        }
    }
    return zeta;
}

// The recurrence of N_G nodes from that of k nodes and the continued zeta (Fox,
// Laurent & Passalacqua 2023, eq. 26): a_i = zeta_{2i} + zeta_{2i+1}, shifted
// back by lower, for i >= k, and b_i = zeta_{2i-1} zeta_{2i} for i > k.
[[nodiscard]] Recurrence recurrence_from_zetas(const Recurrence& recurrence, const std::vector<double>& zeta, const std::size_t k,
                                               const std::size_t continued_count, const double lower)
{
    Recurrence continued{std::vector<double>(continued_count, 0.0), std::vector<double>(continued_count, 0.0), {}};
    const auto known = static_cast<std::ptrdiff_t>(k);
    std::copy(recurrence.a.begin(), recurrence.a.begin() + known, continued.a.begin());
    std::copy(recurrence.b.begin(), recurrence.b.begin() + known + 1, continued.b.begin());
    for (std::size_t i = k; i < continued_count; ++i)
    {
        continued.a[i] = zeta[2 * i] + zeta[(2 * i) + 1] + lower;
    }
    for (std::size_t i = k + 1; i < continued_count; ++i)
    {
        continued.b[i] = zeta[(2 * i) - 1] * zeta[2 * i];
    }
    return continued;
}

// The gamma continuation of GQMOM (Fox, Laurent & Passalacqua 2023,
// Definition 2, eqs. 23, 25 and 26).
[[nodiscard]] Recurrence gamma_continuation(const Recurrence& recurrence, const std::size_t k, const std::size_t continued_count, const double lower)
{
    std::vector<double> zeta = known_zetas(recurrence, k, continued_count, lower);
    const double alpha = (zeta[1] / zeta[2]) - 1.0;
    for (std::size_t i = k + 1; i <= continued_count; ++i)
    {
        zeta[(2 * i) - 1] = (static_cast<double>(i) + alpha) / (static_cast<double>(k) + alpha) * zeta[(2 * k) - 1];
        zeta[2 * i] = static_cast<double>(i) / static_cast<double>(k) * zeta[2 * k];
    }
    return recurrence_from_zetas(recurrence, zeta, k, continued_count, lower);
}

// The lognormal continuation of GQMOM (Fox, Laurent & Passalacqua 2023,
// Definition 3, eqs. 23, 26, 29 and 30). eta^2 = M_2 M_0 / M_1^2 of the
// shifted moments is 1 + zeta_2 / zeta_1.
[[nodiscard]] Recurrence lognormal_continuation(const Recurrence& recurrence, const std::size_t k, const std::size_t continued_count,
                                                const double lower)
{
    std::vector<double> zeta = known_zetas(recurrence, k, continued_count, lower);
    const double eta_squared = 1.0 + (zeta[2] / zeta[1]);
    const double base_factor = std::pow(eta_squared, static_cast<double>(k)) - 1.0;
    for (std::size_t i = k + 1; i <= continued_count; ++i)
    {
        const auto steps = static_cast<double>(i - k);
        zeta[(2 * i) - 1] = std::pow(eta_squared, 2.0 * steps) * zeta[(2 * k) - 1];
        zeta[2 * i] = std::pow(eta_squared, steps) * (std::pow(eta_squared, static_cast<double>(i)) - 1.0) / base_factor * zeta[2 * k];
    }
    return recurrence_from_zetas(recurrence, zeta, k, continued_count, lower);
}

// The beta continuation of GQMOM (Fox, Laurent & Passalacqua 2023,
// Definition 4, eqs. 23, 26 and 31 to 35) on the interval of the bounds.
// Dividing the shifted zeta by the width of the interval maps the recurrence
// to (0, 1), where zeta_1 = p_1 and zeta_i = p_i (1 - p_{i-1}) with the
// canonical moments p_i (eq. 31). alpha and beta from p_1 and p_2 (eq. 33) give
// the canonical moments p^J_i of the beta number density (eq. 32), and every
// continued p_i follows from the last known p of the same parity (eqs. 34 and
// 35). The two cases of eqs. 34 and 35 agree where the choice between them
// changes, so rounding near that point does not change the result
// discontinuously.
[[nodiscard]] Recurrence beta_continuation(const Recurrence& recurrence, const std::size_t k, const std::size_t continued_count,
                                           const AbscissaBounds bounds)
{
    std::vector<double> zeta = known_zetas(recurrence, k, continued_count, bounds.lower);
    const double width = bounds.upper - bounds.lower;
    std::vector<double> p(zeta.size(), 0.0);
    for (std::size_t i = 1; i <= 2 * k; ++i)
    {
        p[i] = zeta[i] / width / (1.0 - p[i - 1]);
    }
    const double alpha = (1.0 - p[1] - (2.0 * p[2]) + (p[1] * p[2])) / p[2];
    const double beta = (p[1] - p[2] - (p[1] * p[2])) / p[2];
    const auto jacobi = [alpha, beta](const std::size_t j)
    {
        const auto i = static_cast<double>((j + 1) / 2);
        return j % 2 == 1 ? (beta + i) / ((2.0 * i) + alpha + beta) : i / ((2.0 * i) + 1.0 + alpha + beta);
    };
    for (std::size_t j = (2 * k) + 1; j <= 2 * continued_count; ++j)
    {
        const std::size_t last = j % 2 == 1 ? (2 * k) - 1 : 2 * k;
        const double p_jacobi = jacobi(j);
        const double last_jacobi = jacobi(last);
        if (p[last] <= last_jacobi || last_jacobi >= p_jacobi)
        {
            p[j] = p_jacobi / last_jacobi * p[last];
        }
        else
        {
            p[j] = ((p[last] * (1.0 - p_jacobi)) + p_jacobi - last_jacobi) / (1.0 - last_jacobi);
        }
        zeta[j] = width * p[j] * (1.0 - p[j - 1]);
    }
    return recurrence_from_zetas(recurrence, zeta, k, continued_count, bounds.lower);
}

[[nodiscard]] bool is_supported_continuation(const GeneralizedContinuation continuation)
{
    switch (continuation)
    {
    case GeneralizedContinuation::normal:
    case GeneralizedContinuation::gamma:
    case GeneralizedContinuation::lognormal:
    case GeneralizedContinuation::beta:
        return true;
    }
    return false;
}

// The HyQMOM closure for a_n, the mean of a_0, ..., a_{n-1} (Fox & Laurent
// 2022, eq. 4.16, eq. 4.13 in standardized moments).
[[nodiscard]] double hyperbolic_closure(const std::vector<double>& a, const std::size_t n)
{
    const auto end = static_cast<std::ptrdiff_t>(n);
    return std::accumulate(a.begin(), a.begin() + end, 0.0) / static_cast<double>(n);
}

// m_{2n+1} of the HyQMOM quadrature, the Gauss quadrature of J_{n+1}, which
// the closure determines (Fox & Laurent 2022, Section 7): m_0 times entry
// (0, 0) of J_{n+1}^{2n+1}. A diagonal scaling that leaves e_1 unchanged takes
// J_{n+1} to the tridiagonal matrix with the same diagonal, b_1, ..., b_n
// above it and ones below it, so the powers act on e_1 through that matrix,
// without square roots.
[[nodiscard]] double hyperbolic_closure_moment(const double zeroth_moment, const Recurrence& recurrence, const std::size_t n)
{
    std::vector<double> power(n + 1, 0.0);
    power[0] = 1.0;
    std::vector<double> next(n + 1, 0.0);
    for (std::size_t order = 1; order <= (2 * n) + 1; ++order)
    {
        for (std::size_t i = 0; i <= n; ++i)
        {
            next[i] = recurrence.a[i] * power[i];
            if (i > 0)
            {
                next[i] += power[i - 1];
            }
            if (i < n)
            {
                next[i] += recurrence.b[i + 1] * power[i + 1];
            }
        }
        power.swap(next);
    }
    return zeroth_moment * power[0];
}

// The hyperbolic reconstruction of m_0, ..., m_{2n}, n = node_count, up to its
// last step. For moments that give the n + 1 nodes of J_{n+1}, outcome is
// empty and recurrence holds a_0, ..., a_n, the last the closure, and b_1,
// ..., b_n. Otherwise outcome is the result of the reconstruction: a failure,
// or in adaptive mode the classical quadrature of moments on or below the
// boundary.
struct HyperbolicRecurrence
{
    std::optional<Quadrature1D> outcome;
    Recurrence recurrence;
};

[[nodiscard]] HyperbolicRecurrence hyperbolic_recurrence(const MomentSet& moments, const std::size_t node_count, const WheelerOptions& options)
{
    const auto finished = [](Quadrature1D quadrature) { return HyperbolicRecurrence{std::move(quadrature), {}}; };
    if (const std::optional<Quadrature1D> failure = input_failure_through_m2n(moments, node_count, options))
    {
        return finished(*failure);
    }
    const std::size_t moment_count = (2 * node_count) + 1;

    // The recurrence of m_0, ..., m_{2n} at coordinate scale zero, with the
    // closure for a_n. Its n + 1 entries give the Jacobi matrix J_{n+1}.
    Recurrence recurrence = compute_recurrence(moments, moment_count, 0.0);
    recurrence.a[node_count] = hyperbolic_closure(recurrence.a, node_count);

    // In adaptive mode, moments on or below the boundary take the classical
    // selection from the first k with sigma(k, k) <= tau_k, which is HyQMOM's
    // own limit there. Interior moments give the quadrature of J_{n+1} without
    // the acceptance tests: the rounding test already excludes nodes made of
    // rounding, and near the boundary the extra node tends to a_n, with a
    // weight that vanishes on the boundary. Fixed mode requires the n + 1
    // nodes.
    if (options.adaptive)
    {
        const std::size_t usable_count = supported_count(recurrence);
        if (usable_count <= node_count)
        {
            return finished(adaptive_selection(moments, recurrence, usable_count, options));
        }
    }
    std::optional<Quadrature1D> failure =
        options.adaptive ? recurrence_failure(moments, recurrence, node_count + 1) : fixed_count_failure(moments, recurrence, node_count + 1);
    return HyperbolicRecurrence{std::move(failure), std::move(recurrence)};
}

// The generalized reconstruction of m_0, ..., m_{2N}, N = node_count, up to its
// continuation. k is the base node count found by the classical selection, 0
// if none was found. For moments in the interior of moment space, outcome is
// empty and recurrence holds the recurrence of m_0, ..., m_{2N}. Otherwise
// outcome is the result of the reconstruction: a failure, or the classical
// k-node quadrature of moments on the boundary or, in adaptive mode, below it.
struct GeneralizedBase
{
    std::optional<internal::GeneralizedQuadrature> outcome;
    Recurrence recurrence;
    std::size_t k = 0;
};

[[nodiscard]] GeneralizedBase generalized_base(const MomentSet& moments, const std::size_t node_count, const WheelerOptions& options,
                                               const GeneralizedOptions& generalized)
{
    internal::performance::record_wheeler_call();
    const auto finished = [](Quadrature1D quadrature, const std::size_t k)
    { return GeneralizedBase{internal::GeneralizedQuadrature{std::move(quadrature), k}, {}, k}; };
    if (generalized.continuation == GeneralizedContinuation::normal && !std::isfinite(generalized.nu))
    {
        return finished(Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "The generalized parameter nu must be finite."}, 0);
    }
    if (generalized.node_count <= node_count)
    {
        return finished(Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "The generalized node count must exceed node_count."}, 0);
    }
    if (!is_supported_continuation(generalized.continuation))
    {
        return finished(Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "Unsupported generalized continuation."}, 0);
    }
    if (const std::optional<Quadrature1D> failure = input_failure_through_m2n(moments, node_count, options))
    {
        return finished(*failure, 0);
    }

    // The recurrence of m_0, ..., m_{2N} at coordinate scale zero. Its first N
    // entries are those of m_0, ..., m_{2N-1}, so the classical quadrature of
    // k <= N nodes below is the one compute_quadrature(moments, N, options)
    // returns, and entry N adds sigma(N, N) with tau_N.
    Recurrence recurrence = compute_recurrence(moments, (2 * node_count) + 1, 0.0);
    Quadrature1D classical;
    if (options.adaptive)
    {
        classical = adaptive_selection(moments, recurrence, std::min(supported_count(recurrence), node_count), options);
    }
    else
    {
        if (const std::optional<Quadrature1D> failure = fixed_count_failure(moments, recurrence, node_count))
        {
            return finished(*failure, 0);
        }
        classical = gauss_rule(moments, recurrence, node_count);
    }
    if (!classical.ok())
    {
        return finished(std::move(classical), 0);
    }

    // sigma(k, k) = m0 b_1 ... b_k. Within its rounding bound, b_k is zero, the
    // moments lie on the boundary of moment space, and GQMOM is the classical
    // k-node quadrature. Below the bound, b_k is negative and the moments lie
    // below the boundary. Fixed mode refuses these moments, and a recurrence
    // that overflowed at k, while adaptive mode keeps the classical selection,
    // which stopped at this k: an overflowed sigma(k, k) of +-inf has an
    // infinite bound and so lies within it. A NaN does not, and its b_k fails
    // the checks of the continued recurrence.
    const std::size_t k = classical.weights.size();
    const SigmaDiagonal& diagonal = recurrence.diagonal[k];
    if (!options.adaptive && (!std::isfinite(diagonal.sigma) || !std::isfinite(diagonal.tolerance)))
    {
        return finished(one_node_result(moments, QuadratureStatus::not_realizable, "Wheeler recurrence produced non-finite values."), 0);
    }
    if (!options.adaptive && diagonal.sigma < -diagonal.tolerance)
    {
        return finished(negative_b_result(moments, recurrence.b[k]), 0);
    }
    if (diagonal.sigma <= diagonal.tolerance)
    {
        return finished(std::move(classical), k);
    }
    return GeneralizedBase{std::nullopt, std::move(recurrence), k};
}

// The result of a reconstruction with exact_base_node_count whose classical
// selection found fewer base nodes than requested.
[[nodiscard]] internal::GeneralizedQuadrature fewer_base_nodes_than_requested()
{
    return internal::GeneralizedQuadrature{
        Quadrature1D{{}, {}, QuadratureStatus::not_realizable, "The generalized reconstruction found fewer base nodes than requested."}, 0};
}

// The continuation of the recurrence of a generalized base without outcome to
// its N_G = generalized.node_count nodes. The gamma and lognormal continuations
// read the lower bound and the beta continuation both bounds. The normal
// continuation reads neither.
[[nodiscard]] internal::GeneralizedQuadrature continued_generalized_quadrature(const MomentSet& moments, const GeneralizedBase& base,
                                                                               const GeneralizedOptions& generalized, const AbscissaBounds bounds)
{
    // The continuation was checked before the moments were read, so
    // is_supported_continuation must accept exactly the cases of this switch.
    Recurrence continued;
    switch (generalized.continuation)
    {
    case GeneralizedContinuation::normal:
        continued = normal_continuation(base.recurrence, base.k, generalized.node_count, generalized.nu);
        break;
    case GeneralizedContinuation::gamma:
        continued = gamma_continuation(base.recurrence, base.k, generalized.node_count, bounds.lower);
        break;
    case GeneralizedContinuation::lognormal:
        continued = lognormal_continuation(base.recurrence, base.k, generalized.node_count, bounds.lower);
        break;
    case GeneralizedContinuation::beta:
        continued = beta_continuation(base.recurrence, base.k, generalized.node_count, bounds);
        break;
    }
    if (const std::optional<Quadrature1D> failure = recurrence_failure(moments, continued, generalized.node_count))
    {
        return internal::GeneralizedQuadrature{*failure, 0};
    }
    Quadrature1D quadrature = gauss_rule(moments, continued, generalized.node_count);
    const std::size_t base_node_count = quadrature.ok() ? base.k : 0;
    return internal::GeneralizedQuadrature{std::move(quadrature), base_node_count};
}

[[nodiscard]] Quadrature1D compute_quadrature_at_scale(const MomentSet& moments, const std::size_t max_nodes, const WheelerOptions& options,
                                                       const double coordinate_scale)
{
    internal::performance::record_wheeler_call();
    // Wheeler inversion path:
    // raw moments -> recurrence coefficients -> Jacobi matrix -> eigenproblem
    // -> quadrature abscissas and weights.
    // input_failure_through_m2n repeats these input checks, and
    // compute_wheeler_recurrence those other than the option and cutoff
    // checks; a change to them belongs there as well.
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
    HyperbolicRecurrence reconstruction = hyperbolic_recurrence(moments, node_count, options);
    if (reconstruction.outcome.has_value())
    {
        return std::move(*reconstruction.outcome);
    }
    return gauss_rule(moments, reconstruction.recurrence, node_count + 1);
}

GeneralizedQuadrature compute_generalized_quadrature(const MomentSet& moments, const std::size_t node_count, const WheelerOptions& options,
                                                     const GeneralizedOptions& generalized, const bool exact_base_node_count)
{
    GeneralizedBase base = generalized_base(moments, node_count, options, generalized);
    if (exact_base_node_count && base.k != 0 && base.k < node_count)
    {
        return fewer_base_nodes_than_requested();
    }
    if (base.outcome.has_value())
    {
        return std::move(*base.outcome);
    }
    return continued_generalized_quadrature(moments, base, generalized, AbscissaBounds{});
}

GeneralizedQuadrature compute_bounded_generalized_quadrature(const MomentSet& moments, const std::size_t node_count, const WheelerOptions& options,
                                                             const GeneralizedOptions& generalized, const AbscissaBounds bounds,
                                                             const bool exact_base_node_count)
{
    if (!std::isfinite(bounds.lower) || !std::isfinite(bounds.upper) || !(bounds.lower < bounds.upper))
    {
        return GeneralizedQuadrature{
            Quadrature1D{{}, {}, QuadratureStatus::invalid_input, "Bounded-support interval must have finite lower < upper values."}, 0};
    }
    // As for bounded classical QMOM, a reconstruction is accepted when it
    // succeeds with every node inside the bounds. Fixed mode tries node_count
    // only. When a reconstruction is rejected, adaptive mode tries again with
    // one node fewer than the base node count k of the rejected
    // reconstruction, or than the node count tried if k was not found, down to
    // one node. The continuation runs once per base node count. Invalid input
    // does not depend on the node count and is returned at once. With
    // exact_base_node_count, the base node count must be node_count, so there
    // is a single attempt.
    std::size_t count = node_count;
    while (true)
    {
        GeneralizedBase base = generalized_base(moments, count, options, generalized);
        if (exact_base_node_count && base.k != 0 && base.k < count)
        {
            return fewer_base_nodes_than_requested();
        }
        GeneralizedQuadrature candidate =
            base.outcome.has_value() ? std::move(*base.outcome) : continued_generalized_quadrature(moments, base, generalized, bounds);
        Quadrature1D& quadrature = candidate.quadrature;
        if (quadrature.status == QuadratureStatus::invalid_input)
        {
            return candidate;
        }
        if (quadrature.ok() && !abscissas_inside_bounds(quadrature.abscissas, bounds))
        {
            quadrature.status = QuadratureStatus::not_realizable;
            quadrature.message = "Bounded-support Wheeler inversion returned nodes outside the configured interval.";
            candidate.base_node_count = 0;
        }
        count = (base.k > 0 ? base.k : count) - 1;
        if (quadrature.ok() || !options.adaptive || exact_base_node_count || count == 0)
        {
            return candidate;
        }
    }
}

HyperbolicCharacteristicSpeeds compute_hyperbolic_characteristic_speeds(const MomentSet& moments, const std::size_t node_count,
                                                                        const WheelerOptions& options)
{
    // One Wheeler inversion of the moments, recorded as for
    // compute_hyperbolic_quadrature.
    performance::record_wheeler_call();
    const auto failure = [](const QuadratureStatus status, std::string message)
    { return HyperbolicCharacteristicSpeeds{{}, {}, 0.0, status, std::move(message)}; };

    // The speeds are those of the quadrature that compute_hyperbolic_quadrature
    // reconstructs from the same moments, decided on the same recurrence. Each
    // speed is paired with its weight, so that sorting the speeds keeps the
    // pairs together.
    HyperbolicRecurrence reconstruction = hyperbolic_recurrence(moments, node_count, options);
    std::vector<std::pair<double, double>> weighted_speeds;
    double closure_moment = 0.0;
    if (reconstruction.outcome.has_value())
    {
        const Quadrature1D& quadrature = *reconstruction.outcome;
        if (!quadrature.ok())
        {
            return failure(quadrature.status, quadrature.message);
        }

        // A classical quadrature of k <= n nodes has the speeds of the k-point
        // system on its boundary, beta_k = 0, where R_{k+1} = (x - a_k) Q_k
        // (Fox & Laurent 2022, Remark 4.9, for k = n): the roots of Q_k, its k
        // abscissas, each twice, and the closure a_k, the trace of J_k over k
        // and so the arithmetic mean of the abscissas. The weights are the
        // limit of those of the k-point system as beta_k tends to zero. The
        // copy of an abscissa from J_k has (k + 1) / (2k + 1) of its weight,
        // the copy from K_{k+1} has k / (2k + 1), and the mean, which
        // decouples from K_{k+1}, has weight zero. The eigenvalue of K_{k+1}
        // that tends to an abscissa below the mean approaches it from below,
        // and one that tends to an abscissa above the mean from above (Fox &
        // Laurent 2022, Theorem 4.8), so the copy from K_{k+1} comes first
        // below the mean. At an abscissa equal to the mean, always the case
        // for one node, two eigenvalues of K_{k+1} approach it from both
        // sides and share its weight, so only the sum of the weights at the
        // three equal speeds is continuous there.
        const std::vector<double>& abscissas = quadrature.abscissas;
        const auto k = static_cast<double>(abscissas.size());
        const double mean = std::accumulate(abscissas.begin(), abscissas.end(), 0.0) / k;
        const double j_copy_fraction = (k + 1.0) / ((2.0 * k) + 1.0);
        const double k_copy_fraction = k / ((2.0 * k) + 1.0);
        weighted_speeds.reserve((2 * abscissas.size()) + 1);
        for (std::size_t i = 0; i < abscissas.size(); ++i)
        {
            const double j_copy_weight = j_copy_fraction * quadrature.weights[i];
            const double k_copy_weight = k_copy_fraction * quadrature.weights[i];
            if (abscissas[i] < mean)
            {
                weighted_speeds.emplace_back(abscissas[i], k_copy_weight);
                weighted_speeds.emplace_back(abscissas[i], j_copy_weight);
            }
            else
            {
                weighted_speeds.emplace_back(abscissas[i], j_copy_weight);
                weighted_speeds.emplace_back(abscissas[i], k_copy_weight);
            }
            closure_moment += quadrature.weights[i] * std::pow(abscissas[i], static_cast<double>((2 * node_count) + 1));
        }
        weighted_speeds.emplace_back(mean, 0.0);
    }
    else
    {
        // n + 1 nodes: J_n uses a_0, ..., a_{n-1} and b_1, ..., b_{n-1}.
        // K_{n+1} adds the closure a_n and b_n, replaced by beta_n =
        // (2n + 1) b_n / n (Fox & Laurent 2022, eqs. 4.13, 4.16 and 4.18).
        // The closure moment needs b_n itself, so it comes first.
        Recurrence& recurrence = reconstruction.recurrence;
        closure_moment = hyperbolic_closure_moment(moments[0], recurrence, node_count);
        const EigenDecomposition j_eigen = symmetric_jacobi_eigen(make_jacobi_matrix(recurrence.a, recurrence.b, node_count), {});
        recurrence.b[node_count] = static_cast<double>((2 * node_count) + 1) * recurrence.b[node_count] / static_cast<double>(node_count);
        const EigenDecomposition k_eigen = symmetric_jacobi_eigen(make_jacobi_matrix(recurrence.a, recurrence.b, node_count + 1), {});
        if (!j_eigen.converged || !k_eigen.converged)
        {
            return failure(QuadratureStatus::not_realizable, "Jacobi matrix eigen solve did not converge.");
        }

        // The speeds with these weights have the moments m_0, ..., m_{2n+1} of
        // the quadrature of J_{n+1}: the Gauss quadratures of J_n and K_{n+1}
        // share the moments up to m_{2n-1}, and the factors (n + 1) / (2n + 1)
        // and n / (2n + 1) sum to one and match m_{2n} and m_{2n+1}, because
        // beta_n = (2n + 1) b_n / n (docs/developer/qmom_families.md).
        const auto n = static_cast<double>(node_count);
        const std::vector<double> j_weights = compute_weights(moments[0], j_eigen.eigenvectors);
        const std::vector<double> k_weights = compute_weights(moments[0], k_eigen.eigenvectors);
        weighted_speeds.reserve((2 * node_count) + 1);
        for (std::size_t i = 0; i < node_count; ++i)
        {
            weighted_speeds.emplace_back(j_eigen.eigenvalues[i], (n + 1.0) / ((2.0 * n) + 1.0) * j_weights[i]);
        }
        for (std::size_t i = 0; i <= node_count; ++i)
        {
            weighted_speeds.emplace_back(k_eigen.eigenvalues[i], n / ((2.0 * n) + 1.0) * k_weights[i]);
        }
    }

    // A stable sort keeps equal speeds in the order in which they were added.
    std::stable_sort(weighted_speeds.begin(), weighted_speeds.end(), [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });
    std::vector<double> speeds;
    std::vector<double> weights;
    speeds.reserve(weighted_speeds.size());
    weights.reserve(weighted_speeds.size());
    for (const auto& [speed, weight] : weighted_speeds)
    {
        speeds.push_back(speed);
        weights.push_back(weight);
    }
    if (!all_finite(speeds))
    {
        return failure(QuadratureStatus::not_realizable, "Characteristic speeds exceed the range of a double.");
    }
    if (!std::isfinite(closure_moment))
    {
        return failure(QuadratureStatus::not_realizable, "Closure moment exceeds the range of a double.");
    }
    return HyperbolicCharacteristicSpeeds{std::move(speeds), std::move(weights), closure_moment, QuadratureStatus::success, ""};
}

BoundedQuadratureSelection select_bounded_quadrature(const std::vector<double>& moment_sequence, const BoundedQuadratureSelectionRequest& request,
                                                     const AbscissaBounds interval, const WheelerOptions& wheeler_options,
                                                     const double coordinate_scale)
{
    // Both selection branches below assume at least one requested node. The
    // adaptive branch counts candidates down to one, and from a zero initial
    // count its unsigned count would wrap around. Rejecting the malformed
    // request here keeps the bounded path's status for a zero node count the
    // same as the unbounded Wheeler path's.
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
        // An adaptive inversion returns fewer nodes than node_count when the
        // moments carry fewer points or a larger rule fails the rmin or eabs
        // test. Each recurrence entry depends only on the moments up to its
        // own order, so every count from that number of nodes up to node_count
        // gives the same rule, and the bounds alone decide on it.
        const std::size_t inverted_node_count = candidate.quadrature.weights.size();
        if ((request.adaptive ? inverted_node_count > node_count : inverted_node_count != node_count) ||
            candidate.quadrature.abscissas.size() != inverted_node_count)
        {
            candidate.quadrature.status = QuadratureStatus::not_realizable;
            candidate.message = "Bounded-support Wheeler inversion returned an unexpected node count.";
            return candidate;
        }
        candidate.inverted_node_count = inverted_node_count;

        if (!abscissas_inside_bounds(candidate.quadrature.abscissas, interval))
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
    std::size_t candidate_count = first_count;
    while (true)
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
            selection.accepted_fallback_node_count = selection.quadrature.weights.size();
            break;
        }
        if (!rejected_candidate.has_value())
        {
            rejected_candidate = candidate;
        }
        // A rule rejected for its bounds is also the rule of every smaller
        // count down to its own number of nodes, so the next candidate has one
        // node less than the rule.
        const std::size_t next_count = (candidate.inverted_node_count > 0 ? candidate.inverted_node_count : candidate_count) - 1;
        if (next_count == 0)
        {
            selection.quadrature = std::move(candidate.quadrature);
            selection.selected_node_count_ceiling = first_count;
            selection.diagnostic_candidate = std::move(*rejected_candidate);
            selection.quadrature.message = candidate.message;
            return selection;
        }
        candidate_count = next_count;
    }

    if (request.allow_growth && request.adaptive && wheeler_options.adaptive && request.previous_active_node_count.has_value() &&
        selection.quadrature.ok() && selection.quadrature.weights.size() == *request.previous_active_node_count &&
        *request.previous_active_node_count < request.configured_node_count)
    {
        const std::size_t growth_count = *request.previous_active_node_count + 1;
        BoundedQuadratureCandidate growth = evaluate(growth_count, "growth");
        selection.diagnostic_candidate = growth;
        // The trial counts only if it gets the one node more it tried for.
        if (growth.accepted && growth.inverted_node_count == growth_count)
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
