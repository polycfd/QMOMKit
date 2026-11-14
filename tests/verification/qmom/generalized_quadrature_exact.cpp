// GQMOM reconstructions with the normal and gamma continuations, against Gauss
// rules of the continued recurrence computed in exact rational arithmetic and
// diagonalized at 50 digits.
//
// What is verified
//   That compute_qmom(...) with the generalized family returns the N_G nodes and
//   weights of Definition 1 of Fox, Laurent & Passalacqua (2023), eqs. 21, 22
//   and 36, or of Definition 2, eqs. 23, 25 and 26, its base node count, and
//   moments m_0, ..., m_{2k} reproduced to first order in the errors of the
//   nodes and weights (eq. 14):
//   - hermite: Gaussian moments with nu = 1 and m0 = 2, N_G = 3N and 101, give
//     the Gauss-Hermite rule of the Gaussian (eq. 20): mean 0 and standard
//     deviation 1 for N = 1, 3, 5, 8, and mean 3 and 1/2 for N = 1, 3, 5.
//     Their moments are exact in double and their recurrence coefficients
//     small dyadic numbers, so these cases mainly test the continuation and
//     the eigensolve; the trimodal cases, whose moments are rounded, test
//     rounding in the recurrence. The shifted Gaussian at N = 8 is left out:
//     sigma(8, 8) exceeds tau_8 only 470-fold, below the thousandfold margin
//     the script requires, as it prints;
//   - chebyshev: N = 1 and nu = 0, N_G = 3, 10, 101, give the Chebyshev rule of
//     the second kind (footnote 2);
//   - figures: standardized Gaussian moments with N = 5, N_G = 101 and
//     nu = 0, 2 (Fig. 1), and the trimodal case of section 4.1 with N_G = 101,
//     N = 3, 6, 9 and nu = 1 (Fig. 5), and N = 3 with nu = -1;
//   - boundary: the moments of three dyadic points at N = 3 give their
//     classical quadrature in both modes. Those of two points at N = 4 give it
//     in adaptive mode, and fixed mode is not_realizable;
//   - reduction: a point of weight 1e-8 makes the classical 3-node rule fail
//     rmin = 1e-6, and adaptive mode continues the 2-node rule to N_G = 9;
//   - laguerre: with the gamma continuation and the default bounds, the
//     moments of the gamma number density m0 beta^(alpha + 1) xi^alpha
//     exp(-beta xi) / Gamma(alpha + 1) give the Gauss rule of the generalized
//     Laguerre polynomials (eq. 24): alpha = 10, beta = 1 and m0 = 1 with
//     N = 2, 4, 6 and N_G = 101 (Fig. 2), and alpha = 1/2, beta = 2 and m0 = 2
//     with N = 1, 3, 5 and N_G = 10, 101;
//   - shifted_laguerre: the gamma number density with alpha = beta = m0 = 1,
//     shifted to [lower, inf) for lower = -1 and 4, with the bounds
//     [lower, the largest double], N = 1, 3, 5 and N_G = 101;
//   - gamma_mixture: 1/3 of the gamma number density with alpha = 2 and
//     beta = 1 and 2/3 of that with alpha = 9 and beta = 1/2, moments rounded
//     to double, with the default bounds, N = 3, 4, 5 and N_G = 20, 101. For a
//     gamma density, every pair of known zeta gives the same alpha, so the
//     laguerre cases cannot tell which pair the continuation reads. For this
//     mixture, the alpha of Definition 2, from m_0, m_1 and m_2, differs from
//     that of the last known pair, and its rounded moments test rounding
//     along the zeta path;
//   - gamma_boundary: the boundary cases with dyadic points in (0, inf);
//   - gamma_reduction: the moments of the points -1, 1 and 2 at N = 3 with the
//     default bounds. The classical 3-node rule has the node -1, and the
//     continuation from the 2-node rule has negative nodes, so adaptive mode
//     continues the 1-node rule to N_G = 9, and fixed mode is not_realizable.
//   Every case other than the boundary and reduction ones of either
//   continuation lies in the interior of moment space, and fixed mode returns
//   the identical result.
//
// How the reference is constructed
//   The moments are rounded to double where they are not exact in double, and
//   the references are computed from those doubles, the input the function
//   receives: the recurrence of docs/developer/qmom.md and the continuation in
//   exact rational arithmetic, then the eigenvalues and first eigenvector
//   components of the Jacobi matrix with mpmath at 50 digits, rounded to
//   double. The script generalized_quadrature_exact.py, beside this file,
//   generated generalized_quadrature_exact.inc. It asserts that the closed
//   forms of the hermite, chebyshev, laguerre and shifted_laguerre cases agree
//   with the references to 1e-40, and that the decisions the test expects are
//   not marginal: sigma(k, k) exceeds its rounding bound tau_k a thousandfold
//   where the result is continued, the rmin and eabs tests on the classical
//   rules pass or fail with a factor of ten to spare, and the nodes the test
//   expects inside or outside the bounds lie there by at least a thousandfold
//   the margin of the bounds check. The script computes the gamma
//   continuation as Definition 2 states it, from the moments of the
//   coordinate minus the lower bound.
//
// Error metric and pass criteria
//   Positions relative to the largest |position| of the reference, weights
//   relative to m0, both to one tolerance per case. This test sets each
//   tolerance to the smallest power of ten at least 30 times the largest
//   error observed in the Debug, Release, GCC 12 and GCC 16 builds. The
//   errors depend on the build, for example on whether it fuses
//   multiply-adds, and grow with the conditioning of the case:
//   - 1e-4 for the trimodal case at N = 9, 1.4e-7 in positions and 3.9e-7 in
//     weights. A change of one ulp in its moments moves the nodes by 9.1e-7,
//     which the script prints, so these errors are what the conditioning of
//     the case allows in double precision;
//   - 1e-8 for the trimodal case at N = 6, 6.1e-11;
//   - 1e-10 for the gamma mixture at N = 5, 2.4e-12, and 1e-11 at N = 4,
//     2.2e-13. A change of one ulp in their moments moves the nodes by
//     7.4e-12 and 1.6e-12, which the script prints;
//   - 1e-11 for reduction (rmin_n3), 2.7e-13;
//   - 1e-12 for every other case, whose largest error is 1.2e-14, except
//     3.0e-14 for the trimodal case at N = 3.
//   The moments of the result must lie within the first-order bound
//   sum_i (dw |xi_i|^k + w_i k |xi_i|^(k-1) dxi) of the tolerance, plus
//   gamma_{k+N} sum_i w_i |xi_i|^k for evaluating the sum in double. Weights
//   relative to m0 check the small weights of the tails only in absolute
//   terms, to the tolerance times m0, not to their own relative accuracy.
//
// Reproduce
//   ctest -R verification_generalized_quadrature_exact
//   The references: python3 generalized_quadrature_exact.py, with mpmath,
//   writes generalized_quadrature_exact.inc.
//
// Build modes
//   Both. QMOMKIT_VERIFY_CHECK and QMOMKIT_VERIFY_NEAR report and exit rather
//   than relying on assert(), so NDEBUG does not weaken this study.

#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace
{

struct GeneralizedReference
{
    const char* name;
    std::size_t n;
    std::size_t continued_node_count;
    // Read by the normal continuation only.
    double nu;
    double rmin;
    std::size_t base_node_count;
    std::size_t reproduced_order;
    bool check_fixed;
    bool fixed_refused;
    std::vector<double> moments;
    std::vector<double> abscissas;
    std::vector<double> weights;
    qmomkit::GeneralizedContinuation continuation = qmomkit::GeneralizedContinuation::normal;
    std::optional<qmomkit::AbscissaBounds> bounds{};
};

#include "generalized_quadrature_exact.inc"

// The tolerance of a case, relative to the largest |position| and to m0,
// from the observed errors in the header.
double tolerance_for(const GeneralizedReference& reference)
{
    const std::string name = reference.name;
    if (name == "trimodal_n9_nu1")
    {
        return 1.0e-4;
    }
    if (name == "trimodal_n6_nu1")
    {
        return 1.0e-8;
    }
    if (name == "mixture_n5_N20" || name == "mixture_n5_N101")
    {
        return 1.0e-10;
    }
    if (name == "rmin_n3" || name == "mixture_n4_N20" || name == "mixture_n4_N101")
    {
        return 1.0e-11;
    }
    return 1.0e-12;
}

double power(const double value, const std::size_t exponent)
{
    double result = 1.0;
    for (std::size_t i = 0; i < exponent; ++i)
    {
        result *= value;
    }
    return result;
}

void verify_moments(const GeneralizedReference& reference, const qmomkit::QmomResult& result)
{
    const std::size_t nodes = reference.weights.size();
    double scale = 0.0;
    for (const double abscissa : reference.abscissas)
    {
        scale = std::max(scale, std::abs(abscissa));
    }
    const double weight_error = tolerance_for(reference) * reference.moments[0];
    const double abscissa_error = tolerance_for(reference) * scale;
    for (std::size_t k = 0; k <= reference.reproduced_order; ++k)
    {
        double sum = 0.0;
        double absolute_sum = 0.0;
        double bound = 0.0;
        for (std::size_t i = 0; i < nodes; ++i)
        {
            sum += result.weights[i] * power(result.abscissas[i], k);
            const double magnitude = std::abs(reference.abscissas[i]);
            absolute_sum += reference.weights[i] * power(magnitude, k);
            bound += weight_error * power(magnitude, k);
            if (k > 0)
            {
                bound += reference.weights[i] * static_cast<double>(k) * power(magnitude, k - 1) * abscissa_error;
            }
        }
        const double terms = static_cast<double>(k + nodes) * std::numeric_limits<double>::epsilon();
        bound += terms / (1.0 - terms) * absolute_sum;
        QMOMKIT_VERIFY_NEAR(sum, reference.moments[k], bound);
    }
}

void verify_reference(const GeneralizedReference& reference)
{
    qmomkit::QmomOptions adaptive;
    adaptive.family = qmomkit::QmomFamily::generalized;
    adaptive.rmin = reference.rmin;
    adaptive.generalized.node_count = reference.continued_node_count;
    adaptive.generalized.nu = reference.nu;
    adaptive.generalized.continuation = reference.continuation;
    const qmomkit::MomentSet moments(reference.moments);
    const qmomkit::QmomResult result = qmomkit::compute_qmom(moments, reference.n, adaptive, reference.bounds);
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.base_node_count == reference.base_node_count);
    const std::size_t nodes = reference.weights.size();
    QMOMKIT_VERIFY_CHECK(result.abscissas.size() == nodes && result.weights.size() == nodes);

    double scale = 0.0;
    for (const double abscissa : reference.abscissas)
    {
        scale = std::max(scale, std::abs(abscissa));
    }
    for (std::size_t i = 0; i < nodes; ++i)
    {
        QMOMKIT_VERIFY_NEAR(result.abscissas[i], reference.abscissas[i], tolerance_for(reference) * scale);
        QMOMKIT_VERIFY_NEAR(result.weights[i], reference.weights[i], tolerance_for(reference) * reference.moments[0]);
    }
    verify_moments(reference, result);

    if (!reference.check_fixed)
    {
        return;
    }
    qmomkit::QmomOptions fixed = adaptive;
    fixed.adaptive = false;
    const qmomkit::QmomResult fixed_result = qmomkit::compute_qmom(moments, reference.n, fixed, reference.bounds);
    if (reference.fixed_refused)
    {
        QMOMKIT_VERIFY_CHECK(fixed_result.status == qmomkit::QmomStatus::not_realizable && fixed_result.base_node_count == 0);
        return;
    }
    QMOMKIT_VERIFY_CHECK(fixed_result.ok() && fixed_result.abscissas == result.abscissas && fixed_result.weights == result.weights &&
                         fixed_result.base_node_count == result.base_node_count);
}

void verify_references(const std::vector<GeneralizedReference>& references)
{
    for (const GeneralizedReference& reference : references)
    {
        verify_reference(reference);
    }
}

void test_hermite()
{
    verify_references(hermite_references);
}

void test_chebyshev()
{
    verify_references(chebyshev_references);
}

void test_figures()
{
    verify_references(figures_references);
}

void test_boundary()
{
    verify_references(boundary_references);
}

void test_reduction()
{
    verify_references(reduction_references);
}

void test_laguerre()
{
    verify_references(laguerre_references);
}

void test_shifted_laguerre()
{
    verify_references(shifted_laguerre_references);
}

void test_gamma_mixture()
{
    verify_references(gamma_mixture_references);
}

void test_gamma_boundary()
{
    verify_references(gamma_boundary_references);
}

void test_gamma_reduction()
{
    verify_references(gamma_reduction_references);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"hermite", test_hermite},
                                                 {"chebyshev", test_chebyshev},
                                                 {"figures", test_figures},
                                                 {"boundary", test_boundary},
                                                 {"reduction", test_reduction},
                                                 {"laguerre", test_laguerre},
                                                 {"shifted_laguerre", test_shifted_laguerre},
                                                 {"gamma_mixture", test_gamma_mixture},
                                                 {"gamma_boundary", test_gamma_boundary},
                                                 {"gamma_reduction", test_gamma_reduction},
                                             });
}
