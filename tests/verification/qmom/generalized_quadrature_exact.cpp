// GQMOM reconstructions with the normal continuation, against Gauss rules of
// the continued recurrence computed in exact rational arithmetic and
// diagonalized at 50 digits.
//
// What is verified
//   That compute_qmom(...) with the generalized family returns the N_G nodes and
//   weights of Definition 1 of Fox, Laurent & Passalacqua (2023), eqs. 21, 22
//   and 36, its base node count, and moments m_0, ..., m_{2k} reproduced to
//   first order in the errors of the nodes and weights (eq. 14):
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
//     rmin = 1e-6, and adaptive mode continues the 2-node rule to N_G = 9.
//   Every case other than the boundary and reduction ones lies in the interior
//   of moment space, and fixed mode returns the identical result.
//
// How the reference is constructed
//   The moments are rounded to double where they are not exact in double, and
//   the references are computed from those doubles, the input the function
//   receives: the recurrence of docs/developer/qmom.md and the continuation in
//   exact rational arithmetic, then the eigenvalues and first eigenvector
//   components of the Jacobi matrix with mpmath at 50 digits, rounded to
//   double. The script generalized_quadrature_exact.py, beside this file,
//   generated generalized_quadrature_exact.inc. It asserts that the closed
//   forms of the hermite and chebyshev cases agree with the references to
//   1e-40, and that the decisions the test expects are not marginal:
//   sigma(k, k) exceeds its rounding bound tau_k a thousandfold where the
//   result is continued, and the rmin and eabs tests on the classical rules
//   pass or fail with a factor of ten to spare.
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
//   - 1e-11 for the reduction, 2.7e-13;
//   - 1e-12 for every other case, whose largest error is 3.7e-15, except
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
#include <string>
#include <vector>

namespace
{

struct GeneralizedReference
{
    const char* name;
    std::size_t n;
    std::size_t continued_node_count;
    double nu;
    double rmin;
    std::size_t base_node_count;
    std::size_t reproduced_order;
    bool check_fixed;
    bool fixed_refused;
    std::vector<double> moments;
    std::vector<double> abscissas;
    std::vector<double> weights;
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
    if (name == "rmin_n3")
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
    const qmomkit::MomentSet moments(reference.moments);
    const qmomkit::QmomResult result = qmomkit::compute_qmom(moments, reference.n, adaptive);
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
    const qmomkit::QmomResult fixed_result = qmomkit::compute_qmom(moments, reference.n, fixed);
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
                                             });
}
