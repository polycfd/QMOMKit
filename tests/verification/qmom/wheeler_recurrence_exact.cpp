// The Wheeler recurrence that compute_wheeler_recurrence reports, against the
// same recurrence in exact rational arithmetic.
//
// What is verified
//   That a_k, b_k and sigma(k, k) agree with exact references for rules of two
//   to five weighted points, a symmetric rule with a_k = 0, and a rule at
//   number density 1024; that entry 0 holds a_0 = m1/m0, b_0 = 0,
//   sigma(0, 0) = m0 and tolerance 0; and that every rounding bound tau_k is
//   finite and positive and below sigma(k, k), so that the n points of each
//   rule are all supported.
//
// How the reference is constructed
//   Each rule has nodes and weights that are short binary fractions, so its
//   moments m_0, ..., m_{2n-1} are exact in double. Those moments, and a_k,
//   b_k and sigma(k, k) computed from them in exact rational arithmetic with
//   the recurrence of docs/developer/qmom.md, were generated offline with the
//   Python script below and stored here, the moments exactly and the
//   coefficients to 17 significant digits. For the rules with two points and
//   the symmetric rule, every coefficient is exact in double; for the others
//   they are not, so the comparison exercises the rounding of the recurrence.
//
// Error metric and pass criteria
//   b_k and sigma(k, k) relative to the reference; a_k relative to the largest
//   |node| of the rule, since a_k lies between the smallest and largest node
//   and is zero for the symmetric rule. The tolerance, 1e-13 or about 450
//   epsilon, has a margin of about 35 over the largest measured error, 13
//   epsilon (2.9e-15) for the five-point rule; the rules with exact
//   coefficients are reproduced exactly.
//
// Reproduce
//   ctest -R verification_wheeler_recurrence_exact
//   The references: python3 with the script below.
//
// Build modes
//   Both. QMOMKIT_VERIFY_CHECK and QMOMKIT_VERIFY_NEAR report and exit rather
//   than relying on assert(), so NDEBUG does not weaken this study.
//
// Generator of the references:
//
//   from fractions import Fraction as F
//
//   rules = {
//       "two_points": ([F(-1, 2), F(5, 4)], [F(1, 4), F(3, 4)]),
//       "three_points": ([F(-5, 4), F(1, 2), F(2)], [F(1, 4), F(3, 8), F(3, 8)]),
//       "four_points": ([F(-3, 2), F(-1, 4), F(3, 4), F(5, 2)], [F(1, 8), F(3, 8), F(1, 4), F(1, 4)]),
//       "five_points": ([F(-2), F(-3, 4), F(1, 4), F(3, 2), F(13, 4)],
//                       [F(1, 16), F(1, 4), F(3, 8), F(3, 16), F(1, 8)]),
//       "symmetric_three_points": ([F(-1), F(0), F(1)], [F(1, 4), F(1, 2), F(1, 4)]),
//       "three_points_at_number_density_1024": ([F(-5, 4), F(1, 2), F(2)], [F(256), F(384), F(384)]),
//   }
//
//   def recurrence(m, n):
//       sigma = {(-1, l): F(0) for l in range(2 * n + 1)}
//       for l in range(2 * n):
//           sigma[(0, l)] = m[l]
//       a, b = [m[1] / m[0]], [F(0)]
//       for k in range(1, n):
//           for l in range(k, 2 * n - k):
//               sigma[(k, l)] = sigma[(k - 1, l + 1)] - a[k - 1] * sigma[(k - 1, l)] - b[k - 1] * sigma[(k - 2, l)]
//           a.append(sigma[(k, k + 1)] / sigma[(k, k)] - sigma[(k - 1, k)] / sigma[(k - 1, k - 1)])
//           b.append(sigma[(k, k)] / sigma[(k - 1, k - 1)])
//       return a, b, [sigma[(k, k)] for k in range(n)]
//
//   for name, (x, w) in rules.items():
//       n = len(x)
//       m = [sum(wi * xi**k for xi, wi in zip(x, w)) for k in range(2 * n)]
//       assert all(F(float(v)) == v for v in m), name
//       a, b, s = recurrence(m, n)
//       print(name, [repr(float(v)) for v in m])
//       for values in (a, b, s):
//           print([f"{float(v):.17g}" for v in values])
//       print(repr(float(max(abs(v) for v in x))))

#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

struct ExactRecurrence
{
    std::vector<double> moments;
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> sigma;
    double node_scale;
};

void verify_against_exact_recurrence(const ExactRecurrence& reference)
{
    constexpr double tolerance = 1.0e-13;
    const std::size_t n = reference.a.size();
    const qmomkit::MomentSet moments(reference.moments);
    const qmomkit::WheelerRecurrence recurrence = qmomkit::compute_wheeler_recurrence(moments, n);

    QMOMKIT_VERIFY_CHECK(recurrence.ok());
    QMOMKIT_VERIFY_CHECK(recurrence.a.size() == n && recurrence.b.size() == n && recurrence.sigma.size() == n && recurrence.tolerance.size() == n);
    QMOMKIT_VERIFY_CHECK(recurrence.a[0] == moments[1] / moments[0]);
    QMOMKIT_VERIFY_CHECK(recurrence.b[0] == 0.0);
    QMOMKIT_VERIFY_CHECK(recurrence.sigma[0] == moments[0]);
    QMOMKIT_VERIFY_CHECK(recurrence.tolerance[0] == 0.0);
    for (std::size_t k = 0; k < n; ++k)
    {
        QMOMKIT_VERIFY_NEAR(recurrence.a[k], reference.a[k], tolerance * reference.node_scale);
    }
    for (std::size_t k = 1; k < n; ++k)
    {
        QMOMKIT_VERIFY_NEAR(recurrence.b[k], reference.b[k], tolerance * std::abs(reference.b[k]));
        QMOMKIT_VERIFY_NEAR(recurrence.sigma[k], reference.sigma[k], tolerance * std::abs(reference.sigma[k]));
        QMOMKIT_VERIFY_CHECK(std::isfinite(recurrence.tolerance[k]) && recurrence.tolerance[k] > 0.0);
        QMOMKIT_VERIFY_CHECK(recurrence.sigma[k] > recurrence.tolerance[k]);
    }
}

void test_two_points()
{
    verify_against_exact_recurrence({{1.0, 0.8125, 1.234375, 1.43359375}, {0.8125, -0.0625}, {0, 0.57421875}, {1, 0.57421875}, 1.25});
}

void test_three_points()
{
    verify_against_exact_recurrence({{1.0, 0.625, 1.984375, 2.55859375, 6.6337890625, 11.248779296875},
                                     {0.625, 0.20220588235294118, 0.42279411764705882},
                                     {0, 1.59375, 1.0073664576124568},
                                     {1, 1.59375, 1.605490291819853},
                                     2.0});
}

void test_four_points()
{
    verify_against_exact_recurrence({{1.0, 0.53125, 2.0078125, 3.583984375, 10.47900390625, 23.5238037109375, 62.503570556640625, 150.4854965209961},
                                     {0.53125, 0.92757852292020371, -0.083007687316728265, 0.12417916439652453},
                                     {0, 1.7255859375, 1.6083431098146264, 0.6461013428320147},
                                     {1, 1.7255859375, 2.7753342529711373, 1.7931471876523384},
                                     2.5});
}

void test_five_points()
{
    verify_against_exact_recurrence({{1.0, 0.46875, 2.15625, 4.32421875, 15.9755859375, 44.688720703125, 153.48284912109375, 473.90345764160156,
                                      1576.7133293151855, 5031.808705329895},
                                     {0.46875, 1.2422938729198185, 0.32159139483413546, 0.18772682444997923, 0.029637907796066839},
                                     {0, 1.9365234375, 2.9210429345350764, 1.9343992139594544, 0.97550370369442529},
                                     {1, 1.9365234375, 5.6566681046709535, 10.942254335305009, 10.674209630856417},
                                     3.25});
}

void test_symmetric_three_points()
{
    verify_against_exact_recurrence({{1.0, 0.0, 0.5, 0.0, 0.5, 0.0}, {0, 0, 0}, {0, 0.5, 0.5}, {1, 0.5, 0.25}, 1.0});
}

void test_three_points_at_number_density_1024()
{
    verify_against_exact_recurrence({{1024.0, 640.0, 2032.0, 2620.0, 6793.0, 11518.75},
                                     {0.625, 0.20220588235294118, 0.42279411764705882},
                                     {0, 1.59375, 1.0073664576124568},
                                     {1024, 1632, 1644.0220588235295},
                                     2.0});
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"two_points", test_two_points},
                                                 {"three_points", test_three_points},
                                                 {"four_points", test_four_points},
                                                 {"five_points", test_five_points},
                                                 {"symmetric_three_points", test_symmetric_three_points},
                                                 {"three_points_at_number_density_1024", test_three_points_at_number_density_1024},
                                             });
}
