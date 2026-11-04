// HyQMOM characteristic speeds, against the eigenvalues of the Jacobian of the
// closed moment system computed in exact rational arithmetic, and against the
// speeds of the standard normal distribution and of rules on the boundary.
//
// What is verified
//   That the 2n + 1 speeds of the moment system of free transport,
//   d_t m_k + d_x m_{k+1} = 0 for k = 0, ..., 2n, closed by HyQMOM are the
//   eigenvalues of J_n and K_{n+1} (Fox & Laurent 2022, eqs. 4.17 and 4.18),
//   ascending and distinct, and that a state below n + 1 nodes gives the
//   speeds of its classical rule of k nodes, its abscissas each twice and
//   their mean. The cases:
//   - normal: states of normal moments at node counts n = 1 to 3, whose
//     speeds are mu + sigma s for the roots s of He_n and of R_{n+1} of the
//     standard normal distribution: 0 and -+sqrt(3) for n = 1; -+1 and 0,
//     -+sqrt(6) for n = 2; 0, -+sqrt(3) and -+sqrt(5 -+ 3 sqrt(2)) for n = 3,
//     the roots of x^4 - 10 x^2 + 7. Fox & Laurent 2022, section 5.3, prints
//     x^4 - 10 x^2 + 9, a typo: eq. 4.11 with a_k = 0, b_k = k and
//     beta_3 = 7 gives + 7.
//   - boundary: states of k <= n weighted points, two at n = 2, and three,
//     two and one at n = 3, whose speeds are the points, each twice, and
//     their mean.
//   - skewed and bimodal: states of the moments of the skewed and bimodal
//     cases of hyperbolic_quadrature_exact.cpp, n = 1 to 10, which its
//     generator asserts to lie in the interior of moment space.
//
// How the reference is constructed
//   Normal and boundary: mean, standard deviation, points and weights are
//   dyadic, and the script below computes the moments in exact arithmetic and
//   asserts that double holds them exactly. The normal speeds are evaluated in
//   double from the roots above; the boundary speeds are the points and their
//   mean. Skewed and bimodal: from the double moments of each shape, the
//   script below closes m_{2n+1} in exact rational arithmetic, by the
//   recurrence of docs/developer/qmom.md, the closure a_n = (a_0 + ... +
//   a_{n-1}) / n and the reverse step of Fox & Laurent 2022, Algorithm 2, and
//   differentiates it exactly with dual numbers. The speeds are the roots of
//   lambda^{2n+1} - sum_j (d m_{2n+1} / d m_j) lambda^j, the characteristic
//   polynomial of the Jacobian, found with mpmath at 60 digits. The script
//   asserts that these roots are real, that they agree with the eigenvalues of
//   J_n and K_{n+1} to 1e-40, that the roots of R_{n+1} strictly separate
//   those of Q_n (Fox & Laurent 2022, Theorem 4.8), and that the closure
//   equals m_{2n+1} of the J_{n+1} rule. The smallest gap between two speeds
//   is 2.3e-3 of the largest |speed|, at skewed n = 10.
//
// Error metric and pass criteria
//   Speeds relative to the largest |speed| of the reference, as positions in
//   hyperbolic_quadrature_exact.cpp. Every case runs through a state, whose
//   speeds come from the nodes of the reconstruction of the given moments; the
//   map from 2n + 1 moments to the rule worsens in conditioning with n. The
//   largest observed deviation is 1.6e-16 for the boundary cases, 2.1e-16 for
//   the normal ones, 5.4e-15 up to n = 5, 3.6e-14 for n = 6 and 7, 2.7e-13
//   for n = 8 and 1.3e-11 for n = 9 and 10, alike to two digits in Debug,
//   Release and GCC builds. The tolerances are those of
//   hyperbolic_quadrature_exact.cpp: 1e-12 up to n = 7, 1e-11 for n = 8 and
//   5e-10 for n = 9 and 10, a margin of 27 or more over those deviations, and
//   1e-12 for the normal and boundary cases, a safety bound.
//
// Reproduce
//   ctest -R verification_hyperbolic_characteristic_speeds_exact
//   The references: python3 with mpmath and the script below.
//
// Build modes
//   Both. QMOMKIT_VERIFY_CHECK and QMOMKIT_VERIFY_NEAR report and exit rather
//   than relying on assert(), so NDEBUG does not weaken this study.
//
// Generator of the references:
//
//   from fractions import Fraction as F
//   from math import comb
//   import mpmath as mp
//
//   mp.mp.dps = 60
//
//
//   class D:
//       # A dual number over the rationals: a value and its exact gradient in
//       # m_0, ..., m_{2n}.
//       def __init__(self, v, g):
//           self.v, self.g = v, g
//
//       def __add__(s, o):
//           return D(s.v + o.v, [x + y for x, y in zip(s.g, o.g)])
//
//       def __sub__(s, o):
//           return D(s.v - o.v, [x - y for x, y in zip(s.g, o.g)])
//
//       def __mul__(s, o):
//           return D(s.v * o.v, [s.v * y + o.v * x for x, y in zip(s.g, o.g)])
//
//       def __truediv__(s, o):
//           return D(s.v / o.v, [(x * o.v - s.v * y) / (o.v * o.v) for x, y in zip(s.g, o.g)])
//
//
//   def shapes(n):
//       # The asymmetric shapes of hyperbolic_quadrature_exact.cpp: a skewed
//       # geometric rule at number density 2500, and two binomial clusters of
//       # masses 3/4 and 1/4, shifted off the origin.
//       p = 2 * n + 2
//       skewed = ([F(j - n - 1, 2) for j in range(p + 1)], [2500 * F(3, 4)**j / 4 for j in range(p + 1)])
//       bimodal = ([1 + c + F(j, 2) - F(n, 4) for c in (-2, 2) for j in range(n + 1)],
//                  [mass * F(comb(n, j), 2**n) for mass in (F(3, 4), F(1, 4)) for j in range(n + 1)])
//       return {"skewed": skewed, "bimodal": bimodal}
//
//
//   def closure(m, n):
//       # m_{2n+1} of the HyQMOM closure and its exact gradient: the recurrence of
//       # m_0, ..., m_{2n}, docs/developer/qmom.md, the closure a_n = mean(a_0, ...,
//       # a_{n-1}), then the reverse step of Fox & Laurent 2022, Algorithm 2.
//       size = 2 * n + 1
//       zero = D(F(0), [F(0)] * size)
//       s = {(-1, l): zero for l in range(size + 1)}
//       for l in range(size):
//           s[(0, l)] = D(m[l], [F(int(j == l)) for j in range(size)])
//       a, b = [s[(0, 1)] / s[(0, 0)]], [zero]
//       for k in range(1, n + 1):
//           for l in range(k, 2 * n - k + 1):
//               s[(k, l)] = s[(k - 1, l + 1)] - a[k - 1] * s[(k - 1, l)] - b[k - 1] * s[(k - 2, l)]
//           b.append(s[(k, k)] / s[(k - 1, k - 1)])
//           if k < n:
//               a.append(s[(k, k + 1)] / s[(k, k)] - s[(k - 1, k)] / s[(k - 1, k - 1)])
//       total = zero
//       for value in a:
//           total = total + value
//       a.append(total / D(F(n), [F(0)] * size))
//       s[(n, n + 1)] = s[(n, n)] * (a[n] + s[(n - 1, n)] / s[(n - 1, n - 1)])
//       for k in range(n - 1, -1, -1):
//           s[(k, 2 * n - k + 1)] = s[(k + 1, 2 * n - k)] + a[k] * s[(k, 2 * n - k)] + b[k] * s[(k - 1, 2 * n - k)]
//       return s[(0, size)], [v.v for v in a], [v.v for v in b]
//
//
//   def q(v):
//       return mp.mpf(v.numerator) / v.denominator
//
//
//   def eigenvalues(a, off, size):
//       J = mp.matrix(size, size)
//       for i in range(size):
//           J[i, i] = q(a[i])
//       for i in range(1, size):
//           J[i - 1, i] = J[i, i - 1] = mp.sqrt(q(off[i]))
//       E, Q = mp.eigsy(J)
//       return sorted(E[i] for i in range(size)), E, Q
//
//
//   def speeds(m, n):
//       # The eigenvalues of the Jacobian of the closed system, rows e_{k+1} for
//       # k < 2n and the gradient of m_{2n+1}: the roots of
//       # lambda^{2n+1} - sum_j (d m_{2n+1} / d m_j) lambda^j.
//       m_next, a, b = closure(m, n)
//       coefficients = [mp.mpf(1)] + [-q(c) for c in reversed(m_next.g)]
//       roots = mp.polyroots(coefficients, maxsteps=4000, extraprec=800)
//       assert all(abs(mp.im(r)) < mp.mpf(10)**-40 for r in roots)
//       roots = sorted(mp.re(r) for r in roots)
//       scale = max(abs(r) for r in roots)
//       # They are the eigenvalues of J_n and K_{n+1} (eqs. 4.17 and 4.18), which
//       # interlace (Theorem 4.8): K_{n+1} at even positions, J_n at odd ones.
//       beta = b[:n] + [F(2 * n + 1, n) * b[n]]
//       j_roots, _, _ = eigenvalues(a, b, n)
//       k_roots, _, _ = eigenvalues(a, beta, n + 1)
//       assert max(abs(x - y) for x, y in zip(roots[1::2], j_roots)) < mp.mpf(10)**-40 * scale
//       assert max(abs(x - y) for x, y in zip(roots[0::2], k_roots)) < mp.mpf(10)**-40 * scale
//       assert all(x < y for x, y in zip(roots, roots[1:]))
//       # The closure is the moment m_{2n+1} of the HyQMOM rule, the Gauss rule
//       # of J_{n+1}.
//       _, E, Q = eigenvalues(a, b, n + 1)
//       rule = sum(q(m[0]) * Q[0, i]**2 * E[i]**(2 * n + 1) for i in range(n + 1))
//       assert abs(rule - q(m_next.v)) < mp.mpf(10)**-40 * max(1, abs(rule))
//       gap = min(y - x for x, y in zip(roots, roots[1:])) / scale
//       return roots, gap
//
//
//   fmt = lambda values: ", ".join(repr(float(v)) for v in values)
//   for name in ["skewed", "bimodal"]:
//       print(f"const std::vector<SpeedReference> {name}_references{{")
//       for n in range(1, 11):
//           x, w = shapes(n)[name]
//           m = [F(float(sum(wi * xi**k for xi, wi in zip(x, w)))) for k in range(2 * n + 1)]
//           roots, gap = speeds(m, n)
//           print(f"    {{{n}, {{{fmt(m)}}}, {{{fmt(roots)}}}}},  // gap {mp.nstr(gap, 3)}")
//       print("};")
//
//
//   def exact(values):
//       # Moments that double holds exactly.
//       assert all(F(float(v)) == v for v in values)
//       return values
//
//
//   def normal_moments(m0, mu, sigma, n):
//       # m_k = m0 E[(mu + sigma Z)^k], with E[Z^j] = (j - 1)!! for even j.
//       double_factorial = lambda j: 1 if j <= 0 else j * double_factorial(j - 2)
//       return exact([m0 * sum(comb(k, j) * mu**(k - j) * sigma**j * double_factorial(j - 1) for j in range(0, k + 1, 2)) for k in range(2 * n + 1)])
//
//
//   print("const std::vector<NormalReference> normal_references{")
//   for m0, mu, sigma in [(F(2), F(3, 4), F(3, 2)), (F(1, 2), F(-5, 4), F(1, 2))]:
//       for n in range(1, 4):
//           print(f"    {{{n}, {float(mu)!r}, {float(sigma)!r}, {{{fmt(normal_moments(m0, mu, sigma, n))}}}}},")
//   print("};")
//
//   print("const std::vector<BoundaryReference> boundary_references{")
//   for n, x, w in [(2, [F(-5, 4), F(11, 4)], [F(1, 4), F(3, 4)]), (3, [F(-3, 2), F(1, 4), F(11, 4)], [F(1, 4), F(1, 2), F(1, 4)]),
//                   (3, [F(-5, 4), F(11, 4)], [F(1, 4), F(3, 4)]), (3, [F(1, 2)], [F(2)])]:
//       m = exact([sum(wi * xi**k for xi, wi in zip(x, w)) for k in range(2 * n + 1)])
//       print(f"    {{{n}, {{{fmt(m)}}}, {{{fmt(sorted(x + x + [sum(x) / len(x)]))}}}}},")
//   print("};")

#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

struct NormalReference
{
    std::size_t n;
    double mean;
    double deviation;
    std::vector<double> moments;
};

struct BoundaryReference
{
    std::size_t n;
    std::vector<double> moments;
    std::vector<double> speeds;
};

struct SpeedReference
{
    std::size_t n;
    std::vector<double> moments;
    std::vector<double> speeds;
};

const std::vector<SpeedReference> skewed_references{
    {1, {1906.73828125, -529.78515625, 1005.859375}, {-1.4401683060260546, -0.27784891165172854, 0.8844704827225974}},
    {2,
     {2166.290283203125, -1167.9840087890625, 2410.9268188476562, -2073.6122131347656, 4378.595352172852},
     {-1.7794211153500687, -1.1794289710374264, -0.43774974713851716, 0.745040746570402, 1.5655885257880493}},
    {3,
     {2312.288284301758, -2000.8468627929688, 4558.477401733398, -6054.644584655762, 14285.075664520264, -21167.715787887573, 50763.39900493622},
     {-2.2089659497558096, -1.8329521297645959, -1.3208884648106287, -0.36732578387504483, 0.5493558707623952, 1.5343108075312626,
      2.0925424023262056}},
    {4,
     {2394.4121599197388, -2975.1452803611755, 7545.082718133926, -13756.671771407127, 36502.71091610193, -73473.83892163634, 199253.2003019005,
      -419099.66898616403, 1150354.7028289177},
     {-2.6571205656212133, -2.406485170577068, -2.0108808951844437, -1.2500792800533713, -0.49749807634945425, 0.5035280765586867, 1.3762971021691042,
      2.195006424697774, 2.591664998268534}},
    {5,
     {2440.606839954853, -4046.9658002257347, 11461.991518735886, -26657.965388149023, 79605.64850829542, -202169.53322757035, 619497.4381173961,
      -1637902.138066711, 5093283.868140425, -13770086.759769183, 43210961.905376025},
     {-3.1136925688757566, -2.9459138257940114, -2.613218118493308, -1.9890939593595145, -1.3700026615453482, -0.47677372676087465,
      0.3538374479637336, 1.314281569595781, 2.092885860473, 2.791574263778195, 3.08307922622917}},
    {6,
     {2466.591347474605, -5183.747588889673, 16389.16230353061, -46312.88371543633, 155164.4275210856, -474039.91843115364, 1634010.4200513815,
      -5177983.017803363, 18144957.834031165, -58709631.144430175, 207889409.35737365, -681497625.6810303, 2430001846.586147},
     {-3.5776236310309915, -3.468310803145069, -3.1671778688838765, -2.6456810705345384, -2.124657249850629, -1.3386732867571225, -0.5987846256777711,
      0.3331896261963086, 1.1632062407151762, 2.057507622364052, 2.7397556255335394, 3.3534806744236327, 3.5720463988330273}},
    {7,
     {2481.207632954465, -6362.754202273209, 22392.097685224144, -74333.7589562725, 277976.93596203317, -989115.4775384894, 3811165.6601780397,
      -14020855.882794509, 54979243.89775543, -206230121.2693944, 817845465.3737432, -3106289585.837921, 12413355918.631016, -47549220884.27788,
      191045680946.67596},
     {-4.0495256612745845, -3.98139452964148, -3.689971143453025, -3.2521266564068787, -2.8036908001357377, -2.1079654933633987, -1.4541224003464062,
      -0.5902645099086131, 0.19998245202637974, 1.1254454030104954, 1.9204990063146425, 2.7459258616348183, 3.3415872267876447, 3.895064762797971,
      4.060595420792988}},
    {8,
     {2489.4292935368867, -7568.7095920102365, 29521.675532544123, -112382.35900838731, 466216.6479112017, -1888635.4139253213, 8077486.650411985,
      -33738317.72896783, 146921459.1302053, -624931776.2437028, 2753546379.719804, -11852465691.534496, 52648303122.93557, -228497572831.86853,
      1020884555935.2437, -4457064011188.101, 19998508763676.832},
     {-4.5296048785805425, -4.489120655413785, -4.190900408006874, -3.8262407067061783, -3.4322476333047254, -2.8108861593471337, -2.2303617275329115,
      -1.438018862287728, -0.7084137874907576, 0.1870398298052856, 0.9869479668592022, 1.8842932929899119, 2.6310077616393754, 3.392765063229097,
      3.913110352990816, 4.423933721860695, 4.549698568072851}},
    {9,
     {2494.0539776144988, -8791.622156698508, 37815.80683595393, -162164.01832356126, 741526.2405970474, -3363407.7871257607, 15854708.854400141,
      -73977018.9637701, 355105941.01243275, -1685530848.8588831, 8188030951.409323, -39310149685.55352, 192558006507.68842, -931893190046.4677,
      4592457896088.602, -22355798268955.23, 110668743551186.84, -541095464909048.2, 2687799548115842.5},
     {-5.016820175826822, -4.993685431674, -4.676722447784404, -4.378545344787504, -4.025370217771874, -3.4660838695814045, -2.9467259326797097,
      -2.2199043605621736, -1.5516494948405655, -0.7052594554371885, 0.06166414666263719, 0.9614308850526975, 1.7482107742919235, 2.608061097796523,
      3.3028558503986445, 4.0083883440468275, 4.463410463196715, 4.944442238491908, 5.0398627047375495}},
    {10,
     {2496.6553624081557,     -10025.084781938833,   47301.83022268322,    -225421.90071622323, 1129078.9270932479, -5662472.1712942775,
      29222292.60182534,      -150465233.23406377,   790600099.830338,     -4137367824.2340574, 22000210196.705147, -116397507813.75946,
      624158911120.998,       -3328095482910.079,    17956490710833.46,    -96299410310567.08,  521994013148592.2,  -2811681587349518.5,
      1.5295371439514244e+16, -8.26664632916014e+16, 4.509586135115074e+17},
     {-5.509252307091799, -5.496367709027718,  -5.153221439350712,  -4.915659914004084,  -4.5922732159350295, -4.0862495320785435,
      -3.617905753426561, -2.9488613271861492, -2.3369718110156477, -1.5431822067185081, -0.8211547738229005, 0.0517151289031325,
      0.8310461483683911, 1.719123339442227,   2.480471132733749,   3.2998088562801398,  3.9435368493223373,  4.600400175434689,
      4.9983595876088485, 5.459244182977154,   5.531333676033897}},
};
const std::vector<SpeedReference> bimodal_references{
    {1, {1.0, 0.0, 3.0625}, {-3.031088913245535, 0.0, 3.031088913245535}},
    {2, {1.0, 0.0, 3.125, 6.0, 23.28125}, {-1.2658808658067997, -1.0516162655934158, 0.6847074631082078, 2.9716162655934157, 3.461173402698592}},
    {3,
     {1.0, 0.0, 3.1875, 6.0, 24.45703125, 71.25, 245.798583984375},
     {-2.5295742873647553, -1.2767419172736196, -1.006264146579963, 0.11360992199006094, 2.555585838850039, 3.19068020464995, 3.683650207583201}},
    {4,
     {1.0, 0.0, 3.25, 6.0, 25.65625, 75.0, 268.9140625, 893.8125, 3162.970703125},
     {-1.9212451064758034, -1.6495402017466672, -1.010204580363288, -0.6848151789715347, 0.49090589338277313, 2.532325146448094, 3.0493687168124834,
      3.5768533075965108, 4.109703918301839}},
    {5,
     {1.0, 0.0, 3.3125, 6.0, 26.87890625, 78.75, 293.157470703125, 993.0703125, 3640.607925415039, 13130.302734375, 48949.67107105255},
     {-2.6842245757298424, -1.930515529518236, -1.6153097614496106, -0.9783514553808118, -0.6053652391493506, 0.26542135964086744, 2.198340757445132,
      2.847678554224984, 3.3322504772786123, 3.8669400631789226, 4.2597159321799305}},
    {6,
     {1.0, 0.0, 3.375, 6.0, 28.125, 82.5, 318.55078125, 1097.25, 4161.005859375, 15403.59375, 59430.741943359375, 229576.48828125, 906680.3247070312},
     {-2.5403358710605217, -2.246492101146379, -1.7501889039687621, -1.3652880045997124, -0.8119682322128231, -0.39201495636528155,
      0.4549525240171007, 2.156563889454139, 2.743149025025261, 3.249806409733509, 3.762150666941274, 4.191410925215742, 4.668557980599157}},
    {7,
     {1.0, 0.0, 3.4375, 6.0, 29.39453125, 86.25, 345.115966796875, 1206.3515625, 4726.183242797852, 17913.134765625, 71396.32638454437,
      283950.08615112305, 1159949.9033668637, 4771716.448802948, 19964867.77930809},
     {-3.012321430049305, -2.5356643197065454, -2.1997777562360183, -1.689041592940346, -1.3008710067354003, -0.7537476188291211,
      -0.31666952746409227, 0.41472792320003193, 1.912584404529286, 2.5906932570413677, 3.0727899672346504, 3.5945754965049224, 4.016627926007771,
      4.490167841020274, 4.812449978473776}},
    {8,
     {1.0, 0.0, 3.5, 6.0, 30.6875, 90.0, 372.875, 1320.375, 5338.197265625, 20670.0, 84973.3291015625, 347091.29296875, 1463676.3436279297,
      6206511.03515625, 26844734.510131836, 117140157.34643555, 518339002.90265656},
     {-3.104165257450218, -2.8276403077246584, -2.4159481519934674, -2.0249441519730267, -1.5762736694489283, -1.1445521000298373,
      -0.6481325028147144, -0.17293860548300768, 0.5123026187197909, 1.8590484947479764, 2.504456684685326, 2.9991386580821677, 3.504876857328948,
      3.9400383514696644, 4.401331001333805, 4.788280413868706, 5.16503701671719}},
    {9,
     {1.0, 0.0, 3.5625, 6.0, 32.00390625, 93.75, 401.849853515625, 1439.3203125, 5999.143569946289, 23685.263671875, 100294.4586725235,
      419856.6303405762, 1824597.5784035325, 7956866.686420441, 35505583.545894235, 159746474.70682776, 730247516.4005091, 3368991707.626599,
      15730191242.099619},
     {-3.4242240077698645, -3.110566374251313, -2.769772597882723, -2.343475223646785, -1.9607740538296614, -1.4996243161645428, -1.0811883944477407,
      -0.581580708954464, -0.112374663476971, 0.5364983208974893, 1.6980354261533492, 2.3934195854843847, 2.875477498059148, 3.383708041246774,
      3.813262934713826, 4.272778828483964, 4.6550033573187815, 5.078530955930239, 5.339542400079352}},
    {10,
     {1.0,
      0.0,
      3.625,
      6.0,
      33.34375,
      97.5,
      432.0625,
      1563.1875,
      6711.15625,
      26970.0,
      117498.33569335938,
      503140.6875,
      2250003.0205688477,
      10070424.799804688,
      46279566.01956177,
      214231184.19396973,
      1009708406.0938034,
      4802088140.427246,
      23145922502.42427,
      112567773713.92671,
      553265153998.2307},
     {-3.615429020090759,  -3.3889549730164634, -3.0239037292616753, -2.65592591887067,     -2.2604372630430443, -1.8470660699953396,
      -1.4214935122603811, -0.9764179684333214, -0.5109935105147398, -0.017071430943982953, 0.5891359046598842,  1.6476866383631135,
      2.318118778664144,   2.8112620158468244,  3.303815987736173,   3.7427439846092754,    4.1890523305415375,  4.592873617927588,
      4.994731812154773,   5.361743761366499,   5.6353632439529635}},
};
const std::vector<NormalReference> normal_references{
    {1, 0.75, 1.5, {2.0, 1.5, 5.625}},
    {2, 0.75, 1.5, {2.0, 1.5, 5.625, 10.96875, 46.1953125}},
    {3, 0.75, 1.5, {2.0, 1.5, 5.625, 10.96875, 46.1953125, 133.365234375, 619.72119140625}},
    {1, -1.25, 0.5, {0.5, -0.625, 0.90625}},
    {2, -1.25, 0.5, {0.5, -0.625, 0.90625, -1.4453125, 2.486328125}},
    {3, -1.25, 0.5, {0.5, -0.625, 0.90625, -1.4453125, 2.486328125, -4.55322265625, 8.7994384765625}},
};
const std::vector<BoundaryReference> boundary_references{
    {2, {1.0, 1.75, 6.0625, 15.109375, 43.50390625}, {-1.25, -1.25, 0.75, 2.75, 2.75}},
    {3, {1.0, 0.4375, 2.484375, 4.36328125, 15.5654296875, 37.421142578125, 110.97528076171875}, {-1.5, -1.5, 0.25, 0.25, 0.5, 2.75, 2.75}},
    {3, {1.0, 1.75, 6.0625, 15.109375, 43.50390625, 117.1943359375, 325.336181640625}, {-1.25, -1.25, 0.75, 2.75, 2.75}},
    {3, {2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125}, {0.5, 0.5, 0.5}},
};

// The speeds of the standard normal distribution, ascending.
std::vector<double> standard_normal_speeds(const std::size_t n)
{
    switch (n)
    {
    case 1:
        return {-std::sqrt(3.0), 0.0, std::sqrt(3.0)};
    case 2:
        return {-std::sqrt(6.0), -1.0, 0.0, 1.0, std::sqrt(6.0)};
    default:
        return {-std::sqrt(5.0 + (3.0 * std::sqrt(2.0))), -std::sqrt(3.0), -std::sqrt(5.0 - (3.0 * std::sqrt(2.0))), 0.0,
                std::sqrt(5.0 - (3.0 * std::sqrt(2.0))),  std::sqrt(3.0),  std::sqrt(5.0 + (3.0 * std::sqrt(2.0)))};
    }
}

double tolerance_for(const std::size_t n)
{
    if (n <= 7)
    {
        return 1.0e-12;
    }
    return n == 8 ? 1.0e-11 : 5.0e-10;
}

void verify_speeds(const std::vector<double>& actual, const std::vector<double>& expected, const double tolerance)
{
    QMOMKIT_VERIFY_CHECK(actual.size() == expected.size());
    const double scale = std::max(std::abs(expected.front()), std::abs(expected.back()));
    for (std::size_t i = 0; i < expected.size(); ++i)
    {
        QMOMKIT_VERIFY_NEAR(actual[i], expected[i], tolerance * scale);
    }
}

void verify_strictly_ascending(const std::vector<double>& speeds)
{
    for (std::size_t i = 1; i < speeds.size(); ++i)
    {
        QMOMKIT_VERIFY_CHECK(speeds[i - 1] < speeds[i]);
    }
}

qmomkit::QmomState hyperbolic_state(const std::vector<double>& moments, const std::size_t n)
{
    qmomkit::QmomOptions options;
    options.family = qmomkit::QmomFamily::hyperbolic;
    return qmomkit::QmomCase{qmomkit::MomentSet(moments), n, options}.initial_state();
}

void test_normal()
{
    for (const NormalReference& reference : normal_references)
    {
        const qmomkit::QmomState state = hyperbolic_state(reference.moments, reference.n);
        QMOMKIT_VERIFY_CHECK(state.active_node_count() == reference.n + 1);
        const qmomkit::HyqmomCharacteristicSpeeds result = qmomkit::compute_hyqmom_characteristic_speeds(state);
        QMOMKIT_VERIFY_CHECK(result.ok());
        std::vector<double> expected = standard_normal_speeds(reference.n);
        for (double& speed : expected)
        {
            speed = reference.mean + (reference.deviation * speed);
        }
        verify_speeds(result.speeds, expected, 1.0e-12);
        verify_strictly_ascending(result.speeds);
    }
}

void test_boundary()
{
    for (const BoundaryReference& reference : boundary_references)
    {
        const qmomkit::QmomState state = hyperbolic_state(reference.moments, reference.n);
        QMOMKIT_VERIFY_CHECK((2 * state.active_node_count()) + 1 == reference.speeds.size());
        const qmomkit::HyqmomCharacteristicSpeeds result = qmomkit::compute_hyqmom_characteristic_speeds(state);
        QMOMKIT_VERIFY_CHECK(result.ok());
        verify_speeds(result.speeds, reference.speeds, 1.0e-12);
    }
}

void verify_against_references(const std::vector<SpeedReference>& references)
{
    for (const SpeedReference& reference : references)
    {
        const qmomkit::QmomState state = hyperbolic_state(reference.moments, reference.n);
        QMOMKIT_VERIFY_CHECK(state.active_node_count() == reference.n + 1);
        const qmomkit::HyqmomCharacteristicSpeeds result = qmomkit::compute_hyqmom_characteristic_speeds(state);
        QMOMKIT_VERIFY_CHECK(result.ok());
        verify_speeds(result.speeds, reference.speeds, tolerance_for(reference.n));
        verify_strictly_ascending(result.speeds);
    }
}

void test_skewed()
{
    verify_against_references(skewed_references);
}

void test_bimodal()
{
    verify_against_references(bimodal_references);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(argc, argv,
                                             {
                                                 {"normal", test_normal},
                                                 {"boundary", test_boundary},
                                                 {"skewed", test_skewed},
                                                 {"bimodal", test_bimodal},
                                             });
}
