// HyQMOM characteristic speeds, their weights and the closure moment, against
// the eigenvalues of the Jacobian of the closed moment system computed in exact
// rational arithmetic, and against the speeds of the standard normal
// distribution and of quadratures on the boundary.
//
// What is verified
//   That the 2n + 1 speeds of the moment system of free transport,
//   d_t m_k + d_x m_{k+1} = 0 for k = 0, ..., 2n, closed by HyQMOM are the
//   eigenvalues of J_n and K_{n+1} (Fox & Laurent 2022, eqs. 4.17 and 4.18),
//   ascending and distinct, and that moments that the reconstruction takes to
//   a classical quadrature of k <= n nodes give its abscissas each twice and
//   their mean. The speeds are those of the internal function
//   compute_hyperbolic_characteristic_speeds(...) for the given moments, the
//   moments the references are computed from; a state applies it to its
//   stored moments, bit for bit (tests/unit/public_api.cpp). The cases
//   weights and closure_moment verify, for all moments below, that weights[i]
//   is the weight of speeds[i] and that closure_moment is m_{2n+1}. For n + 1
//   nodes, the weights are the unique weights with which the 2n + 1 speeds
//   have the moments m_0, ..., m_{2n}, and the closure moment is the HyQMOM
//   closure. For a classical quadrature of k nodes, the two copies of each
//   abscissa share its weight in the ratio k + 1 to k, with the share k first
//   below the mean, the mean has weight zero, and the closure moment is
//   m_{2n+1} of the quadrature. The moments:
//   - normal: normal moments at node counts n = 1 to 3, whose speeds are
//     mu + sigma s for the roots s of He_n and of R_{n+1} of the standard
//     normal distribution: 0 and -+sqrt(3) for n = 1; -+1 and 0,
//     -+sqrt(6) for n = 2; 0, -+sqrt(3) and -+sqrt(5 -+ 3 sqrt(2)) for n = 3,
//     the roots of x^4 - 10 x^2 + 7. Fox & Laurent 2022, section 5.3, prints
//     x^4 - 10 x^2 + 9, a typo: eq. 4.11 with a_k = 0, b_k = k and
//     beta_3 = 7 gives + 7.
//   - boundary: the moments of k <= n weighted points, two at n = 2, and
//     three, two and one at n = 3, whose speeds are the points, each twice,
//     and their mean.
//   - skewed and bimodal: the moments of the skewed and bimodal cases of
//     hyperbolic_quadrature_exact.cpp, n = 1 to 10, which
//     hyperbolic_quadrature_exact.py asserts to lie in the interior of moment
//     space.
//
// How the reference is constructed
//   Normal and boundary: mean, standard deviation, points and weights are
//   dyadic, and the script hyperbolic_characteristic_speeds_exact.py, beside
//   this file, computes the moments in exact arithmetic and asserts that double
//   holds them exactly. The normal speeds are evaluated in double from the
//   roots above; the boundary speeds are the points and their mean. Skewed and
//   bimodal: from the double moments of each shape, the script closes m_{2n+1}
//   in exact rational arithmetic, by the recurrence of docs/developer/qmom.md,
//   the closure a_n = (a_0 + ... + a_{n-1}) / n and the reverse step of Fox &
//   Laurent 2022, Algorithm 2, and differentiates it exactly with dual numbers.
//   The speeds are the roots of lambda^{2n+1} - sum_j (d m_{2n+1} / d m_j)
//   lambda^j, the characteristic polynomial of the Jacobian, found with mpmath
//   at 60 digits. The script asserts that these roots are real, that they agree
//   with the eigenvalues of J_n and K_{n+1} to 1e-40, that the roots of R_{n+1}
//   strictly separate those of Q_n (Fox & Laurent 2022, Theorem 4.8), and that
//   the closure equals m_{2n+1} of the J_{n+1} quadrature. The smallest gap
//   between two speeds is 2.3e-3 of the largest |speed|, at skewed n = 10.
//   Weights: for the normal, skewed and bimodal cases, the script solves the
//   Vandermonde system of the 60-digit speeds for m_0, ..., m_{2n} at 60
//   digits. It asserts that these weights are positive, that the speeds with
//   them have the closure as m_{2n+1} to 1e-30 of the sum of w |speed|^{2n+1},
//   and that they are (n + 1) / (2n + 1) times the Gauss weights of J_n and
//   n / (2n + 1) times those of K_{n+1} to 1e-30 of m_0. For the boundary
//   cases, the weights follow in exact arithmetic from the points, their
//   weights and the sharing above, a convention of QMOMKit rather than an
//   independent reference. Closure moment: for the skewed and bimodal cases, the exact
//   closure above. For the normal cases, the normal m_{2n+1}, which the script
//   asserts to equal the exact closure, since the moments are symmetric about
//   their mean. For the boundary cases, m_{2n+1} of the points.
//
// Error metric and pass criteria
//   Speeds relative to the largest |speed| of the reference, as positions in
//   hyperbolic_quadrature_exact.cpp. The map from 2n + 1 moments to the
//   speeds worsens in conditioning with n, and the tolerances are those of
//   hyperbolic_quadrature_exact.cpp: 1e-12 up to n = 7, 1e-11 for n = 8 and
//   5e-10 for n = 9 and 10. For the normal and boundary cases, 1e-12 is a
//   safety bound. The deviations depend on the build, for example on whether
//   it fuses multiply-adds. The largest errors observed in the Debug, Release,
//   GCC 12 and GCC 16 builds are 2.5e-13 up to n = 7, 1.1e-12 for n = 8 and
//   6.7e-11 for n = 9 and 10, all for the skewed case, and 2.0e-16 for the
//   normal and boundary cases. Weights relative to m_0, their sum, with the
//   tolerances of the speeds. The largest errors observed in the same builds
//   are 2.6e-13 up to n = 7, 1.0e-12 for n = 8 and 2.8e-11 for n = 9 and 10.
//   Closure moment relative to the sum of w |speed|^{2n+1} of the reference,
//   with 1e-13 for every case, since its errors did not grow with n in these
//   builds: the largest observed is 1.1e-15.
//
// Reproduce
//   ctest -R verification_hyperbolic_characteristic_speeds_exact
//   The references: python3 hyperbolic_characteristic_speeds_exact.py, with
//   mpmath, prints the arrays stored below, with the smallest gap between two
//   speeds of each skewed and bimodal case as a comment.
//
// Build modes
//   Both. QMOMKIT_VERIFY_CHECK and QMOMKIT_VERIFY_NEAR report and exit rather
//   than relying on assert(), so NDEBUG does not weaken this study.

#include "../verification_helpers.hpp"
#include "internal/hyperbolic_quadrature.hpp"
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
    std::vector<double> weights;
    double closure_moment;
};

struct BoundaryReference
{
    std::size_t n;
    std::vector<double> moments;
    std::vector<double> speeds;
    std::vector<double> weights;
    double closure_moment;
};

struct SpeedReference
{
    std::size_t n;
    std::vector<double> moments;
    std::vector<double> speeds;
    std::vector<double> weights;
    double closure_moment;
};

const std::vector<SpeedReference> skewed_references{
    {1,
     {1906.73828125, -529.78515625, 1005.859375},
     {-1.4401683060260546, -0.27784891165172854, 0.8844704827225974},
     {317.7897135416667, 1271.1588541666667, 317.7897135416667},
     -756.6319509695087},
    {2,
     {2166.290283203125, -1167.9840087890625, 2410.9268188476562, -2073.6122131347656, 4378.595352172852},
     {-1.7794211153500687, -1.1794289710374264, -0.43774974713851716, 0.745040746570402, 1.5655885257880493},
     {199.9675544602513, 867.3428882423719, 576.4917018271277, 432.431281679503, 90.05685699387094},
     -4609.840459923223},
    {3,
     {2312.288284301758, -2000.8468627929688, 4558.477401733398, -6054.644584655762, 14285.075664520264, -21167.715787887573, 50763.39900493622},
     {-2.2089659497558096, -1.8329521297645959, -1.3208884648106287, -0.36732578387504483, 0.5493558707623952, 1.5343108075312626,
      2.0925424023262056},
     {157.94188301881198, 647.2808964472697, 545.7885874216718, 521.168247022134, 251.92999728905625, 152.85844756017227, 35.320225542641964},
     -80088.61790710024},
    {4,
     {2394.4121599197388, -2975.1452803611755, 7545.082718133926, -13756.671771407127, 36502.71091610193, -73473.83892163634, 199253.2003019005,
      -419099.66898616403, 1150354.7028289177},
     {-2.6571205656212133, -2.406485170577068, -2.0108808951844437, -1.2500792800533713, -0.49749807634945425, 0.5035280765586867, 1.3762971021691042,
      2.195006424697774, 2.591664998268534},
     {145.6135252299034, 524.8581215949076, 479.4665410708299, 512.1687744716234, 317.6252396301015, 232.4144931460305, 105.34752981545313,
      60.78758852062679, 16.13034644026262},
     -2484796.316623015},
    {5,
     {2440.606839954853, -4046.9658002257347, 11461.991518735886, -26657.965388149023, 79605.64850829542, -202169.53322757035, 619497.4381173961,
      -1637902.138066711, 5093283.868140425, -13770086.759769183, 43210961.905376025},
     {-3.1136925688757566, -2.9459138257940114, -2.613218118493308, -1.9890939593595145, -1.3700026615453482, -0.47677372676087465,
      0.3538374479637336, 1.314281569595781, 2.092885860473, 2.791574263778195, 3.08307922622917},
     {150.95489500027392, 451.85326249322185, 411.43310688478687, 477.66714138294896, 341.7943555500081, 275.30256632202565, 151.5019822259773,
      100.07228298349995, 45.640916713457734, 26.34484133913258, 8.041489059520222},
     -118493057.91330819},
    {6,
     {2466.591347474605, -5183.747588889673, 16389.16230353061, -46312.88371543633, 155164.4275210856, -474039.91843115364, 1634010.4200513815,
      -5177983.017803363, 18144957.834031165, -58709631.144430175, 207889409.35737365, -681497625.6810303, 2430001846.586147},
     {-3.5776236310309915, -3.468310803145069, -3.1671778688838765, -2.6456810705345384, -2.124657249850629, -1.3386732867571225, -0.5987846256777711,
      0.3331896261963086, 1.1632062407151762, 2.057507622364052, 2.7397556255335394, 3.3534806744236327, 3.5720463988330273},
     {169.47204041521195, 406.25380760991, 346.21068737057107, 438.6283010098268, 344.22166335391563, 296.200797989955, 184.41831439779878,
      130.9994727348434, 69.30127884406114, 43.955973276162744, 20.55697579519043, 12.126219096396916, 4.245815580760939},
     -8036322555.409902},
    {7,
     {2481.207632954465, -6362.754202273209, 22392.097685224144, -74333.7589562725, 277976.93596203317, -989115.4775384894, 3811165.6601780397,
      -14020855.882794509, 54979243.89775543, -206230121.2693944, 817845465.3737432, -3106289585.837921, 12413355918.631016, -47549220884.27788,
      191045680946.67596},
     {-4.0495256612745845, -3.98139452964148, -3.689971143453025, -3.2521266564068787, -2.8036908001357377, -2.1079654933633987, -1.4541224003464062,
      -0.5902645099086131, 0.19998245202637974, 1.1254454030104954, 1.9204990063146425, 2.7459258616348183, 3.3415872267876447, 3.895064762797971,
      4.060595420792988},
     {197.07732225190554, 377.0070696389355, 284.9649495687306, 401.8220005267361, 335.2239608203936, 303.2852237358402, 206.5465369156689,
      155.14875963623933, 90.42361889640489, 60.33100753404626, 31.771414703401167, 19.886862083128506, 9.5573315892684, 5.8298144207889635,
      2.331760632977387},
     -736205157513.6132},
    {8,
     {2489.4292935368867, -7568.7095920102365, 29521.675532544123, -112382.35900838731, 466216.6479112017, -1888635.4139253213, 8077486.650411985,
      -33738317.72896783, 146921459.1302053, -624931776.2437028, 2753546379.719804, -11852465691.534496, 52648303122.93557, -228497572831.86853,
      1020884555935.2437, -4457064011188.101, 19998508763676.832},
     {-4.5296048785805425, -4.489120655413785, -4.190900408006874, -3.8262407067061783, -3.4322476333047254, -2.8108861593471337, -2.2303617275329115,
      -1.438018862287728, -0.7084137874907576, 0.1870398298052856, 0.9869479668592022, 1.8842932929899119, 2.6310077616393754, 3.392765063229097,
      3.913110352990816, 4.423933721860695, 4.549698568072851},
     {227.2908132486695, 357.999812342179, 231.5618954053202, 369.3395996433083, 320.1405560766656, 301.8842313814065, 220.11843435881127,
      173.02931298407267, 108.80844036944393, 75.7602753151015, 42.91842007336078, 27.77031211376765, 14.787754748056686, 9.253960065384765,
      4.554512879686469, 2.895651556660825, 1.3153109749910836},
     -87694302012495.61},
    {9,
     {2494.0539776144988, -8791.622156698508, 37815.80683595393, -162164.01832356126, 741526.2405970474, -3363407.7871257607, 15854708.854400141,
      -73977018.9637701, 355105941.01243275, -1685530848.8588831, 8188030951.409323, -39310149685.55352, 192558006507.68842, -931893190046.4677,
      4592457896088.602, -22355798268955.23, 110668743551186.84, -541095464909048.2, 2687799548115842.5},
     {-5.016820175826822, -4.993685431674, -4.676722447784404, -4.378545344787504, -4.025370217771874, -3.4660838695814045, -2.9467259326797097,
      -2.2199043605621736, -1.5516494948405655, -0.7052594554371885, 0.06166414666263719, 0.9614308850526975, 1.7482107742919235, 2.608061097796523,
      3.3028558503986445, 4.0083883440468275, 4.463410463196715, 4.944442238491908, 5.0398627047375495},
     {253.38202650471257, 345.57902693556036, 190.89938974539396, 341.62648692713486, 301.5883298453427, 295.42071632242676, 227.21078892518796,
      185.45830133675568, 124.06709627129003, 89.74594733538459, 54.029185510100376, 36.01637866708152, 20.24266910108437, 12.926177453495349,
      7.0064149116689, 4.412043104806346, 2.2133275146426756, 1.4749101355117848, 0.7547610669180227},
     -1.3185608254453832e+16},
    {10,
     {2496.6553624081557,     -10025.084781938833,   47301.83022268322,    -225421.90071622323, 1129078.9270932479, -5662472.1712942775,
      29222292.60182534,      -150465233.23406377,   790600099.830338,     -4137367824.2340574, 22000210196.705147, -116397507813.75946,
      624158911120.998,       -3328095482910.079,    17956490710833.46,    -96299410310567.08,  521994013148592.2,  -2811681587349518.5,
      1.5295371439514244e+16, -8.26664632916014e+16, 4.509586135115074e+17},
     {-5.509252307091799, -5.496367709027718,  -5.153221439350712,  -4.915659914004084,  -4.5922732159350295, -4.0862495320785435,
      -3.617905753426561, -2.9488613271861492, -2.3369718110156477, -1.5431822067185081, -0.8211547738229005, 0.0517151289031325,
      0.8310461483683911, 1.719123339442227,   2.480471132733749,   3.2998088562801398,  3.9435368493223373,  4.600400175434689,
      4.9983595876088485, 5.459244182977154,   5.531333676033897},
     {272.2293114669295, 337.4404402304291, 164.98323635250028, 318.5264605855846,  280.6655275510083,  286.072006150773,   229.5153378273488,
      193.4010049785687, 136.1809276894156, 101.8992390958499,  64.65932068220015,  44.40084220415069,  26.147676395960712, 17.011921591385146,
      9.595028802541767, 6.106616997844715, 3.376543192274911,  2.1470263273025396, 1.0933275285189314, 0.7662983376218174, 0.43726841994655724},
     -2.443761851595815e+18},
};
const std::vector<SpeedReference> bimodal_references{
    {1, {1.0, 0.0, 3.0625}, {-3.031088913245535, 0.0, 3.031088913245535}, {0.16666666666666666, 0.6666666666666666, 0.16666666666666666}, 0.0},
    {2,
     {1.0, 0.0, 3.125, 6.0, 23.28125},
     {-1.2658808658067997, -1.0516162655934158, 0.6847074631082078, 2.9716162655934157, 3.461173402698592},
     {0.2383761955600233, 0.4431684585802658, 0.09279859863339958, 0.15683154141973416, 0.06882520580657714},
     69.1974},
    {3,
     {1.0, 0.0, 3.1875, 6.0, 24.45703125, 71.25, 245.798583984375},
     {-2.5295742873647553, -1.2767419172736196, -1.006264146579963, 0.11360992199006094, 2.555585838850039, 3.19068020464995, 3.683650207583201},
     {0.010749992001207222, 0.32659364754776465, 0.3029640195393387, 0.11836393446922339, 0.08070444975467322, 0.1264709894115834,
      0.03415296727620943},
     788.3011020350955},
    {4,
     {1.0, 0.0, 3.25, 6.0, 25.65625, 75.0, 268.9140625, 893.8125, 3162.970703125},
     {-1.9212451064758034, -1.6495402017466672, -1.010204580363288, -0.6848151789715347, 0.49090589338277313, 2.532325146448094, 3.0493687168124834,
      3.5768533075965108, 4.109703918301839},
     {0.04478042094146327, 0.14297331019456302, 0.25254622192623366, 0.27080250516655163, 0.045658739440936735, 0.08217628482700444,
      0.09263422981879008, 0.059603455367436435, 0.008824832317020696},
     11097.68850801814},
    {5,
     {1.0, 0.0, 3.3125, 6.0, 26.87890625, 78.75, 293.157470703125, 993.0703125, 3640.607925415039, 13130.302734375, 48949.67107105255},
     {-2.6842245757298424, -1.930515529518236, -1.6153097614496106, -0.9783514553808118, -0.6053652391493506, 0.26542135964086744, 2.198340757445132,
      2.847678554224984, 3.3322504772786123, 3.8669400631789226, 4.2597159321799305},
     {0.0024002432201221256, 0.07315898927982754, 0.13524905443229834, 0.2811632191187097, 0.19950388215177364, 0.061761304661939985,
      0.04461025571241405, 0.0984553887314828, 0.06729592354699147, 0.030915643662585465, 0.005486095481854898},
     183033.68233187468},
    {6,
     {1.0, 0.0, 3.375, 6.0, 28.125, 82.5, 318.55078125, 1097.25, 4161.005859375, 15403.59375, 59430.741943359375, 229576.48828125, 906680.3247070312},
     {-2.5403358710605217, -2.246492101146379, -1.7501889039687621, -1.3652880045997124, -0.8119682322128231, -0.39201495636528155,
      0.4549525240171007, 2.156563889454139, 2.743149025025261, 3.249806409733509, 3.762150666941274, 4.191410925215742, 4.668557980599157},
     {0.005768715349243987, 0.03008386971936519, 0.09829181327240419, 0.20010067001026438, 0.20892102122082068, 0.17087555935977447,
      0.038954633272117706, 0.04655052892115204, 0.07338864580324661, 0.07797190010141987, 0.03506480614873622, 0.01287901034956251,
      0.0011488264718921453},
     3609350.60213757},
    {7,
     {1.0, 0.0, 3.4375, 6.0, 29.39453125, 86.25, 345.115966796875, 1206.3515625, 4726.183242797852, 17913.134765625, 71396.32638454437,
      283950.08615112305, 1159949.9033668637, 4771716.448802948, 19964867.77930809},
     {-3.012321430049305, -2.5356643197065454, -2.1997777562360183, -1.689041592940346, -1.3008710067354003, -0.7537476188291211,
      -0.31666952746409227, 0.41472792320003193, 1.912584404529286, 2.5906932570413677, 3.0727899672346504, 3.5945754965049224, 4.016627926007771,
      4.490167841020274, 4.812449978473776},
     {0.0007125261780804627, 0.013174719286583897, 0.038127105666010505, 0.12796195343837447, 0.17385427986763335, 0.21711494719544808,
      0.1337065171463266, 0.04704668569109789, 0.030899021290147433, 0.0712463488479682, 0.06802649834409379, 0.05123995079301652,
      0.020493247216788503, 0.005548728080844282, 0.0008474709575860406},
     84286662.13421974},
    {8,
     {1.0, 0.0, 3.5, 6.0, 30.6875, 90.0, 372.875, 1320.375, 5338.197265625, 20670.0, 84973.3291015625, 347091.29296875, 1463676.3436279297,
      6206511.03515625, 26844734.510131836, 117140157.34643555, 518339002.90265656},
     {-3.104165257450218, -2.8276403077246584, -2.4159481519934674, -2.0249441519730267, -1.5762736694489283, -1.1445521000298373,
      -0.6481325028147144, -0.17293860548300768, 0.5123026187197909, 1.8590484947479764, 2.504456684685326, 2.9991386580821677, 3.504876857328948,
      3.9400383514696644, 4.401331001333805, 4.788280413868706, 5.16503701671719},
     {0.0007958474239861942, 0.005437885327369055, 0.022524823945033783, 0.069908644744228, 0.12656371086327584, 0.19562937641603464,
      0.17013248219967092, 0.12236997116970755, 0.037753195268465906, 0.03273712199949091, 0.055973240954114864, 0.07216651356229971,
      0.04820455974802166, 0.02888666654514099, 0.008440204773265087, 0.0022755849416115025, 0.00020017011828338575},
     2314551519.924856},
    {9,
     {1.0, 0.0, 3.5625, 6.0, 32.00390625, 93.75, 401.849853515625, 1439.3203125, 5999.143569946289, 23685.263671875, 100294.4586725235,
      419856.6303405762, 1824597.5784035325, 7956866.686420441, 35505583.545894235, 159746474.70682776, 730247516.4005091, 3368991707.626599,
      15730191242.099619},
     {-3.4242240077698645, -3.110566374251313, -2.769772597882723, -2.343475223646785, -1.9607740538296614, -1.4996243161645428, -1.0811883944477407,
      -0.581580708954464, -0.112374663476971, 0.5364983208974893, 1.6980354261533492, 2.3934195854843847, 2.875477498059148, 3.383708041246774,
      3.813262934713826, 4.272778828483964, 4.6550033573187815, 5.078530955930239, 5.339542400079352},
     {0.0002118286522616021, 0.002286096128337332, 0.008901791771846348, 0.03629955351741642, 0.07312054669825159, 0.14665021495887853,
      0.16732792970468907, 0.1720155961681827, 0.10114186709886515, 0.042465235065671854, 0.026131126764522304, 0.05442661315229061,
      0.06040338465648201, 0.056281766401769386, 0.03181572072652736, 0.014956782582040735, 0.0044863212311685614, 0.0009339314990966368,
      0.00014369322170179187},
     74089665129.77542},
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
      4.994731812154773,   5.361743761366499,   5.6353632439529635},
     {0.00014709599789973154, 0.0009636873057217259,  0.004521315306432679,  0.01786423472672289, 0.04286865427605846,  0.0956476056560478,
      0.1345452302005446,     0.17495402472961544,    0.14174633296466185,   0.09808229066226982, 0.037895577311828244, 0.028323217162521404,
      0.044720986779576075,   0.06207945668745603,    0.05069168600469929,   0.03810641747135874, 0.01722114121754575,  0.007400870037467353,
      0.0017867618561584935,  0.00038771937034260827, 4.5694275071015505e-05},
     2741747798830.0146},
};
const std::vector<NormalReference> normal_references{
    {1, 0.75, 1.5, {2.0, 1.5, 5.625}, {0.3333333333333333, 1.3333333333333333, 0.3333333333333333}, 10.96875},
    {2, 0.75, 1.5, {2.0, 1.5, 5.625, 10.96875, 46.1953125}, {0.06666666666666667, 0.6, 0.6666666666666666, 0.6, 0.06666666666666667}, 133.365234375},
    {3,
     0.75,
     1.5,
     {2.0, 1.5, 5.625, 10.96875, 46.1953125, 133.365234375, 619.72119140625},
     {0.012255205375272136, 0.19047619047619047, 0.4163162231961564, 0.7619047619047619, 0.4163162231961564, 0.19047619047619047,
      0.012255205375272136},
     2265.2215576171875},
    {1, -1.25, 0.5, {0.5, -0.625, 0.90625}, {0.08333333333333333, 0.3333333333333333, 0.08333333333333333}, -1.4453125},
    {2,
     -1.25,
     0.5,
     {0.5, -0.625, 0.90625, -1.4453125, 2.486328125},
     {0.016666666666666666, 0.15, 0.16666666666666666, 0.15, 0.016666666666666666},
     -4.55322265625},
    {3,
     -1.25,
     0.5,
     {0.5, -0.625, 0.90625, -1.4453125, 2.486328125, -4.55322265625, 8.7994384765625},
     {0.003063801343818034, 0.047619047619047616, 0.1040790557990391, 0.19047619047619047, 0.1040790557990391, 0.047619047619047616,
      0.003063801343818034},
     -17.829132080078125},
};
const std::vector<BoundaryReference> boundary_references{
    {2, {1.0, 1.75, 6.0625, 15.109375, 43.50390625}, {-1.25, -1.25, 0.75, 2.75, 2.75}, {0.1, 0.15, 0.0, 0.45, 0.3}, 117.1943359375},
    {3,
     {1.0, 0.4375, 2.484375, 4.36328125, 15.5654296875, 37.421142578125, 110.97528076171875},
     {-1.5, -1.5, 0.25, 0.25, 0.5, 2.75, 2.75},
     {0.10714285714285714, 0.14285714285714285, 0.21428571428571427, 0.2857142857142857, 0.0, 0.14285714285714285, 0.10714285714285714},
     293.0791778564453},
    {3,
     {1.0, 1.75, 6.0625, 15.109375, 43.50390625, 117.1943359375, 325.336181640625},
     {-1.25, -1.25, 0.75, 2.75, 2.75},
     {0.1, 0.15, 0.0, 0.45, 0.3},
     890.8598022460938},
    {3, {2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125}, {0.5, 0.5, 0.5}, {1.3333333333333333, 0.6666666666666666, 0.0}, 0.015625},
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

// The node count of the quadrature that the hyperbolic family reconstructs
// from the moments.
std::size_t node_count_of(const std::vector<double>& moments, const std::size_t n)
{
    return qmomkit::internal::compute_hyperbolic_quadrature(qmomkit::MomentSet(moments), n, {}).weights.size();
}

// The speeds, weights and closure moment of the moments with the default
// options, those of a default hyperbolic case.
qmomkit::internal::HyperbolicCharacteristicSpeeds result_of(const std::vector<double>& moments, const std::size_t n)
{
    qmomkit::internal::HyperbolicCharacteristicSpeeds result =
        qmomkit::internal::compute_hyperbolic_characteristic_speeds(qmomkit::MomentSet(moments), n, {});
    QMOMKIT_VERIFY_CHECK(result.status == qmomkit::QuadratureStatus::success);
    return result;
}

std::vector<double> speeds_of(const std::vector<double>& moments, const std::size_t n)
{
    return result_of(moments, n).speeds;
}

// The weights relative to m_0, their sum.
void verify_weights(const std::vector<double>& actual, const std::vector<double>& expected, const double m0, const double tolerance)
{
    QMOMKIT_VERIFY_CHECK(actual.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
    {
        QMOMKIT_VERIFY_NEAR(actual[i], expected[i], tolerance * m0);
    }
}

// The closure moment m_{2n+1} relative to the reference sum of w |speed|^{2n+1}.
void verify_closure_moment(const double actual, const double expected, const std::vector<double>& speeds, const std::vector<double>& weights,
                           const std::size_t n, const double tolerance)
{
    double absolute = 0.0;
    for (std::size_t i = 0; i < speeds.size(); ++i)
    {
        absolute += weights[i] * std::pow(std::abs(speeds[i]), static_cast<double>((2 * n) + 1));
    }
    QMOMKIT_VERIFY_NEAR(actual, expected, tolerance * absolute);
}

// The reference speeds of a normal case: those of the standard normal
// distribution, shifted by the mean and scaled by the deviation.
std::vector<double> normal_speeds(const NormalReference& reference)
{
    std::vector<double> speeds = standard_normal_speeds(reference.n);
    for (double& speed : speeds)
    {
        speed = reference.mean + (reference.deviation * speed);
    }
    return speeds;
}

void test_normal()
{
    for (const NormalReference& reference : normal_references)
    {
        QMOMKIT_VERIFY_CHECK(node_count_of(reference.moments, reference.n) == reference.n + 1);
        const std::vector<double> speeds = speeds_of(reference.moments, reference.n);
        verify_speeds(speeds, normal_speeds(reference), 1.0e-12);
        verify_strictly_ascending(speeds);
    }
}

void test_boundary()
{
    for (const BoundaryReference& reference : boundary_references)
    {
        QMOMKIT_VERIFY_CHECK((2 * node_count_of(reference.moments, reference.n)) + 1 == reference.speeds.size());
        verify_speeds(speeds_of(reference.moments, reference.n), reference.speeds, 1.0e-12);
    }
}

void verify_against_references(const std::vector<SpeedReference>& references)
{
    for (const SpeedReference& reference : references)
    {
        QMOMKIT_VERIFY_CHECK(node_count_of(reference.moments, reference.n) == reference.n + 1);
        const std::vector<double> speeds = speeds_of(reference.moments, reference.n);
        verify_speeds(speeds, reference.speeds, tolerance_for(reference.n));
        verify_strictly_ascending(speeds);
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

void test_weights()
{
    for (const NormalReference& reference : normal_references)
    {
        verify_weights(result_of(reference.moments, reference.n).weights, reference.weights, reference.moments[0], 1.0e-12);
    }
    for (const BoundaryReference& reference : boundary_references)
    {
        verify_weights(result_of(reference.moments, reference.n).weights, reference.weights, reference.moments[0], 1.0e-12);
    }
    for (const std::vector<SpeedReference>* references : {&skewed_references, &bimodal_references})
    {
        for (const SpeedReference& reference : *references)
        {
            verify_weights(result_of(reference.moments, reference.n).weights, reference.weights, reference.moments[0], tolerance_for(reference.n));
        }
    }
}

void test_closure_moment()
{
    // Unlike the errors of the speeds and the weights, those of the closure
    // moment did not grow with n in the observed builds (see above).
    const double closure_tolerance = 1.0e-13;
    for (const NormalReference& reference : normal_references)
    {
        verify_closure_moment(result_of(reference.moments, reference.n).closure_moment, reference.closure_moment, normal_speeds(reference),
                              reference.weights, reference.n, closure_tolerance);
    }
    for (const BoundaryReference& reference : boundary_references)
    {
        verify_closure_moment(result_of(reference.moments, reference.n).closure_moment, reference.closure_moment, reference.speeds, reference.weights,
                              reference.n, closure_tolerance);
    }
    for (const std::vector<SpeedReference>* references : {&skewed_references, &bimodal_references})
    {
        for (const SpeedReference& reference : *references)
        {
            verify_closure_moment(result_of(reference.moments, reference.n).closure_moment, reference.closure_moment, reference.speeds,
                                  reference.weights, reference.n, closure_tolerance);
        }
    }
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
                                                 {"weights", test_weights},
                                                 {"closure_moment", test_closure_moment},
                                             });
}
