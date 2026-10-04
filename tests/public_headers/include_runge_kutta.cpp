#include <qmomkit/runge_kutta.hpp>

// Standard headers follow the public header, so they cannot supply an include it is missing.
#include <functional>
#include <type_traits>
#include <vector>

// The right-hand side fills a buffer the integrator supplies rather than
// returning one, which is what keeps the stage evaluations of an integration
// free of heap allocation. Pinning the shape here means a change back to a
// returning form cannot pass unnoticed.
static_assert(std::is_same_v<qmomkit::OdeRhs, std::function<void(double, const std::vector<double>&, std::vector<double>&)>>);

// A returning callable must not silently bind: std::function would happily
// discard the return value if the parameter lists matched, so the extra buffer
// parameter is what makes the old shape a compile error at every call site.
static_assert(!std::is_convertible_v<std::vector<double> (*)(double, const std::vector<double>&), qmomkit::OdeRhs>);

using SspRk3Step = qmomkit::SspRk3Result (*)(const std::vector<double>&, double, double, const qmomkit::OdeRhs&, const qmomkit::SspRk3Options&);
using IntegrateEmbeddedRk = qmomkit::AdaptiveRkResult (*)(const std::vector<double>&, double, double, const qmomkit::OdeRhs&,
                                                          qmomkit::EmbeddedRkMethod, const qmomkit::AdaptiveRkOptions&);

static_assert(std::is_same_v<decltype(&qmomkit::ssp_rk3_step), SspRk3Step>);
static_assert(std::is_same_v<decltype(&qmomkit::integrate_embedded_rk), IntegrateEmbeddedRk>);
