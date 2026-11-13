#include "physics.hpp"

#include <array>
#include <cmath>
#include <utility>
#include <vector>

namespace gqmom_1d_absolute_value_loss
{
namespace
{

// Loss-rate coefficient c of the rate c |xi|.
constexpr double loss_coefficient = 1.0;
constexpr double pi = 3.14159265358979323846;

// The moments of order 0, 1 and 2 of a normal density of unit number density
// with the given mean and standard deviation, restricted to xi > 0 or to
// xi < 0. Integrating xi^(k-1) (xi - mean) by parts gives
// I_k = mean I_(k-1) + sd^2 ((k - 1) I_(k-2) +- [k = 1] density(0)).
std::array<double, 3> half_line_moments(const double mean, const double standard_deviation, const bool positive)
{
    const double variance = standard_deviation * standard_deviation;
    const double density_at_zero = std::exp(-mean * mean / (2.0 * variance)) / std::sqrt(2.0 * pi * variance);
    const double z = mean / (standard_deviation * std::sqrt(2.0));
    const double sign = positive ? 1.0 : -1.0;
    std::array<double, 3> moments{};
    moments[0] = 0.5 * std::erfc(-sign * z);
    moments[1] = (mean * moments[0]) + (sign * variance * density_at_zero);
    moments[2] = (mean * moments[1]) + (variance * moments[0]);
    return moments;
}

} // namespace

qmomkit::QmomMomentSource absolute_value_loss()
{
    return [](double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> source(context.size(), 0.0);
        // An element at xi is lost at the rate c |xi|, so a node of weight w
        // adds -c |xi| w xi^k to dM_k/dt. |xi| is not a polynomial, so the rate
        // is not closed in the moments. Summing over the nodes of the current
        // stage quadrature closes it.
        context.for_each_node(
            [&source](const double weight, const double abscissa)
            {
                double power = 1.0;
                for (double& rate : source)
                {
                    rate -= loss_coefficient * std::abs(abscissa) * weight * power;
                    power *= abscissa;
                }
            });
        return qmomkit::QmomMomentRate{std::move(source)};
    };
}

ExactMoments exact_moments(const double number_density, const double mean, const double standard_deviation, const double time)
{
    // On each half line, n(0, xi) exp(-+a xi) with a = c t is a normal density of
    // the same standard deviation with the mean shifted by -+a sd^2, scaled by
    // exp(-+a mean + a^2 sd^2 / 2).
    const double a = loss_coefficient * time;
    const double variance = standard_deviation * standard_deviation;
    const double positive_scale = std::exp((-a * mean) + (0.5 * a * a * variance));
    const double negative_scale = std::exp((a * mean) + (0.5 * a * a * variance));
    const std::array<double, 3> positive = half_line_moments(mean - (a * variance), standard_deviation, true);
    const std::array<double, 3> negative = half_line_moments(mean + (a * variance), standard_deviation, false);
    return ExactMoments{number_density * ((positive_scale * positive[0]) + (negative_scale * negative[0])),
                        number_density * ((positive_scale * positive[2]) + (negative_scale * negative[2]))};
}

} // namespace gqmom_1d_absolute_value_loss
