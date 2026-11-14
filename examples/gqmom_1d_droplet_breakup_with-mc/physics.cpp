#include "physics.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace gqmom_1d_droplet_breakup_with_mc
{

double breakup_frequency(const double volume)
{
    return breakup_coefficient * std::max(0.0, volume - critical_volume);
}

qmomkit::QmomMomentSource binary_breakup_source()
{
    return [](double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> source(context.size(), 0.0);
        // A droplet of volume V breaks into two droplets of volume V/2 at the
        // frequency g(V), so a node of weight w adds g(V) w (2 (V/2)^k - V^k) to
        // dM_k/dt. g is zero below the critical volume, so the rate is not closed
        // in the moments. Summing over the nodes of the current stage quadrature
        // closes it.
        context.for_each_node(
            [&source](const double weight, const double volume)
            {
                const double frequency = breakup_frequency(volume);
                double power = 1.0;
                for (std::size_t order = 0; order < source.size(); ++order)
                {
                    source[order] += frequency * weight * (std::pow(2.0, 1.0 - static_cast<double>(order)) - 1.0) * power;
                    power *= volume;
                }
            });
        return qmomkit::QmomMomentRate{std::move(source)};
    };
}

} // namespace gqmom_1d_droplet_breakup_with_mc
