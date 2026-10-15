#include "physics.hpp"

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace qmom_1d_breakage_advection
{
namespace
{

// Fragmentation rate a(V) = a_ref * (V / V_ref)^p: larger particles break
// more often.
constexpr double reference_fragmentation_rate = 0.5;
constexpr double reference_volume = 1.0;
constexpr double rate_exponent = 1.5;

} // namespace

qmomkit::PhaseSpaceVelocity advection_velocity()
{
    // Volumes relax toward V = 1: smaller particles grow, larger ones shrink.
    return [](double, const std::vector<double>& volumes, std::vector<double>& velocities, std::size_t)
    {
        for (std::size_t node = 0; node < volumes.size(); ++node)
        {
            velocities[node] = 1.0 - volumes[node];
        }
    };
}

qmomkit::QmomMomentSource binary_fragmentation_source()
{
    return [](double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> source(context.size(), 0.0);
        // Binary fragmentation V -> V/2 at the rate a(V) replaces one particle
        // of volume V by two of volume V/2, so a node of weight w adds
        // a(V) * w * (2^(1-k) - 1) * V^k to dM_k/dt. The exact rate involves
        // the moment of order k + p, which is not transported; summing over
        // the nodes of the current stage quadrature closes it.
        context.for_each_node(
            [&source](const double weight, const double volume)
            {
                const double rate = reference_fragmentation_rate * std::pow(volume / reference_volume, rate_exponent);
                for (std::size_t order = 0; order < source.size(); ++order)
                {
                    const double daughter_factor = std::pow(2.0, 1.0 - static_cast<double>(order));
                    source[order] += rate * weight * (daughter_factor - 1.0) * std::pow(volume, static_cast<double>(order));
                }
            });
        return qmomkit::QmomMomentRate{std::move(source)};
    };
}

} // namespace qmom_1d_breakage_advection
