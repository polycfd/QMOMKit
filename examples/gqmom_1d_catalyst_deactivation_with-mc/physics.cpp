#include "physics.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace gqmom_1d_catalyst_deactivation_with_mc
{

qmomkit::PhaseSpaceVelocity coking_velocity()
{
    return [](double, const std::vector<double>& coked_fractions, std::vector<double>& velocities, std::size_t)
    {
        for (std::size_t node = 0; node < coked_fractions.size(); ++node)
        {
            velocities[node] = coking_rate * (1.0 - coked_fractions[node]);
        }
    };
}

double removal_frequency(const double coked_fraction)
{
    return removal_coefficient * std::max(0.0, coked_fraction - critical_fraction);
}

qmomkit::QmomMomentSource removal_source()
{
    return [](double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> source(context.size(), 0.0);
        // A particle of coked fraction theta leaves at the frequency g(theta), so
        // a node of weight w adds -g(theta) w theta^k to dM_k/dt. g is zero below
        // the critical fraction, so the rate is not closed in the moments. Summing
        // over the nodes of the current stage quadrature closes it.
        context.for_each_node(
            [&source](const double weight, const double coked_fraction)
            {
                const double frequency = removal_frequency(coked_fraction);
                double power = 1.0;
                for (double& rate : source)
                {
                    rate -= frequency * weight * power;
                    power *= coked_fraction;
                }
            });
        return qmomkit::QmomMomentRate{std::move(source)};
    };
}

void apply_removal(qmomkit::monte_carlo::Ensemble& particles, const double step_size, std::mt19937_64& generator)
{
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    qmomkit::monte_carlo::Ensemble remaining;
    remaining.reserve(particles.size());
    for (const qmomkit::monte_carlo::Sample& particle : particles)
    {
        // Particles below the critical fraction have zero frequency and stay.
        const double frequency = removal_frequency(particle[0]);
        if (frequency > 0.0 && uniform(generator) < 1.0 - std::exp(-frequency * step_size))
        {
            continue;
        }
        remaining.push_back(particle);
    }
    particles = std::move(remaining);
}

} // namespace gqmom_1d_catalyst_deactivation_with_mc
