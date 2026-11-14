#include "particles_mc.hpp"

#include "physics.hpp"

#include <cmath>
#include <cstddef>

namespace gqmom_1d_droplet_breakup_with_mc
{

void apply_binary_breakup(qmomkit::monte_carlo::Ensemble& droplets, const double step_size, std::mt19937_64& generator)
{
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const std::size_t initial_count = droplets.size();
    for (std::size_t index = 0; index < initial_count; ++index)
    {
        // Droplets below the critical volume have zero frequency and never break.
        const double frequency = breakup_frequency(droplets[index][0]);
        if (frequency == 0.0 || uniform(generator) >= 1.0 - std::exp(-frequency * step_size))
        {
            continue;
        }
        droplets[index][0] *= 0.5;
        droplets.push_back(qmomkit::monte_carlo::Sample{droplets[index][0]});
    }
}

} // namespace gqmom_1d_droplet_breakup_with_mc
