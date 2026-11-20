#include "particles_mc.hpp"

#include "monte_carlo/moments.hpp"
#include "physics.hpp"

#include <cmath>
#include <stdexcept>

namespace qmom_1d_breakage_advection_with_mc
{

std::vector<Particle> make_particles(const qmomkit::monte_carlo::Ensemble& ensemble)
{
    std::vector<Particle> particles;
    particles.reserve(ensemble.size());
    for (std::size_t sample_index = 0; sample_index < ensemble.size(); ++sample_index)
    {
        particles.push_back(Particle{sample_index, ensemble[sample_index][0]});
    }
    return particles;
}

qmomkit::monte_carlo::Ensemble to_ensemble(const std::vector<Particle>& particles)
{
    qmomkit::monte_carlo::Ensemble ensemble;
    ensemble.reserve(particles.size());
    for (const Particle& particle : particles)
    {
        ensemble.push_back(qmomkit::monte_carlo::Sample{particle.volume});
    }
    return ensemble;
}

void update_particle_volumes(std::vector<Particle>& particles, const qmomkit::monte_carlo::Ensemble& ensemble)
{
    for (std::size_t index = 0; index < particles.size(); ++index)
    {
        particles[index].volume = ensemble[index][0];
    }
}

void apply_binary_breakage(std::vector<Particle>& particles, const double step_size, std::mt19937_64& generator, std::size_t& next_particle_id)
{
    const double breakage_probability = 1.0 - std::exp(-fragmentation_rate * step_size);
    std::bernoulli_distribution breaks(breakage_probability);

    const std::size_t initial_particle_count = particles.size();
    particles.reserve(particles.size() + initial_particle_count);
    for (std::size_t particle_index = 0; particle_index < initial_particle_count; ++particle_index)
    {
        if (!breaks(generator))
        {
            continue;
        }

        particles[particle_index].volume *= 0.5;
        particles.push_back(Particle{next_particle_id, particles[particle_index].volume});
        ++next_particle_id;
    }
}

std::vector<double> compute_particle_moments(const std::vector<Particle>& particles, const double particle_weight, const std::size_t maximum_order)
{
    const double number_density = particle_weight * static_cast<double>(particles.size());
    const qmomkit::monte_carlo::OneDimensionalMomentEstimate estimate =
        qmomkit::monte_carlo::compute_ensemble_moments(to_ensemble(particles), number_density, maximum_order);
    if (!estimate.ok())
    {
        throw std::runtime_error(estimate.message);
    }
    return estimate.moments.values();
}

} // namespace qmom_1d_breakage_advection_with_mc
