#include "particles_mc.hpp"

#include "physics.hpp"

#include <cmath>
#include <random>
#include <stdexcept>

namespace qmom_1d_breakage_advection_mc_local_rates
{

std::vector<Particle> make_particles(const qmomkit::monte_carlo::Ensemble& ensemble)
{
    std::vector<Particle> particles;
    particles.reserve(ensemble.size());
    for (std::size_t sample_index = 0; sample_index < ensemble.size(); ++sample_index)
    {
        if (ensemble[sample_index].size() != 1)
        {
            throw std::runtime_error("The local-rate MC reference expects one-dimensional samples.");
        }
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
        ensemble.push_back(qmomkit::monte_carlo::Sample{particle.diameter});
    }
    return ensemble;
}

void update_particle_diameters(std::vector<Particle>& particles, const qmomkit::monte_carlo::Ensemble& ensemble)
{
    if (particles.size() != ensemble.size())
    {
        throw std::runtime_error("Advection MC evolution changed the particle count unexpectedly.");
    }

    for (std::size_t index = 0; index < particles.size(); ++index)
    {
        if (ensemble[index].size() != 1)
        {
            throw std::runtime_error("The evolved local-rate MC ensemble is not one-dimensional.");
        }
        particles[index].diameter = ensemble[index][0];
    }
}

std::size_t apply_binary_breakage(std::vector<Particle>& particles, const double step_size, std::mt19937_64& generator, std::size_t& next_particle_id)
{
    std::bernoulli_distribution breaks(breakage_probability(step_size));
    std::size_t event_count = 0;

    const std::size_t initial_particle_count = particles.size();
    particles.reserve(particles.size() + initial_particle_count);
    for (std::size_t particle_index = 0; particle_index < initial_particle_count; ++particle_index)
    {
        if (!breaks(generator))
        {
            continue;
        }

        particles[particle_index].diameter = daughter_diameter(particles[particle_index].diameter);
        particles.push_back(Particle{next_particle_id, particles[particle_index].diameter});
        ++next_particle_id;
        ++event_count;
    }
    return event_count;
}

double expected_nucleation_particle_count(const double step_size, const double particle_weight)
{
    if (step_size < 0.0 || !std::isfinite(step_size))
    {
        throw std::invalid_argument("Nucleation sampling requires a non-negative finite step size.");
    }
    if (!(particle_weight > 0.0) || !std::isfinite(particle_weight))
    {
        throw std::invalid_argument("Nucleation sampling requires a positive finite particle weight.");
    }
    return nucleation_rate() * step_size / particle_weight;
}

std::size_t apply_nucleation(std::vector<Particle>& particles, const double step_size, const double particle_weight, std::mt19937_64& generator,
                             std::size_t& next_particle_id)
{
    std::poisson_distribution<std::size_t> births(expected_nucleation_particle_count(step_size, particle_weight));
    const std::size_t birth_count = births(generator);

    particles.reserve(particles.size() + birth_count);
    for (std::size_t birth = 0; birth < birth_count; ++birth)
    {
        particles.push_back(Particle{next_particle_id, nucleation_diameter()});
        ++next_particle_id;
    }
    return birth_count;
}

std::vector<double> compute_particle_moments(const std::vector<Particle>& particles, const double particle_weight, const std::size_t maximum_order)
{
    std::vector<double> moments(maximum_order + 1, 0.0);
    for (const Particle& particle : particles)
    {
        double power = 1.0;
        for (std::size_t order = 0; order < moments.size(); ++order)
        {
            moments[order] += particle_weight * power;
            power *= particle.diameter;
        }
    }
    return moments;
}

} // namespace qmom_1d_breakage_advection_mc_local_rates
