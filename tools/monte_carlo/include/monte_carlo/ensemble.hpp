#pragma once

#include <cstddef>
#include <vector>

namespace qmomkit::monte_carlo
{

/**
 * One realization drawn from a modeled phase-space distribution.
 *
 * Coordinates are stored in QMOMKit phase-space dimension order. The sample is
 * a mathematical object; a physical application may interpret it as a particle,
 * parcel, droplet, or another entity, but the Monte Carlo companion package
 * assigns no application-specific meaning to the coordinates.
 */
using Sample = std::vector<double>;

/**
 * Finite stochastic realization of a distribution.
 *
 * Unless an API states otherwise, samples in an ensemble are equally weighted
 * when estimating moments. The total number density is supplied separately and
 * scales the ensemble average into QMOMKit's absolute raw-moment convention.
 */
using Ensemble = std::vector<Sample>;

struct SamplingOptions
{
    // Number of independent samples to draw for the finite realization.
    std::size_t sample_count = 0;
    // Seed for the deterministic pseudo-random engine used by this package.
    unsigned long long seed = 0;
};

} // namespace qmomkit::monte_carlo
