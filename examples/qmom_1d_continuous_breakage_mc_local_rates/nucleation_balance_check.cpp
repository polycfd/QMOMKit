#include "particles_mc.hpp"
#include "physics.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

int main()
{
    const double tolerance = 1.0e-14;
    const double particle_weight = 1.0e-4;
    const double source_step = 0.25;
    const std::vector<double> moments(6, 0.0);

    const qmomkit::MomentBoundaryContribution boundary = qmom_1d_continuous_breakage_mc_local_rates::boundary_contribution();
    const qmomkit::MomentBoundaryRate boundary_rate = boundary(0.0, moments);
    if (boundary_rate.values.size() != moments.size())
    {
        std::cerr << "Nucleation boundary contribution returned the wrong moment count.\n";
        return EXIT_FAILURE;
    }

    const double expected_birth_particles =
        qmom_1d_continuous_breakage_mc_local_rates::expected_nucleation_particle_count(source_step, particle_weight);
    for (std::size_t order = 0; order < boundary_rate.values.size(); ++order)
    {
        const double expected_rate = qmom_1d_continuous_breakage_mc_local_rates::nucleation_rate() *
                                     std::pow(qmom_1d_continuous_breakage_mc_local_rates::nucleation_diameter(), static_cast<double>(order));
        const double expected_mc_increment = expected_birth_particles * particle_weight *
                                             std::pow(qmom_1d_continuous_breakage_mc_local_rates::nucleation_diameter(), static_cast<double>(order));
        const double expected_moment_increment = expected_rate * source_step;
        if (std::abs(boundary_rate.values[order] - expected_rate) > tolerance)
        {
            std::cerr << "QMOM nucleation boundary rate does not match J_nuc d_nuc^k.\n";
            return EXIT_FAILURE;
        }
        if (std::abs(expected_mc_increment - expected_moment_increment) > tolerance)
        {
            std::cerr << "MC nucleation intensity does not match the QMOM moment increment.\n";
            return EXIT_FAILURE;
        }
    }

    const double diameter = 2.0;
    const qmomkit::MomentSource source = qmom_1d_continuous_breakage_mc_local_rates::binary_fragmentation_source();
    const std::vector<double> single_particle_moments{1.0, diameter, diameter * diameter};
    const std::vector<double> breakage_rates = source(0.0, single_particle_moments);
    if (breakage_rates.size() != single_particle_moments.size() || std::abs(breakage_rates[1]) > tolerance)
    {
        std::cerr << "Binary breakage should conserve the first raw moment for one particle.\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
