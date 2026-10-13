#include "utilities.hpp"

#include "monte_carlo/moments.hpp"

#include <iomanip>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cqmom_5d_particles_breakage_heat_transfer
{

std::vector<std::size_t> maximum_orders(const std::vector<std::size_t>& extents)
{
    std::vector<std::size_t> orders;
    orders.reserve(extents.size());
    for (const std::size_t extent : extents)
    {
        orders.push_back(extent - 1);
    }
    return orders;
}

qmomkit::RawMomentTensor ensemble_moments(const qmomkit::monte_carlo::Ensemble& ensemble, const double particle_weight,
                                          const std::vector<std::size_t>& extents)
{
    const qmomkit::monte_carlo::TensorMomentEstimate estimate =
        qmomkit::monte_carlo::compute_ensemble_moments(ensemble, particle_weight * static_cast<double>(ensemble.size()), maximum_orders(extents));
    if (!estimate.ok())
    {
        throw std::runtime_error(estimate.message);
    }
    return estimate.moments;
}

qmomkit::RawMomentRequest printed_moments()
{
    return {{{0, 0, 0, 0, 0},
             {1, 0, 0, 0, 0},
             {0, 1, 0, 0, 0},
             {0, 0, 1, 0, 0},
             {0, 0, 0, 1, 0},
             {0, 0, 0, 0, 1},
             {2, 0, 0, 0, 0},
             {0, 2, 0, 0, 0},
             {0, 0, 2, 0, 0},
             {0, 0, 0, 2, 0},
             {0, 0, 0, 0, 2},
             {1, 1, 0, 0, 0}}};
}

void print_moments(std::ostream& stream, const std::string_view run, const std::vector<double>& initial_values,
                   const std::vector<double>& final_values)
{
    const qmomkit::RawMomentRequest request = printed_moments();
    stream << run << ", raw moments M[i,j,k,l,m] of V^i T^j u^k v^l w^m:\n";
    stream << "                " << std::setw(14) << "initial" << std::setw(14) << "final" << '\n';
    for (std::size_t entry = 0; entry < request.indices.size(); ++entry)
    {
        std::ostringstream label;
        label << "M[" << request.indices[entry][0] << ',' << request.indices[entry][1] << ',' << request.indices[entry][2] << ','
              << request.indices[entry][3] << ',' << request.indices[entry][4] << ']';
        stream << "  " << std::left << std::setw(14) << label.str() << std::right << std::setw(14) << initial_values[entry] << std::setw(14)
               << final_values[entry] << '\n';
    }
}

} // namespace cqmom_5d_particles_breakage_heat_transfer
