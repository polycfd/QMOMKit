#include "utilities.hpp"

#include <cstddef>
#include <iomanip>
#include <ostream>

namespace cqmom_3d_droplets_condensation_nucleation
{

qmomkit::RawMomentRequest printed_moments()
{
    return {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {2, 0, 0}, {0, 2, 0}, {0, 0, 2}}};
}

void print_moments(std::ostream& stream, const std::string_view run, const std::vector<double>& initial_values,
                   const std::vector<double>& final_values)
{
    const qmomkit::RawMomentRequest request = printed_moments();
    stream << run << ", raw moments M[i,j,k] of r^i w^j T^k:\n";
    stream << "            " << std::setw(14) << "initial" << std::setw(14) << "final" << '\n';
    for (std::size_t entry = 0; entry < request.indices.size(); ++entry)
    {
        const qmomkit::MultiIndex& order = request.indices[entry];
        stream << "  M[" << order[0] << ',' << order[1] << ',' << order[2] << "]  " << std::setw(14) << initial_values[entry] << std::setw(14)
               << final_values[entry] << '\n';
    }
}

} // namespace cqmom_3d_droplets_condensation_nucleation
