#include "utilities.hpp"

#include <cstddef>
#include <iomanip>
#include <ostream>

namespace gqmom_1d_droplet_breakup_with_mc
{

void print_moments(std::ostream& stream, const std::string_view run, const std::vector<double>& initial_moments,
                   const std::vector<double>& final_moments)
{
    stream << run << ", raw moments of the volume distribution:\n";
    stream << "      " << std::setw(14) << "initial" << std::setw(14) << "final" << '\n';
    for (std::size_t order = 0; order < initial_moments.size(); ++order)
    {
        stream << "  M" << order << "  " << std::setw(14) << initial_moments[order] << std::setw(14) << final_moments[order] << '\n';
    }
}

} // namespace gqmom_1d_droplet_breakup_with_mc
