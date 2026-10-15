#include "utilities.hpp"

#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>

namespace qmom_1d_breakage_advection_with_mc
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

void write_particle_history(const std::filesystem::path& path, const std::vector<std::tuple<double, std::size_t, double>>& history)
{
    std::ofstream file(path);
    if (!file)
    {
        throw std::runtime_error("Could not write " + path.string());
    }
    file << "# time particle_id volume\n";
    for (const auto& [time, particle_id, volume] : history)
    {
        file << std::setprecision(12) << time << std::setprecision(17) << ' ' << particle_id << ' ' << volume << '\n';
    }
}

} // namespace qmom_1d_breakage_advection_with_mc
