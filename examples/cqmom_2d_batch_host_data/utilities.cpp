#include "utilities.hpp"

#include "physics.hpp"

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace cqmom_2d_batch_host_data
{

std::filesystem::path location_output_file(const std::filesystem::path& directory, const std::size_t location, const std::string& name)
{
    return directory / ("location_" + std::to_string(location) + "_" + name);
}

void print_location_summary(const LocalHostData& host_data, const std::vector<qmomkit::CqmomState>& states, const double initial_time,
                            const double initial_number_density, const double initial_mean_velocity, const double time)
{
    std::cout << std::setprecision(6);
    std::cout << "location  x      lambda   M[0,0]        exact M[0,0]  mean u        exact mean u\n";
    for (std::size_t location = 0; location < states.size(); ++location)
    {
        const qmomkit::CqmomState& state = states[location];
        const double number_density = state.moment({0, 0}).value();
        const double mean_velocity = state.moment({0, 1}).value() / number_density;
        std::cout << std::left << std::setw(10) << location << std::setw(7) << host_data.position[location] << std::setw(9)
                  << host_data.fragmentation_rate[location] << std::setw(14) << number_density << std::setw(14)
                  << exact_number_density(host_data.fragmentation_rate[location], initial_time, initial_number_density, time) << std::setw(14)
                  << mean_velocity << exact_mean_velocity(host_data.position[location], initial_time, initial_mean_velocity, time) << '\n';
    }
    std::cout << std::right;
}

} // namespace cqmom_2d_batch_host_data
