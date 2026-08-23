#include "utilities.hpp"

#include "monte_carlo/moments.hpp"
#include "qmomkit/multi_index.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <system_error>
#include <utility>

namespace cqmom_thermal_spherical_breakup
{
namespace
{

constexpr std::size_t max_console_moments = 20;

void print_multi_index(const qmomkit::MultiIndex& order)
{
    std::cout << '[';
    for (std::size_t dimension = 0; dimension < order.size(); ++dimension)
    {
        if (dimension > 0)
        {
            std::cout << ',';
        }
        std::cout << order[dimension];
    }
    std::cout << ']';
}

} // namespace

std::vector<std::size_t> maximum_orders(const std::vector<std::size_t>& extents)
{
    std::vector<std::size_t> orders;
    orders.reserve(extents.size());
    for (const std::size_t extent : extents)
    {
        if (extent == 0)
        {
            throw std::invalid_argument("Moment extents must be non-empty.");
        }
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

double maximum_tensor_difference(const qmomkit::RawMomentTensor& left, const qmomkit::RawMomentTensor& right)
{
    if (left.extents != right.extents || left.values.size() != right.values.size())
    {
        throw std::invalid_argument("Cannot compare moment tensors with different extents.");
    }
    double difference = 0.0;
    for (std::size_t index = 0; index < left.values.size(); ++index)
    {
        difference = std::max(difference, std::abs(left.values[index] - right.values[index]));
    }
    return difference;
}

qmomkit::SelectedRawMomentResult select_moments(const qmomkit::RawMomentTensor& moments, const qmomkit::RawMomentRequest& request)
{
    std::vector<double> values;
    values.reserve(request.indices.size());
    for (const qmomkit::MultiIndex& index : request.indices)
    {
        if (index.size() != moments.extents.size())
        {
            throw std::invalid_argument("Requested raw moment rank must match the Monte Carlo tensor.");
        }
        for (std::size_t dimension = 0; dimension < index.size(); ++dimension)
        {
            if (index[dimension] >= moments.extents[dimension])
            {
                throw std::out_of_range("Requested raw moment is outside the Monte Carlo tensor.");
            }
        }
        values.push_back(moments.values[qmomkit::flatten_multi_index(index, moments.extents)]);
    }
    return qmomkit::make_selected_raw_moment_result(request, std::move(values));
}

void print_selected_moments(const char* label, const qmomkit::SelectedRawMomentResult& moments)
{
    if (!moments.ok())
    {
        throw std::runtime_error("Could not select tracked raw moments: " + moments.message);
    }
    std::cout << std::setprecision(12) << label << '\n';
    const std::size_t printed_count = std::min(max_console_moments, moments.values.size());
    for (std::size_t index = 0; index < printed_count; ++index)
    {
        std::cout << "M";
        print_multi_index(moments.indices[index]);
        std::cout << " = " << moments.values[index] << '\n';
    }
    if (moments.values.size() > printed_count)
    {
        std::cout << "... " << (moments.values.size() - printed_count) << " more tracked moments written to the history file\n";
    }
}

void ensure_output_directory(const std::filesystem::path& path)
{
    std::error_code error;
    std::filesystem::create_directories(path, error);
    if (error)
    {
        throw std::runtime_error("Could not create output directory " + path.string() + ": " + error.message());
    }
}

} // namespace cqmom_thermal_spherical_breakup
