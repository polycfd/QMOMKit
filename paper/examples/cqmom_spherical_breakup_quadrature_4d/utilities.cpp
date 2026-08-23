#include "utilities.hpp"

#include "monte_carlo/moments.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/multi_index.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <system_error>
#include <utility>

namespace cqmom_spherical_breakup_quadrature_4d
{
namespace
{

constexpr std::size_t max_console_moments = 20;

const std::vector<qmomkit::MultiIndex>& paper_moment_indices()
{
    static const std::vector<qmomkit::MultiIndex> indices{{0, 0, 0, 0}, {3, 0, 0, 0}, {1, 0, 0, 0}, {2, 0, 0, 0}, {0, 1, 0, 0},
                                                          {0, 2, 0, 0}, {2, 1, 0, 0}, {2, 2, 0, 0}, {6, 0, 0, 0}};
    return indices;
}

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

std::vector<double> paper_moments(const qmomkit::CqmomQuadrature& quadrature)
{
    const qmomkit::SelectedRawMomentResult result = qmomkit::reconstruct_raw_moments(quadrature, {paper_moment_indices()});
    if (!result.ok())
    {
        throw std::runtime_error("CQMOM paper moment reconstruction failed: " + result.message);
    }
    return result.values;
}

std::vector<double> paper_moments(const qmomkit::monte_carlo::Ensemble& ensemble, const double particle_weight)
{
    std::array<double, 9> sums{};
    for (const qmomkit::monte_carlo::Sample& sample : ensemble)
    {
        if (sample.size() != 4)
        {
            throw std::invalid_argument("4D MC samples must be ordered as (d,u,v,w).");
        }
        const double d = sample[0];
        const double u = sample[1];
        const double d2 = d * d;
        sums[0] += 1.0;
        sums[1] += d2 * d;
        sums[2] += d;
        sums[3] += d2;
        sums[4] += u;
        sums[5] += u * u;
        sums[6] += d2 * u;
        sums[7] += d2 * u * u;
        sums[8] += d2 * d2 * d2;
    }
    std::vector<double> values(sums.begin(), sums.end());
    for (double& value : values)
    {
        value *= particle_weight;
    }
    return values;
}

const std::vector<std::string>& paper_moment_labels()
{
    static const std::vector<std::string> labels{"M_0000", "M_3000", "M_1000", "M_2000", "M_0100", "M_0200", "M_2100", "M_2200", "M_6000"};
    return labels;
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

} // namespace cqmom_spherical_breakup_quadrature_4d
