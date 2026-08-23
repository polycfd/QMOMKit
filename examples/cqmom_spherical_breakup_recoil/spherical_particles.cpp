#include "spherical_particles.hpp"

#include "monte_carlo/moments.hpp"
#include "qmomkit/multi_index.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>

namespace cqmom_spherical_breakup_recoil
{
namespace
{

constexpr double pi = 3.14159265358979323846;
constexpr double gas_density = 1.2;
constexpr double gas_viscosity = 1.8e-5;
constexpr double particle_density = 2000.0;
constexpr double gas_u = 2.0;
constexpr double gas_v = 1.5;
constexpr double gas_w = 0.8;
constexpr double reference_diameter = 1.0e-4;
constexpr double reference_speed = 2.5;
constexpr double breakup_frequency_scale = 0.1;
constexpr double recoil_standard_deviation = 0.001;
constexpr double daughter_diameter_factor = 0.79370052598409973738;

double integer_power(const double value, const std::size_t exponent)
{
    double result = 1.0;
    for (std::size_t i = 0; i < exponent; ++i)
    {
        result *= value;
    }
    return result;
}

double binomial(const std::size_t n, std::size_t k)
{
    k = std::min(k, n - k);
    double result = 1.0;
    for (std::size_t i = 1; i <= k; ++i)
    {
        result *= static_cast<double>(n - k + i) / static_cast<double>(i);
    }
    return result;
}

double shifted_gaussian_moment(const double mean, const std::size_t order)
{
    double value = 0.0;
    for (std::size_t gaussian_order = 0; gaussian_order <= order; gaussian_order += 2)
    {
        double centered_moment = 1.0;
        for (std::size_t factor = 1; factor < gaussian_order; factor += 2)
        {
            centered_moment *= static_cast<double>(factor) * recoil_standard_deviation * recoil_standard_deviation;
        }
        value += binomial(order, gaussian_order) * integer_power(mean, order - gaussian_order) * centered_moment;
    }
    return value;
}

double breakup_frequency(const double d, const double u, const double v, const double w)
{
    const double du = gas_u - u;
    const double dv = gas_v - v;
    const double dw = gas_w - w;
    return breakup_frequency_scale * d / reference_diameter * (du * du + dv * dv + dw * dw) / (reference_speed * reference_speed);
}

double parent_monomial(const std::vector<double>& coordinates, const qmomkit::MultiIndex& order)
{
    double value = 1.0;
    for (std::size_t dimension = 0; dimension < order.size(); ++dimension)
    {
        value *= integer_power(coordinates[dimension], order[dimension]);
    }
    return value;
}

double expected_daughter_monomial(const std::vector<double>& coordinates, const qmomkit::MultiIndex& order)
{
    return integer_power(daughter_diameter_factor * coordinates[0], order[0]) * shifted_gaussian_moment(coordinates[1], order[1]) *
           shifted_gaussian_moment(coordinates[2], order[2]) * shifted_gaussian_moment(coordinates[3], order[3]);
}

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

double tensor_value(const qmomkit::RawMomentTensor& moments, const qmomkit::MultiIndex& index)
{
    return moments.values[qmomkit::flatten_multi_index(index, moments.extents)];
}

} // namespace

std::optional<std::string> read_configuration_text(const std::filesystem::path& path)
{
    std::ifstream input(path);
    if (!input)
    {
        std::cerr << "Could not open configuration file: " << path << '\n';
        return std::nullopt;
    }
    return std::string(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
}

std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
sample_initial_population(const qmomkit::config::InitialConditionSpecification& initial_condition, const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    return std::visit(
        [&hierarchy](const auto& specification) -> std::optional<qmomkit::monte_carlo::TensorSampledInitialization>
        {
            using Specification = std::decay_t<decltype(specification)>;
            if constexpr (std::is_same_v<Specification, qmomkit::config::IndependentProductInitialization>)
            {
                const std::vector<std::size_t> orders = specification.maximum_orders.value_or(maximum_orders(hierarchy.moment_extents()));
                return qmomkit::monte_carlo::sample_initialization(specification.distribution, orders,
                                                                   {initialization_sample_count, initialization_seed});
            }
            std::cerr << "This example requires a sampleable independent-product initial condition.\n";
            return std::nullopt;
        },
        initial_condition);
}

qmomkit::PhaseSpaceVelocity particle_velocity()
{
    return [](double, const std::vector<double>& coordinates)
    {
        if (coordinates.size() != 4 || !(coordinates[0] > 0.0))
        {
            throw std::invalid_argument("Particle coordinates must be ordered as positive (d,u,v,w).");
        }
        const double du = gas_u - coordinates[1];
        const double dv = gas_v - coordinates[2];
        const double dw = gas_w - coordinates[3];
        const double slip = std::sqrt(du * du + dv * dv + dw * dw);
        const double reynolds = gas_density * coordinates[0] * slip / gas_viscosity;
        const double drag = reynolds > 1.0e-14 ? 3.0 * gas_density / (4.0 * particle_density * coordinates[0]) *
                                                     (24.0 / reynolds * (1.0 + 0.15 * std::pow(reynolds, 0.687))) * slip
                                               : 18.0 * gas_viscosity / (particle_density * coordinates[0] * coordinates[0]);
        return std::vector<double>{0.0, drag * du, drag * dv, drag * dw};
    };
}

qmomkit::CqmomMomentSource breakup_source()
{
    return [](double, const qmomkit::CqmomSourceContext& context)
    {
        if (context.dimension_count() != 4)
        {
            throw std::invalid_argument("The spherical recoil source requires phase space (d,u,v,w).");
        }
        std::vector<double> rates(context.size(), 0.0);
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const qmomkit::MultiIndex& order = context.index(entry);
            context.for_each_terminal_node(
                [&](const double weight, const std::vector<double>& coordinates)
                {
                    rates[entry] += weight * breakup_frequency(coordinates[0], coordinates[1], coordinates[2], coordinates[3]) *
                                    (2.0 * expected_daughter_monomial(coordinates, order) - parent_monomial(coordinates, order));
                });
        }
        return qmomkit::CqmomMomentRate{std::move(rates)};
    };
}

void apply_breakup(qmomkit::monte_carlo::Ensemble& ensemble, const double interval, std::mt19937_64& generator)
{
    qmomkit::monte_carlo::Ensemble completed;
    completed.reserve(2 * ensemble.size());
    std::normal_distribution<double> recoil(0.0, recoil_standard_deviation);
    for (const qmomkit::monte_carlo::Sample& sample : ensemble)
    {
        std::vector<std::pair<qmomkit::monte_carlo::Sample, double>> pending{{sample, interval}};
        while (!pending.empty())
        {
            auto [parent, remaining] = std::move(pending.back());
            pending.pop_back();
            if (parent.size() != 4)
            {
                throw std::runtime_error("MC particles must be ordered as (d,u,v,w).");
            }
            std::exponential_distribution<double> wait(breakup_frequency(parent[0], parent[1], parent[2], parent[3]));
            const double event_time = wait(generator);
            if (event_time >= remaining)
            {
                completed.push_back(std::move(parent));
                continue;
            }
            const double d = daughter_diameter_factor * parent[0];
            const double ru = recoil(generator);
            const double rv = recoil(generator);
            const double rw = recoil(generator);
            const double remaining_after_event = remaining - event_time;
            pending.push_back({{d, parent[1] + ru, parent[2] + rv, parent[3] + rw}, remaining_after_event});
            pending.push_back({{d, parent[1] - ru, parent[2] - rv, parent[3] - rw}, remaining_after_event});
        }
    }
    ensemble = std::move(completed);
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

qmomkit::RawMomentRequest diagnostic_request()
{
    return qmomkit::RawMomentRequest{{{0, 0, 0, 0}, {3, 0, 0, 0}, {3, 1, 0, 0}, {3, 0, 1, 0}, {3, 0, 0, 1}, {3, 2, 0, 0}}};
}

Diagnostics diagnostics(const qmomkit::SelectedRawMomentResult& moments)
{
    if (!moments.ok() || moments.values.size() != 6)
    {
        throw std::runtime_error("Could not evaluate spherical-particle diagnostics.");
    }
    return {moments.values[0],
            pi / 6.0 * moments.values[1],
            pi / 6.0 * moments.values[2],
            pi / 6.0 * moments.values[3],
            pi / 6.0 * moments.values[4],
            pi / 6.0 * moments.values[5]};
}

Diagnostics diagnostics(const qmomkit::RawMomentTensor& moments)
{
    return {tensor_value(moments, {0, 0, 0, 0}),
            pi / 6.0 * tensor_value(moments, {3, 0, 0, 0}),
            pi / 6.0 * tensor_value(moments, {3, 1, 0, 0}),
            pi / 6.0 * tensor_value(moments, {3, 0, 1, 0}),
            pi / 6.0 * tensor_value(moments, {3, 0, 0, 1}),
            pi / 6.0 * tensor_value(moments, {3, 2, 0, 0})};
}

void print_diagnostics(const char* label, const Diagnostics& values)
{
    std::cout << label << " N=" << values.number_density << " volume=" << values.volume << " volume_momentum=(" << values.momentum_u << ','
              << values.momentum_v << ',' << values.momentum_w << ") volume_u2=" << values.velocity_square_u << '\n';
}

void ensure_output_directory(const std::filesystem::path& path)
{
    std::filesystem::create_directories(path);
}

} // namespace cqmom_spherical_breakup_recoil
