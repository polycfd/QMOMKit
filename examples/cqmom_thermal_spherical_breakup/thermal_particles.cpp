#include "thermal_particles.hpp"

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

namespace cqmom_thermal_spherical_breakup
{
namespace
{

constexpr double pi = 3.14159265358979323846;
constexpr double gas_density = 1.2;
constexpr double gas_viscosity = 1.8e-5;
constexpr double gas_thermal_conductivity = 0.026;
constexpr double gas_prandtl_number = 0.71;
constexpr double gas_temperature = 300.0;
constexpr double particle_density = 2000.0;
constexpr double particle_heat_capacity = 1000.0;
constexpr double gas_u = 2.0;
constexpr double gas_v = 1.5;
constexpr double gas_w = 0.8;
constexpr double reference_diameter = 1.0e-4;
constexpr double reference_speed = 2.5;
constexpr double breakup_frequency_scale = 0.1;

double diameter_from_volume(const double volume)
{
    if (!(volume > 0.0))
    {
        throw std::invalid_argument("Particle volume must remain positive.");
    }
    return std::cbrt(6.0 * volume / pi);
}

double breakup_frequency(const double volume, const double u, const double v, const double w)
{
    const double d = diameter_from_volume(volume);
    const double du = gas_u - u;
    const double dv = gas_v - v;
    const double dw = gas_w - w;
    return breakup_frequency_scale * d / reference_diameter * (du * du + dv * dv + dw * dw) / (reference_speed * reference_speed);
}

double breakup_volume_coefficient()
{
    return breakup_frequency_scale * std::cbrt(6.0 / pi) / (reference_diameter * reference_speed * reference_speed);
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
        if (coordinates.size() != 5)
        {
            throw std::invalid_argument("Thermal particle coordinates must be ordered as (V,T,u,v,w).");
        }
        const double d = diameter_from_volume(coordinates[0]);
        const double du = gas_u - coordinates[2];
        const double dv = gas_v - coordinates[3];
        const double dw = gas_w - coordinates[4];
        const double slip = std::sqrt(du * du + dv * dv + dw * dw);
        const double reynolds = gas_density * d * slip / gas_viscosity;
        const double drag = reynolds > 1.0e-14 ? 3.0 * gas_density / (4.0 * particle_density * d) *
                                                     (24.0 / reynolds * (1.0 + 0.15 * std::pow(reynolds, 0.687))) * slip
                                               : 18.0 * gas_viscosity / (particle_density * d * d);
        const double nusselt = 2.0 + 0.6 * std::sqrt(reynolds) * std::cbrt(gas_prandtl_number);
        const double temperature_rate =
            6.0 * gas_thermal_conductivity * nusselt / (particle_density * particle_heat_capacity * d * d) * (gas_temperature - coordinates[1]);
        return std::vector<double>{0.0, temperature_rate, drag * du, drag * dv, drag * dw};
    };
}

qmomkit::CqmomMomentSource breakup_source(std::size_t& generalized_request_count)
{
    return [&generalized_request_count](double, const qmomkit::CqmomSourceContext& context)
    {
        if (context.dimension_count() != 5)
        {
            throw std::invalid_argument("The thermal breakup source requires phase space (V,T,u,v,w).");
        }
        std::vector<double> rates(context.size(), 0.0);
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const qmomkit::MultiIndex& order = context.index(entry);
            const double p = static_cast<double>(order[0]);
            const double daughter_factor = std::pow(2.0, 1.0 - p) - 1.0;
            if (daughter_factor == 0.0)
            {
                continue;
            }
            const auto generalized = [&](const std::size_t du, const std::size_t dv, const std::size_t dw)
            {
                ++generalized_request_count;
                return context.generalized_moment({p + 1.0 / 3.0, static_cast<double>(order[1]), static_cast<double>(order[2] + du),
                                                   static_cast<double>(order[3] + dv), static_cast<double>(order[4] + dw)});
            };
            const double scalar_speed_squared = gas_u * gas_u + gas_v * gas_v + gas_w * gas_w;
            const double slip_squared_moment = scalar_speed_squared * generalized(0, 0, 0) - 2.0 * gas_u * generalized(1, 0, 0) -
                                               2.0 * gas_v * generalized(0, 1, 0) - 2.0 * gas_w * generalized(0, 0, 1) + generalized(2, 0, 0) +
                                               generalized(0, 2, 0) + generalized(0, 0, 2);
            rates[entry] = breakup_volume_coefficient() * daughter_factor * slip_squared_moment;
        }
        return qmomkit::CqmomMomentRate{std::move(rates)};
    };
}

void apply_breakup(qmomkit::monte_carlo::Ensemble& ensemble, const double interval, std::mt19937_64& generator)
{
    qmomkit::monte_carlo::Ensemble completed;
    completed.reserve(2 * ensemble.size());
    for (const qmomkit::monte_carlo::Sample& sample : ensemble)
    {
        std::vector<std::pair<qmomkit::monte_carlo::Sample, double>> pending{{sample, interval}};
        while (!pending.empty())
        {
            auto [parent, remaining] = std::move(pending.back());
            pending.pop_back();
            if (parent.size() != 5)
            {
                throw std::runtime_error("MC particles must be ordered as (V,T,u,v,w).");
            }
            std::exponential_distribution<double> wait(breakup_frequency(parent[0], parent[2], parent[3], parent[4]));
            const double event_time = wait(generator);
            if (event_time >= remaining)
            {
                completed.push_back(std::move(parent));
                continue;
            }
            parent[0] *= 0.5;
            const double remaining_after_event = remaining - event_time;
            pending.push_back({parent, remaining_after_event});
            pending.push_back({std::move(parent), remaining_after_event});
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
    return qmomkit::RawMomentRequest{
        {{0, 0, 0, 0, 0}, {1, 0, 0, 0, 0}, {1, 1, 0, 0, 0}, {1, 0, 1, 0, 0}, {1, 0, 0, 1, 0}, {1, 0, 0, 0, 1}, {1, 0, 2, 0, 0}}};
}

Diagnostics diagnostics(const qmomkit::SelectedRawMomentResult& moments)
{
    if (!moments.ok() || moments.values.size() != 7)
    {
        throw std::runtime_error("Could not evaluate thermal spherical-particle diagnostics.");
    }
    return {moments.values[0], moments.values[1], moments.values[2], moments.values[3], moments.values[4], moments.values[5], moments.values[6]};
}

Diagnostics diagnostics(const qmomkit::RawMomentTensor& moments)
{
    return {tensor_value(moments, {0, 0, 0, 0, 0}), tensor_value(moments, {1, 0, 0, 0, 0}), tensor_value(moments, {1, 1, 0, 0, 0}),
            tensor_value(moments, {1, 0, 1, 0, 0}), tensor_value(moments, {1, 0, 0, 1, 0}), tensor_value(moments, {1, 0, 0, 0, 1}),
            tensor_value(moments, {1, 0, 2, 0, 0})};
}

void print_diagnostics(const char* label, const Diagnostics& values)
{
    std::cout << label << " N=" << values.number_density << " volume=" << values.volume << " volume_temperature=" << values.volume_temperature
              << " volume_momentum=(" << values.momentum_u << ',' << values.momentum_v << ',' << values.momentum_w
              << ") volume_u2=" << values.volume_u_square << '\n';
}

void ensure_output_directory(const std::filesystem::path& path)
{
    std::filesystem::create_directories(path);
}

} // namespace cqmom_thermal_spherical_breakup
