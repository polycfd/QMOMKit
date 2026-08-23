#include "physics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cqmom_spherical_breakup_recoil
{
namespace
{

constexpr double gas_density = 1.2;
constexpr double gas_viscosity = 1.8e-5;
constexpr double particle_density = 2000.0;
constexpr double gas_u = 2.0;
constexpr double gas_v = 1.5;
constexpr double gas_w = 0.8;
constexpr double breakup_rate = 20.0;
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

} // namespace

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
                    rates[entry] +=
                        weight * breakup_rate * (2.0 * expected_daughter_monomial(coordinates, order) - parent_monomial(coordinates, order));
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
            std::exponential_distribution<double> wait(breakup_rate);
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

} // namespace cqmom_spherical_breakup_recoil
