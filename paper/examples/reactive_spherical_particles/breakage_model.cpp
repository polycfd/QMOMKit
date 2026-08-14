#include "breakage_model.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cqmom_reactive_spherical_particles
{
namespace
{

constexpr std::size_t dimension_count = 5;
// Equal-volume binary fragmentation gives
// d_daughter / d_parent = 2^(-1/3).
constexpr double daughter_diameter_factor = 0.79370052598409973737585281963615413019574666394993;

[[nodiscard]] bool is_nonnegative_finite(const double value)
{
    return std::isfinite(value) && value >= 0.0;
}

void require_nonnegative_finite(const double value, const char* name)
{
    if (!is_nonnegative_finite(value))
    {
        throw std::invalid_argument(std::string("Binary fragmentation parameter '") + name + "' must be supplied as a non-negative finite value.");
    }
}

void require_sample(const qmomkit::monte_carlo::Sample& sample)
{
    if (sample.size() != dimension_count)
    {
        throw std::runtime_error("The 5D MC particle example expects samples ordered as (d, u, v, w, T).");
    }
}

[[nodiscard]] double integer_power(double base, std::size_t exponent)
{
    double value = 1.0;
    for (std::size_t count = 0; count < exponent; ++count)
    {
        value *= base;
    }
    return value;
}

[[nodiscard]] double daughter_diameter(const double parent_diameter)
{
    return daughter_diameter_factor * parent_diameter;
}

[[nodiscard]] double breakage_rate(const double diameter, const BinaryFragmentationParameters& parameters)
{
    if (diameter < parameters.minimum_fragmentation_diameter)
    {
        return 0.0;
    }

    return parameters.fragmentation_rate;
}

[[nodiscard]] double parent_monomial(const std::vector<std::vector<double>>& abscissas, const std::size_t node, const qmomkit::MultiIndex& order)
{
    double value = 1.0;
    for (std::size_t dimension = 0; dimension < dimension_count; ++dimension)
    {
        value *= integer_power(abscissas[dimension][node], order[dimension]);
    }
    return value;
}

[[nodiscard]] double daughter_monomial(const std::vector<std::vector<double>>& abscissas, const std::size_t node, const qmomkit::MultiIndex& order)
{
    const double d = abscissas[0][node];
    const double u = abscissas[1][node];
    const double v = abscissas[2][node];
    const double w = abscissas[3][node];
    const double temperature = abscissas[4][node];
    return integer_power(daughter_diameter(d), order[0]) * integer_power(u, order[1]) * integer_power(v, order[2]) * integer_power(w, order[3]) *
           integer_power(temperature, order[4]);
}

} // namespace

BinaryFragmentationParameters binary_fragmentation_parameters()
{
    return {};
}

void validate_binary_fragmentation_parameters(const BinaryFragmentationParameters& parameters)
{
    require_nonnegative_finite(parameters.fragmentation_rate, "fragmentation_rate");
    require_nonnegative_finite(parameters.minimum_fragmentation_diameter, "minimum_fragmentation_diameter");
}

qmomkit::CqmomMomentSource binary_breakage_source(const qmomkit::CqmomCase& cqmom_case, const BinaryFragmentationParameters& parameters)
{
    validate_binary_fragmentation_parameters(parameters);
    const qmomkit::ConditionalMomentOptions conditional_options = cqmom_case.conditional_options();
    const qmomkit::QmomOptions qmom_options = cqmom_case.qmom_options();
    const bool adaptive_branch_counts = cqmom_case.adaptive_branch_counts();
    const qmomkit::CqmomDiagnosticLevel diagnostic_level = cqmom_case.diagnostic_level();
    const qmomkit::VandermondeSolver vandermonde_solver = cqmom_case.vandermonde_solver();
    const qmomkit::CqmomCoordinateSupport coordinate_support = cqmom_case.coordinate_support();

    return [parameters, conditional_options, qmom_options, adaptive_branch_counts, diagnostic_level, coordinate_support,
            vandermonde_solver](double, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> source(moments.size(), 0.0);
        for (std::size_t index = 0; index < moments.size(); ++index)
        {
            if (moments.index(index).size() != dimension_count)
            {
                throw std::invalid_argument("The 5D binary breakage source received a moment index with the wrong rank.");
            }
        }

        // Reconstruct the active quadrature at the current moment-source stage before
        // accumulating terminal-node contributions to the requested moment rates.
        qmomkit::CqmomMomentLayout layout{moments.hierarchy(), moments.indices()};
        qmomkit::CqmomMomentData moment_data{std::move(layout), moments.values()};
        const qmomkit::CqmomLeanCase stage_case{std::move(moment_data), conditional_options, qmom_options,      adaptive_branch_counts,
                                                diagnostic_level,       vandermonde_solver,  coordinate_support};
        const qmomkit::CqmomResult stage_result = qmomkit::compute_cqmom(stage_case);
        if (!stage_result.ok())
        {
            throw std::runtime_error("Binary fragmentation source could not reconstruct the stage CQMOM quadrature: " + stage_result.message);
        }
        const qmomkit::CqmomQuadrature& stage_quadrature = stage_result.quadrature;
        if (stage_quadrature.number_of_dimensions() != dimension_count || stage_quadrature.abscissas().size() != dimension_count)
        {
            throw std::runtime_error("Binary fragmentation source reconstructed a quadrature with the wrong rank.");
        }

        const std::vector<double>& weights = stage_quadrature.weights();
        const std::vector<std::vector<double>>& abscissas = stage_quadrature.abscissas();
        for (std::size_t moment_index = 0; moment_index < moments.size(); ++moment_index)
        {
            const qmomkit::MultiIndex& order = moments.index(moment_index);
            for (std::size_t node = 0; node < weights.size(); ++node)
            {
                const double beta = breakage_rate(abscissas[0][node], parameters);
                if (beta == 0.0)
                {
                    continue;
                }
                // One parent is removed and two statistically identical daughters
                // are created, giving beta(d) * (2D - P) before quadrature weighting.
                source[moment_index] +=
                    weights[node] * beta * (2.0 * daughter_monomial(abscissas, node, order) - parent_monomial(abscissas, node, order));
            }
        }
        return qmomkit::CqmomMomentRate{std::move(source)};
    };
}

void apply_binary_breakage(qmomkit::monte_carlo::Ensemble& ensemble, const double step_size, std::mt19937_64& generator,
                           const BinaryFragmentationParameters& parameters)
{
    validate_binary_fragmentation_parameters(parameters);
    if (!is_nonnegative_finite(step_size))
    {
        throw std::invalid_argument("The MC binary breakage stage requires a non-negative finite step size.");
    }

    qmomkit::monte_carlo::Ensemble survivors;
    survivors.reserve(ensemble.size());
    std::vector<std::pair<qmomkit::monte_carlo::Sample, double>> pending;
    pending.reserve(2 * ensemble.size());
    for (qmomkit::monte_carlo::Sample sample : ensemble)
    {
        require_sample(sample);
        pending.emplace_back(std::move(sample), step_size);
        // Advance one lineage through the remaining moment-source stage time.
        // Daughters re-enter with only the unused interval, so repeated
        // fragmentation is allowed without restarting the stage clock.
        while (!pending.empty())
        {
            qmomkit::monte_carlo::Sample particle = std::move(pending.back().first);
            const double remaining_time = pending.back().second;
            pending.pop_back();

            const double beta = breakage_rate(particle[0], parameters);
            if (beta == 0.0 || remaining_time <= 0.0)
            {
                survivors.push_back(std::move(particle));
                continue;
            }

            std::exponential_distribution<double> waiting_time(beta);
            const double event_time = waiting_time(generator);
            if (event_time >= remaining_time)
            {
                survivors.push_back(std::move(particle));
                continue;
            }

            const double child_remaining_time = remaining_time - event_time;
            const double daughter_d = daughter_diameter(particle[0]);
            qmomkit::monte_carlo::Sample daughter{daughter_d, particle[1], particle[2], particle[3], particle[4]};
            qmomkit::monte_carlo::Sample second_daughter = daughter;
            pending.emplace_back(std::move(daughter), child_remaining_time);
            pending.emplace_back(std::move(second_daughter), child_remaining_time);
        }
    }

    ensemble = std::move(survivors);
}

} // namespace cqmom_reactive_spherical_particles
