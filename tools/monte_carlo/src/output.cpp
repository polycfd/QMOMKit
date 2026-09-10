#include "monte_carlo/output.hpp"

#include "qmomkit/multi_index.hpp"

#include <cstddef>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace qmomkit::monte_carlo
{

void write_moment_text(std::ostream& output, const OneDimensionalMomentEstimate& estimate, const std::size_t sample_count,
                       const unsigned long long seed, const double number_density)
{
    output << std::setprecision(17);
    output << "# samples=" << sample_count << '\n';
    output << "# seed=" << seed << '\n';
    output << "# number_density=" << number_density << '\n';
    output << "# alpha_0 value standard_error\n";
    if (!estimate.ok())
    {
        return;
    }
    for (std::size_t order = 0; order < estimate.moments.size(); ++order)
    {
        output << order << ' ' << estimate.moments[order] << ' ' << estimate.standard_errors[order] << '\n';
    }
}

void write_moment_text(std::ostream& output, const TensorMomentEstimate& estimate, const std::size_t sample_count, const unsigned long long seed,
                       const double number_density)
{
    output << std::setprecision(17);
    output << "# samples=" << sample_count << '\n';
    output << "# seed=" << seed << '\n';
    output << "# number_density=" << number_density << '\n';
    if (estimate.moments.extents.empty())
    {
        output << "# value standard_error\n";
        return;
    }

    output << '#';
    for (std::size_t dimension = 0; dimension < estimate.moments.extents.size(); ++dimension)
    {
        output << " alpha_" << dimension;
    }
    output << " value standard_error\n";
    if (!estimate.ok())
    {
        return;
    }

    for (std::size_t flat = 0; flat < estimate.moments.values.size(); ++flat)
    {
        const qmomkit::MultiIndex order = qmomkit::unflatten_multi_index(flat, estimate.moments.extents);
        for (const std::size_t component : order)
        {
            output << component << ' ';
        }
        output << estimate.moments.values[flat] << ' ' << estimate.standard_errors[flat] << '\n';
    }
}

void write_moment_history_row(qmomkit::MomentWriter& writer, const double time, const qmomkit::RawMomentTensor& moments, const Ensemble& ensemble,
                              const double number_density, const std::vector<std::vector<double>>& generalized_moment_orders)
{
    const GeneralizedEnsembleMomentEstimate generalized = compute_generalized_ensemble_moments(ensemble, number_density, generalized_moment_orders);
    if (!generalized.ok())
    {
        throw std::invalid_argument("Monte Carlo generalized moment output failed: " + generalized.message);
    }
    writer.write(time, moments, generalized.values);
}

void write_moment_history_row(qmomkit::MomentWriter& writer, const double time, const qmomkit::CqmomMomentData& moments, const Ensemble& ensemble,
                              const double number_density, const std::vector<std::vector<double>>& generalized_moment_orders)
{
    const GeneralizedEnsembleMomentEstimate generalized = compute_generalized_ensemble_moments(ensemble, number_density, generalized_moment_orders);
    if (!generalized.ok())
    {
        throw std::invalid_argument("Monte Carlo generalized moment output failed: " + generalized.message);
    }
    writer.write(time, moments, generalized.values);
}

} // namespace qmomkit::monte_carlo
