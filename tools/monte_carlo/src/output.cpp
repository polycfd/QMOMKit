#include "monte_carlo/output.hpp"

#include "qmomkit/multi_index.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace qmomkit::monte_carlo
{

namespace
{
[[nodiscard]] std::vector<double> generalized_values(const Ensemble& ensemble, const double number_density,
                                                     const qmomkit::OutputOptions& output_options)
{
    const GeneralizedEnsembleMomentEstimate generalized =
        compute_generalized_ensemble_moments(ensemble, number_density, output_options.generalized_moment_orders);
    if (!generalized.ok())
    {
        throw std::invalid_argument("Monte Carlo generalized moment output failed: " + generalized.message);
    }
    return generalized.values;
}

} // namespace

qmomkit::MomentWriter make_moment_history_writer(const std::filesystem::path& path, const std::size_t moment_count,
                                                 const qmomkit::OutputOptions& output_options)
{
    return qmomkit::MomentWriter::for_qmom(path, moment_count, output_options);
}

qmomkit::MomentWriter make_moment_history_writer(const std::filesystem::path& path, std::vector<std::size_t> extents,
                                                 const qmomkit::OutputOptions& output_options)
{
    const std::vector<qmomkit::MultiIndex> indices = qmomkit::enumerate_multi_indices(extents);
    return qmomkit::MomentWriter::for_cqmom(path, std::move(extents), indices, output_options.generalized_moment_orders);
}

qmomkit::MomentWriter make_moment_history_writer(const std::filesystem::path& path, std::vector<std::size_t> extents,
                                                 std::vector<qmomkit::MultiIndex> selected_moments, const qmomkit::OutputOptions& output_options)
{
    return qmomkit::MomentWriter::for_cqmom(path, std::move(extents), std::move(selected_moments), output_options.generalized_moment_orders);
}

qmomkit::MomentWriter make_moment_history_writer(const std::filesystem::path& path, const qmomkit::CqmomMomentLayout& layout,
                                                 const qmomkit::OutputOptions& output_options)
{
    return qmomkit::MomentWriter::for_cqmom(path, layout, output_options.generalized_moment_orders);
}

void write_moment_history_row(qmomkit::MomentWriter& writer, const double time, const qmomkit::MomentSet& moments, const Ensemble& ensemble,
                              const double number_density, const qmomkit::OutputOptions& output_options)
{
    if (output_options.generalized_moment_orders.empty())
    {
        writer.write(time, moments);
        return;
    }
    writer.write(time, moments, generalized_values(ensemble, number_density, output_options));
}

void write_moment_history_row(qmomkit::MomentWriter& writer, const double time, const qmomkit::RawMomentTensor& moments, const Ensemble& ensemble,
                              const double number_density, const qmomkit::OutputOptions& output_options)
{
    if (output_options.generalized_moment_orders.empty())
    {
        writer.write(time, moments);
        return;
    }
    writer.write(time, moments, generalized_values(ensemble, number_density, output_options));
}

void write_moment_history_row(qmomkit::MomentWriter& writer, const double time, const qmomkit::CqmomMomentData& moments, const Ensemble& ensemble,
                              const double number_density, const qmomkit::OutputOptions& output_options)
{
    if (output_options.generalized_moment_orders.empty())
    {
        writer.write(time, moments);
        return;
    }
    writer.write(time, moments, generalized_values(ensemble, number_density, output_options));
}

} // namespace qmomkit::monte_carlo
