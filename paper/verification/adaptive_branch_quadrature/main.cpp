#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

constexpr double initial_time = 0.0;
constexpr double final_time = 8.0;
constexpr double time_step = 0.002;
constexpr std::size_t velocity_dimension = 1;
constexpr double relaxation_rate = 5.0;
constexpr double equilibrium_velocity = 0.0;
constexpr double pulse_start = 4.0;
constexpr double pulse_end = 4.05;
constexpr double source_rate = 20.0;
constexpr double velocity_displacement = 0.10;
constexpr double mass_tolerance = 1.0e-8;
constexpr double mean_tolerance = 5.0e-5;

struct Sample
{
    double time = 0.0;
    double mass = 0.0;
    double mean_velocity = 0.0;
    double velocity_variance = 0.0;
    std::vector<std::size_t> branch_counts;
};

std::size_t validate_case(const qmomkit::CqmomCase& case_)
{
    const qmomkit::PhaseSpaceHierarchy& hierarchy = case_.hierarchy();
    if (hierarchy.number_of_dimensions() != 2 || hierarchy.conditioning_order() != std::vector<std::size_t>{0, 1} || !case_.adaptive_branch_counts())
    {
        throw std::runtime_error("The adaptive paper case requires a 2D (d,u) hierarchy, conditioning order {0,1}, and adaptive branch counts.");
    }
    return hierarchy.node_count(velocity_dimension);
}

bool source_pulse_active(const double time)
{
    return time + 1.0e-12 >= pulse_start && time < pulse_end - 1.0e-12;
}

double binomial_coefficient(const std::size_t n, const std::size_t k)
{
    const std::size_t reduced = std::min(k, n - k);
    double coefficient = 1.0;
    for (std::size_t factor = 1; factor <= reduced; ++factor)
    {
        coefficient *= static_cast<double>(n + 1 - factor);
        coefficient /= static_cast<double>(factor);
    }
    return coefficient;
}

qmomkit::PhaseSpaceVelocity velocity()
{
    return qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>& coordinates)
        {
            if (coordinates.size() != 2)
            {
                throw std::invalid_argument("The adaptive paper study requires phase-space coordinates ordered as (d,u).");
            }
            return std::vector<double>{0.0, -relaxation_rate * (coordinates[velocity_dimension] - equilibrium_velocity)};
        });
}

// During the pulse, each velocity distribution is symmetrically displaced by
// +/- velocity_displacement.  This is the case-specific reinjection reference.
qmomkit::CqmomMomentSource variance_injection_source()
{
    return [](const double time, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> derivative(context.size(), 0.0);
        if (!source_pulse_active(time))
        {
            return qmomkit::CqmomMomentRate{std::move(derivative)};
        }

        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const qmomkit::MultiIndexView index = context.index(entry);
            if (index.size() != 2)
            {
                throw std::invalid_argument("The adaptive paper source requires 2D raw moments ordered as (d,u).");
            }

            double redistribution_change = 0.0;
            for (std::size_t displacement_order = 2; displacement_order <= index[velocity_dimension]; displacement_order += 2)
            {
                qmomkit::MultiIndex lower = index.to_multi_index();
                lower[velocity_dimension] -= displacement_order;
                const std::optional<double> lower_moment = context.value(lower);
                if (!lower_moment.has_value())
                {
                    throw std::invalid_argument("The variance-injection source requires a lower-order velocity moment.");
                }
                redistribution_change += binomial_coefficient(index[velocity_dimension], displacement_order) *
                                         std::pow(velocity_displacement, static_cast<double>(displacement_order)) * *lower_moment;
            }
            derivative[entry] = source_rate * redistribution_change;
        }
        return qmomkit::CqmomMomentRate{std::move(derivative)};
    };
}

bool finite_quadrature(const qmomkit::CqmomQuadrature& quadrature)
{
    const auto finite_vector = [](const std::vector<double>& values)
    { return std::all_of(values.begin(), values.end(), [](const double value) { return std::isfinite(value); }); };
    return finite_vector(quadrature.weights()) && std::all_of(quadrature.abscissas().begin(), quadrature.abscissas().end(), finite_vector);
}

std::vector<std::size_t> velocity_branch_counts(const std::vector<qmomkit::CqmomActiveBranchCount>& active_branches,
                                                const std::size_t configured_size_node_count, const std::size_t configured_velocity_node_count)
{
    std::vector<std::size_t> counts(configured_size_node_count, 0);
    std::vector<bool> found(configured_size_node_count, false);
    for (const qmomkit::CqmomActiveBranchCount& branch : active_branches)
    {
        if (branch.conditioning_level != 1 || branch.phase_space_dimension != velocity_dimension || branch.parent_path.size() != 1)
        {
            continue;
        }

        const std::size_t branch_index = branch.parent_path.node_index_at_level(0);
        if (branch_index >= configured_size_node_count || branch.configured_node_count != configured_velocity_node_count ||
            branch.active_node_count == 0 || branch.active_node_count > configured_velocity_node_count)
        {
            throw std::runtime_error("The adaptive paper study received an invalid velocity-branch count.");
        }
        counts[branch_index] = branch.active_node_count;
        found[branch_index] = true;
    }
    if (!std::all_of(found.begin(), found.end(), [](const bool value) { return value; }))
    {
        throw std::runtime_error("The adaptive paper study could not find every velocity-conditioned size branch.");
    }
    return counts;
}

// Native quadrature reconstruction supplies the mass and velocity statistics
// reported in the paper figure and used for the analytical checks below.
Sample make_sample(const double time, const qmomkit::CqmomQuadrature& quadrature, const std::vector<qmomkit::CqmomActiveBranchCount>& active_branches,
                   const std::size_t configured_size_node_count, const std::size_t configured_velocity_node_count)
{
    const qmomkit::RawMomentRequest request{{{0, 0}, {0, 1}, {0, 2}}};
    const qmomkit::SelectedRawMomentResult moments = qmomkit::reconstruct_raw_moments(quadrature, request);
    if (!moments.ok() || moments.values.size() != 3 ||
        !std::all_of(moments.values.begin(), moments.values.end(), [](const double value) { return std::isfinite(value); }) ||
        moments.values[0] <= 0.0)
    {
        throw std::runtime_error("The adaptive paper study could not reconstruct finite velocity statistics.");
    }

    const double mean = moments.values[1] / moments.values[0];
    const double variance = moments.values[2] / moments.values[0] - mean * mean;
    if (!std::isfinite(mean) || !std::isfinite(variance))
    {
        throw std::runtime_error("The adaptive paper study produced non-finite velocity statistics.");
    }
    return Sample{time, moments.values[0], mean, variance,
                  velocity_branch_counts(active_branches, configured_size_node_count, configured_velocity_node_count)};
}

void write_csv(const std::vector<Sample>& samples)
{
    std::ofstream output("adaptive_branch_quadrature.csv");
    if (!output)
    {
        throw std::runtime_error("Could not open adaptive_branch_quadrature.csv for writing.");
    }
    output << "time,total_mass,mean_velocity,velocity_variance,branch_index,active_nodes\n";
    output << std::scientific << std::setprecision(17);
    for (const Sample& sample : samples)
    {
        for (std::size_t branch = 0; branch < sample.branch_counts.size(); ++branch)
        {
            output << sample.time << ',' << sample.mass << ',' << sample.mean_velocity << ',' << sample.velocity_variance << ',' << branch << ','
                   << sample.branch_counts[branch] << '\n';
        }
    }
}

void print_summary(const std::vector<std::size_t>& initial_counts, const std::vector<std::size_t>& minimum_counts,
                   const std::vector<std::size_t>& maximum_after_reinjection, const double maximum_mass_error, const double maximum_mean_error)
{
    const auto print_counts = [](const char* label, const std::vector<std::size_t>& counts)
    {
        std::cout << label;
        for (const std::size_t count : counts)
        {
            std::cout << ' ' << count;
        }
        std::cout << '\n';
    };
    std::cout << "Adaptive branch-specific quadrature\n";
    print_counts("initial active counts:", initial_counts);
    print_counts("minimum active counts:", minimum_counts);
    print_counts("maximum counts after reinjection:", maximum_after_reinjection);
    std::cout << std::scientific << std::setprecision(6) << "maximum mass error: " << maximum_mass_error << '\n';
    std::cout << "maximum mean-velocity error: " << maximum_mean_error << '\n';
}

} // namespace

int main(const int argc, char** argv)
{
    if (argc > 2)
    {
        std::cerr << "Usage: qmomkit_paper_adaptive_branch_quadrature [configuration.qmom]\n";
        return 2;
    }

    try
    {
        const std::filesystem::path configuration_file = argc == 2 ? argv[1] : "configuration.qmom";
        auto configured = qmomkit::config::load_cqmom(configuration_file);
        if (!configured.ok())
        {
            std::cerr << configured.error() << '\n';
            return 1;
        }
        auto evolution_options_configured = qmomkit::config::load_cqmom_phase_space_evolution_options(configuration_file);
        if (!evolution_options_configured.ok())
        {
            std::cerr << evolution_options_configured.error() << '\n';
            return 1;
        }
        auto output_options_configured = qmomkit::config::load_output_options(configuration_file);
        if (!output_options_configured.ok())
        {
            std::cerr << output_options_configured.error() << '\n';
            return 1;
        }

        const qmomkit::CqmomCase case_ = std::move(configured).value();
        const std::size_t configured_velocity_node_count = validate_case(case_);
        const std::size_t configured_size_node_count = case_.hierarchy().node_count(0);
        const qmomkit::CqmomPhaseSpaceEvolutionOptions evolution_options = std::move(evolution_options_configured).value();
        const qmomkit::OutputOptions output_options = std::move(output_options_configured).value();

        const qmomkit::CqmomResult initial_reconstruction = qmomkit::compute_cqmom(case_);
        if (!initial_reconstruction.ok() || !finite_quadrature(initial_reconstruction.quadrature))
        {
            std::cerr << "Initial CQMOM reconstruction failed: " << initial_reconstruction.message << '\n';
            return 1;
        }

        qmomkit::CqmomState state = case_.initial_state();
        qmomkit::CqmomQuadrature quadrature = initial_reconstruction.quadrature;
        const Sample initial_sample =
            make_sample(initial_time, quadrature, state.active_branch_counts(), configured_size_node_count, configured_velocity_node_count);
        std::vector<Sample> output_samples{initial_sample};
        const std::vector<std::size_t> initial_counts = initial_sample.branch_counts;
        const double initial_mass = initial_sample.mass;
        const double initial_mean_velocity = initial_sample.mean_velocity;
        std::vector<std::size_t> minimum_counts = initial_counts;
        std::vector<std::size_t> maximum_after_reinjection(configured_size_node_count, 0);
        double maximum_mass_error = 0.0;
        double maximum_mean_error = 0.0;

        const qmomkit::PhaseSpaceVelocity advection_velocity = velocity();
        const qmomkit::CqmomMomentSource reinjection_source = variance_injection_source();
        const double output_interval = output_options.interval.value_or(final_time - initial_time);
        double time = initial_time;
        double next_output_time = initial_time + output_interval;
        while (time < final_time - 1.0e-12)
        {
            const double step = std::min(time_step, final_time - time);
            const qmomkit::CqmomPhaseSpaceEvolutionResult advanced =
                qmomkit::advance_cqmom_phase_space(state, time, step, advection_velocity, reinjection_source, evolution_options);
            if (!advanced.ok() || !advanced.quadrature.has_value() || !finite_quadrature(*advanced.quadrature) ||
                !std::all_of(advanced.state.moment_values().begin(), advanced.state.moment_values().end(),
                             [](const double value) { return std::isfinite(value); }))
            {
                std::cerr << "CQMOM phase-space evolution failed: " << advanced.message << '\n';
                return 1;
            }

            state = advanced.state;
            quadrature = *advanced.quadrature;
            time = advanced.final_time;
            const Sample current =
                make_sample(time, quadrature, state.active_branch_counts(), configured_size_node_count, configured_velocity_node_count);
            for (std::size_t branch = 0; branch < current.branch_counts.size(); ++branch)
            {
                minimum_counts[branch] = std::min(minimum_counts[branch], current.branch_counts[branch]);
                if (time >= pulse_start - 1.0e-12)
                {
                    maximum_after_reinjection[branch] = std::max(maximum_after_reinjection[branch], current.branch_counts[branch]);
                }
            }
            maximum_mass_error = std::max(maximum_mass_error, std::abs(current.mass - initial_mass));
            const double exact_mean =
                equilibrium_velocity + (initial_mean_velocity - equilibrium_velocity) * std::exp(-relaxation_rate * current.time);
            maximum_mean_error = std::max(maximum_mean_error, std::abs(current.mean_velocity - exact_mean));

            if (time + 1.0e-12 >= next_output_time || time >= final_time - 1.0e-12)
            {
                output_samples.push_back(current);
                while (next_output_time <= time + 1.0e-12)
                {
                    next_output_time += output_interval;
                }
            }
        }

        write_csv(output_samples);
        print_summary(initial_counts, minimum_counts, maximum_after_reinjection, maximum_mass_error, maximum_mean_error);

        bool collapsed = true;
        bool regrew = true;
        for (std::size_t branch = 0; branch < configured_size_node_count; ++branch)
        {
            collapsed = collapsed && minimum_counts[branch] < initial_counts[branch];
            regrew = regrew && maximum_after_reinjection[branch] > minimum_counts[branch];
        }
        if (!collapsed || !regrew || maximum_mass_error > mass_tolerance || maximum_mean_error > mean_tolerance)
        {
            std::cerr << "The adaptive branch-specific lifecycle verification failed.\n";
            return 1;
        }
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
