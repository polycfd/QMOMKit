#include "../verification/verification_helpers.hpp"
#include "internal/cqmom_active_counts.hpp"
#include "internal/cqmom_case.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

// Historical component-regression snapshot. Its original fixture-generation procedure is not recoverable from repository history, so it is
// not independent numerical verification. The semantic contract is successful adaptive ragged reinversion and preservation of resulting active
// branch counts through zero-rate state evolution; the physical half-step below is diagnostic context only. Update this fixture only for an
// intentional, independently justified change to that contract or reconstruction policy, never by regenerating current output to make it pass.
constexpr double pi = 3.141592653589793238462643383279502884;

struct Recurrence2Diagnostic
{
    double sigma22 = 0.0;
    double sigma23 = 0.0;
    double a0 = 0.0;
    double a1 = 0.0;
    double b1 = 0.0;
    double variance = 0.0;
};

Recurrence2Diagnostic recurrence2(const std::vector<double>& moments)
{
    Recurrence2Diagnostic recurrence;
    recurrence.a0 = moments[1] / moments[0];
    recurrence.sigma22 = moments[2] - recurrence.a0 * moments[1];
    recurrence.sigma23 = moments[3] - recurrence.a0 * moments[2];
    recurrence.b1 = recurrence.sigma22 / moments[0];
    recurrence.a1 = recurrence.sigma23 / recurrence.sigma22 - moments[1] / moments[0];
    recurrence.variance = moments[2] / moments[0] - (moments[1] / moments[0]) * (moments[1] / moments[0]);
    return recurrence;
}

double oblate_spheroid_shape_factor(const double aspect_ratio, const double reynolds_number)
{
    const double e = std::sqrt(std::abs(1.0 - aspect_ratio * aspect_ratio));
    double k0 = 1.0;
    if (std::abs(1.0 - aspect_ratio) >= 0.01)
    {
        const double one_minus_chi_squared = 1.0 - aspect_ratio * aspect_ratio;
        const double atan_value = aspect_ratio < 1.0 ? std::atan(e / aspect_ratio) : 0.0;
        const double term1 = 8.0 / 3.0 / std::pow(aspect_ratio, 1.0 / 3.0);
        const double term2 = 2.0 * aspect_ratio / one_minus_chi_squared +
                             2.0 * (1.0 - 2.0 * aspect_ratio * aspect_ratio) / std::pow(one_minus_chi_squared, 1.5) * atan_value;
        k0 = term1 / term2;
    }

    return k0 + 0.15 * std::pow(aspect_ratio, 95.91) * std::pow(reynolds_number, 0.687) +
           0.2927 * std::pow(1.0 - aspect_ratio, 0.4374) * std::pow(reynolds_number, 0.7512);
}

qmomkit::PhaseSpaceVelocity physical_velocity()
{
    return [](const double time, const std::vector<double>& state)
    {
        if (state.size() != 3)
        {
            throw std::invalid_argument("The reduced non-spherical particle probe expects coordinates ordered as (d, chi, u).");
        }
        const double d = state[0];
        const double chi = state[1];
        const double u = state[2];
        const double fluid_velocity_x = 2.0 * std::cos(time);
        const double relative_velocity_x = fluid_velocity_x - u;
        const double relative_velocity_magnitude = std::abs(relative_velocity_x);
        const double reynolds_number = 1.2 * relative_velocity_magnitude * d / 1.8e-5;
        const double shape_factor = oblate_spheroid_shape_factor(chi, reynolds_number);
        const double drag_rate = 18.0 * 1.8e-5 * shape_factor / (2000.0 * d * d);

        std::vector<double> rate = state;
        rate[0] = 0.0;
        rate[1] = 0.0;
        rate[2] = drag_rate * relative_velocity_x;
        return rate;
    };
}

qmomkit::CqmomCase exact_case()
{
    const std::string text = R"(
[phase_space]
names = d, chi, u
node_counts = 2, 2, 2
conditioning_order = 0, 1, 2
adaptive_branch_counts = true

[cqmom]
diagnostic_level = detailed

[qmom]
adaptive = true

[vandermonde]
solver = dense

[initial]
type = independent_product
number_density = 1.0
maximum_orders = 4, 4, 4

[initial.dimension.0]
family = normal
mean = 130.0e-6
standard_deviation = 2.0e-5

[initial.dimension.1]
family = beta
alpha = 6.0
beta = 14.0

[initial.dimension.2]
family = normal
mean = 1.0
standard_deviation = 1.0
)";
    const auto parsed = qmomkit::config::parse_cqmom(text, "physical_transport_probe");
    QMOMKIT_VERIFY_CHECK(parsed.ok());
    return parsed.value();
}

std::vector<qmomkit::ConditioningPath> captured_terminal_paths()
{
    return {qmomkit::ConditioningPath{{0, 0, 0}}, qmomkit::ConditioningPath{{0, 1, 0}}, qmomkit::ConditioningPath{{1, 0, 0}},
            qmomkit::ConditioningPath{{1, 1, 0}}};
}

std::vector<qmomkit::ConditioningPath> captured_terminal_parent_paths()
{
    return {qmomkit::ConditioningPath{{0, 0}}, qmomkit::ConditioningPath{{0, 1}}, qmomkit::ConditioningPath{{1, 0}},
            qmomkit::ConditioningPath{{1, 1}}};
}

qmomkit::CqmomQuadrature captured_failing_terminal_quadrature()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {0, 1, 2});
    return qmomkit::internal::CqmomConstructionAccess::create_quadrature(
        hierarchy, {0.29500670652258476, 0.203375676586807, 0.29604427987962095, 0.20557333701099686},
        {{0.00010988241902807785, 0.00010988241902807785, 0.00014945058007081771, 0.00014945058007081771},
         {0.21441293841351547, 0.41880958133916973, 0.21711664207959502, 0.4193843251358394},
         {-1.1044079065934338, -1.1181660600053236, -1.1559642207840419, -1.1762606826705884}},
        true, false, captured_terminal_paths());
}

qmomkit::CqmomQuadrature captured_preceding_success_quadrature()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {0, 1, 2});
    return qmomkit::internal::CqmomConstructionAccess::create_quadrature(
        hierarchy, {0.29500670652258476, 0.203375676586807, 0.29604427987962095, 0.20557333701099686},
        {{0.00010988241902807785, 0.00010988241902807785, 0.00014945058007081771, 0.00014945058007081771},
         {0.21441293841351544, 0.41880958133916973, 0.21711664207959502, 0.41938432513583945},
         {-1.1085764283676269, -1.1223128304573988, -1.1600506567321607, -1.1803132705368153}},
        true, false,
        {qmomkit::ConditioningPath{{0, 0, 0}}, qmomkit::ConditioningPath{{0, 1, 0}}, qmomkit::ConditioningPath{{1, 0, 0}},
         qmomkit::ConditioningPath{{1, 1, 0}}});
}

std::vector<qmomkit::CqmomActiveBranchCount> captured_active_branch_counts()
{
    return {
        {0, qmomkit::ConditioningPath{}, 0, 2, 2},       {1, qmomkit::ConditioningPath{{0}}, 1, 2, 2},
        {1, qmomkit::ConditioningPath{{1}}, 1, 2, 2},    {2, qmomkit::ConditioningPath{{0, 0}}, 2, 2, 1},
        {2, qmomkit::ConditioningPath{{0, 1}}, 2, 2, 1}, {2, qmomkit::ConditioningPath{{1, 0}}, 2, 2, 1},
        {2, qmomkit::ConditioningPath{{1, 1}}, 2, 2, 1},
    };
}

bool path_has_prefix(const qmomkit::ConditioningPath& path, const qmomkit::ConditioningPath& prefix)
{
    if (prefix.size() > path.size())
    {
        return false;
    }
    for (std::size_t level = 0; level < prefix.size(); ++level)
    {
        if (path.node_index_at_level(level) != prefix.node_index_at_level(level))
        {
            return false;
        }
    }
    return true;
}

std::vector<double> direct_branch_moments(const qmomkit::CqmomQuadrature& quadrature, const qmomkit::ConditioningPath& branch_path,
                                          const std::size_t dimension, const std::size_t moment_count, double& branch_mass)
{
    branch_mass = 0.0;
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        if (path_has_prefix(quadrature.terminal_paths()[node], branch_path))
        {
            branch_mass += quadrature.weights()[node];
        }
    }
    QMOMKIT_VERIFY_CHECK(branch_mass > 0.0);
    std::vector<double> moments(moment_count, 0.0);
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        if (path_has_prefix(quadrature.terminal_paths()[node], branch_path))
        {
            const double conditional_weight = quadrature.weights()[node] / branch_mass;
            double power = 1.0;
            for (std::size_t order = 0; order < moment_count; ++order)
            {
                moments[order] += conditional_weight * power;
                power *= quadrature.abscissas()[dimension][node];
            }
        }
    }
    return moments;
}

double direct_long_double_variance(const qmomkit::CqmomQuadrature& quadrature, const qmomkit::ConditioningPath& branch_path,
                                   const std::size_t dimension)
{
    long double branch_mass = 0.0L;
    long double mean = 0.0L;
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        if (path_has_prefix(quadrature.terminal_paths()[node], branch_path))
        {
            branch_mass += static_cast<long double>(quadrature.weights()[node]);
            mean += static_cast<long double>(quadrature.weights()[node]) * static_cast<long double>(quadrature.abscissas()[dimension][node]);
        }
    }
    mean /= branch_mass;
    long double variance = 0.0L;
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        if (path_has_prefix(quadrature.terminal_paths()[node], branch_path))
        {
            const long double conditional_weight = static_cast<long double>(quadrature.weights()[node]) / branch_mass;
            const long double delta = static_cast<long double>(quadrature.abscissas()[dimension][node]) - mean;
            variance += conditional_weight * delta * delta;
        }
    }
    return static_cast<double>(variance);
}

void report_tree_spread(const std::string& label, const qmomkit::CqmomQuadrature& quadrature)
{
    for (std::size_t d_branch = 0; d_branch < 2; ++d_branch)
    {
        double min_d = std::numeric_limits<double>::infinity();
        double max_d = -std::numeric_limits<double>::infinity();
        const qmomkit::ConditioningPath prefix{{d_branch}};
        for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
        {
            if (path_has_prefix(quadrature.terminal_paths()[node], prefix))
            {
                min_d = std::min(min_d, quadrature.abscissas()[0][node]);
                max_d = std::max(max_d, quadrature.abscissas()[0][node]);
            }
        }
        std::cout << label << " d_branch=" << d_branch << " d_spread=" << max_d - min_d << '\n';
    }
}

bool active_count_present(const std::vector<qmomkit::CqmomActiveBranchCount>& counts, const qmomkit::ConditioningPath& path,
                          const std::size_t active_node_count)
{
    for (const qmomkit::CqmomActiveBranchCount& count : counts)
    {
        if (count.parent_path == path && count.active_node_count == active_node_count)
        {
            return true;
        }
    }
    return false;
}

bool has_terminal_u_reduction_diagnostic(const std::vector<qmomkit::CqmomDiagnostics>& diagnostics, const qmomkit::ConditioningPath& path)
{
    for (const qmomkit::CqmomDiagnostics& diagnostic : diagnostics)
    {
        if (diagnostic.conditioning_level == 2 && diagnostic.phase_space_dimension == 2 && diagnostic.parent_path == path &&
            diagnostic.requested_node_count == 2 && diagnostic.actual_node_count == 1 && diagnostic.adaptive_node_count_reduced &&
            diagnostic.qmom_status == qmomkit::QmomStatus::success)
        {
            return true;
        }
    }
    return false;
}

void run_captured_physical_regression()
{
    const qmomkit::CqmomQuadrature captured = captured_failing_terminal_quadrature();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(captured.hierarchy());
    const qmomkit::CqmomMomentData moments = qmomkit::reconstruct_cqmom_moments(captured, layout);

    qmomkit::CqmomLeanCase lean_case{moments};
    lean_case.qmom_options.adaptive = true;
    lean_case.adaptive_branch_counts = true;
    lean_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::detailed;
    lean_case.vandermonde_solver = qmomkit::VandermondeSolver::dense;

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(lean_case);
    std::cout << "[physical captured fixture] status=" << static_cast<int>(result.status) << " message=\"" << result.message << "\"\n";
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.quadrature.total_nodes() == 4);
    QMOMKIT_VERIFY_CHECK(result.quadrature.terminal_paths() == captured_terminal_paths());

    for (const qmomkit::ConditioningPath& path : captured_terminal_parent_paths())
    {
        QMOMKIT_VERIFY_CHECK(has_terminal_u_reduction_diagnostic(result.diagnostics, path));
    }

    const qmomkit::ConditioningPath expected_path{{0, 0}};
    const std::vector<qmomkit::CqmomActiveBranchCount> reconstructed_active_counts =
        qmomkit::internal::active_branch_counts_from_quadrature(result.quadrature);
    for (const qmomkit::ConditioningPath& path : captured_terminal_parent_paths())
    {
        QMOMKIT_VERIFY_CHECK(active_count_present(reconstructed_active_counts, path, 1));
    }

    double branch_mass = 0.0;
    const std::vector<double> direct_moments = direct_branch_moments(captured, expected_path, 2, 4, branch_mass);
    qmomkit::WheelerOptions wheeler_options;
    wheeler_options.adaptive = true;
    const qmomkit::Quadrature1D direct = qmomkit::compute_quadrature(qmomkit::MomentSet(direct_moments), 2, wheeler_options);
    qmomkit::WheelerOptions fixed_wheeler_options = wheeler_options;
    fixed_wheeler_options.adaptive = false;
    const qmomkit::Quadrature1D direct_fixed = qmomkit::compute_quadrature(qmomkit::MomentSet(direct_moments), 2, fixed_wheeler_options);
    QMOMKIT_VERIFY_CHECK(direct.ok());
    QMOMKIT_VERIFY_CHECK(direct.weights.size() == 1);

    const double direct_variance =
        direct_moments[2] / direct_moments[0] - (direct_moments[1] / direct_moments[0]) * (direct_moments[1] / direct_moments[0]);
    const double long_double_variance = direct_long_double_variance(captured, expected_path, 2);
    const Recurrence2Diagnostic direct_recurrence = recurrence2(direct_moments);
    const double u = direct_moments[1] / direct_moments[0];
    const double ulp = std::nextafter(u, std::numeric_limits<double>::infinity()) - u;
    QMOMKIT_VERIFY_CHECK(std::abs(direct_variance) <= ulp);
    QMOMKIT_VERIFY_CHECK(long_double_variance == 0.0);
    std::cout << std::setprecision(17);
    std::cout << "[physical captured branch] mass=" << branch_mass << " direct_m0=" << direct_moments[0] << " direct_m1=" << direct_moments[1]
              << " direct_m2=" << direct_moments[2] << " direct_m3=" << direct_moments[3] << " direct_variance=" << direct_variance
              << " long_double_variance=" << long_double_variance << " ulp=" << ulp << " direct_nodes=" << direct.abscissas.size() << '\n';
    std::cout << "[physical recurrence] direct_sigma22=" << direct_recurrence.sigma22 << " direct_b1=" << direct_recurrence.b1
              << " direct_fixed_status=" << static_cast<int>(direct_fixed.status) << " direct_fixed_message=\"" << direct_fixed.message
              << "\" recursive_status=\"" << result.message << "\"\n";

    report_tree_spread("[physical captured tree]", captured);
    report_tree_spread("[physical preceding tree]", captured_preceding_success_quadrature());

    qmomkit::CqmomState state = qmomkit::internal::CqmomConstructionAccess::create_state(lean_case, captured_active_branch_counts());
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };
    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
        qmomkit::advance_cqmom_phase_space(state, 23.022499999999706, 0.0025, zero_velocity, zero_source, options);
    QMOMKIT_VERIFY_CHECK(phase_space_evolved.ok());
    QMOMKIT_VERIFY_CHECK(phase_space_evolved.quadrature.has_value());
    QMOMKIT_VERIFY_CHECK(phase_space_evolved.quadrature->total_nodes() == 4);
    QMOMKIT_VERIFY_CHECK(active_count_present(phase_space_evolved.state.active_branch_counts(), expected_path, 1));

    const qmomkit::CqmomPhaseSpaceEvolutionResult second =
        qmomkit::advance_cqmom_phase_space(phase_space_evolved.state, phase_space_evolved.final_time, 0.0025, zero_velocity, zero_source, options);
    QMOMKIT_VERIFY_CHECK(second.ok());
    QMOMKIT_VERIFY_CHECK(second.quadrature.has_value());
    QMOMKIT_VERIFY_CHECK(second.quadrature->terminal_paths() == phase_space_evolved.quadrature->terminal_paths());
    QMOMKIT_VERIFY_CHECK(active_count_present(second.state.active_branch_counts(), expected_path, 1));
}

} // namespace

int main()
{
    const qmomkit::CqmomCase cqmom_case = exact_case();
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_case);
    QMOMKIT_VERIFY_CHECK(initial.ok());

    qmomkit::NodeAdvectionOptions node_options;
    node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::CqmomNodeAdvectionResult half_step = qmomkit::advance_nodes(initial.quadrature, 0.0, 0.0025, physical_velocity(), node_options);
    QMOMKIT_VERIFY_CHECK(half_step.ok());

    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(cqmom_case.hierarchy());
    const qmomkit::CqmomMomentData moments = qmomkit::reconstruct_cqmom_moments(half_step.quadrature, layout);
    qmomkit::CqmomLeanCase lean_case{moments};
    lean_case.conditional_options = cqmom_case.conditional_options();
    lean_case.qmom_options = cqmom_case.qmom_options();
    lean_case.adaptive_branch_counts = cqmom_case.adaptive_branch_counts();
    lean_case.diagnostic_level = cqmom_case.diagnostic_level();
    lean_case.vandermonde_solver = cqmom_case.vandermonde_solver();

    const qmomkit::CqmomResult reinverted = qmomkit::compute_cqmom(lean_case);
    std::cout << "[physical exact probe] status=" << static_cast<int>(reinverted.status) << " message=\"" << reinverted.message
              << "\" diagnostics=" << reinverted.diagnostics.size() << '\n';
    for (const qmomkit::CqmomDiagnostics& diagnostic : reinverted.diagnostics)
    {
        if (diagnostic.qmom_status != qmomkit::QmomStatus::success)
        {
            std::cout << "[physical exact probe] failed level=" << diagnostic.conditioning_level << " dim=" << diagnostic.phase_space_dimension
                      << " path=" << diagnostic.parent_path.to_string() << " requested=" << diagnostic.requested_node_count
                      << " actual=" << diagnostic.actual_node_count << " message=\"" << diagnostic.message << "\"\n";
        }
    }
    run_captured_physical_regression();
    return 0;
}
