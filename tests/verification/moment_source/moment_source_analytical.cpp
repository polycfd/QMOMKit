// Moment source terms integrated against solutions that can be written down.
//
// What is verified
//   That advancing moments under a source term reproduces the exact solution of
//   the corresponding ordinary differential equation, for sources that act on
//   each moment order independently and for sources that couple orders.
//
// How the reference is constructed
//   The source terms are chosen so the resulting equations are solvable on
//   paper, and the exact moments are evaluated directly from those solutions.
//   The order-dependent case relaxes each moment towards a target at its own
//   rate, dM_k/dt = -lambda_k*(M_k - T_k), whose solution is
//
//       M_k(t) = T_k + (M_k(0) - T_k)*exp(-lambda_k*t),
//
//   and the coupled case uses a linear system with a known closed form. The
//   starting moments come from a chosen set of weighted points, so the initial
//   state is exact as well.
//
// Error metric and pass criteria
//   Maximum deviation of the advanced moments from the closed-form solution, at
//   2.0e-6, and of moments rebuilt from the resulting quadrature at 2.0e-8.
//   Both sit above machine precision because the moments are produced by a time
//   integrator over a finite interval: the bound reflects the accuracy the
//   integrator was asked for, not the accuracy of the arithmetic.
//
// Reproduce
//   ctest -R verification_.*moment_source
//
// Build modes
//   Both. Failures are reported through an exit code, not assert().

#include "qmomkit/configuration.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/moment_source.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace
{

constexpr double moment_tolerance = 2.0e-6;
constexpr double reconstruction_tolerance = 2.0e-8;

struct TimePlan
{
    double t0 = 0.0;
    double dt = 0.0;
    std::size_t steps = 0;
};

struct ErrorMetrics
{
    double moment_abs = 0.0;
    double moment_rel = 0.0;
    double reconstruction_abs = 0.0;
    double reconstruction_rel = 0.0;
    double conservation_abs = 0.0;
};

struct VerificationResult
{
    bool passed = true;
    std::string message;
    ErrorMetrics errors;
    double coarse_error = 0.0;
    double refined_error = 0.0;
    std::vector<double> refinement_errors;
    std::vector<double> observed_orders;
};

struct CaseDefinition
{
    std::string name;
    bool cqmom = false;
    TimePlan one_step;
    TimePlan multi_step;
};

void append_failure(VerificationResult& result, const std::string& message)
{
    result.passed = false;
    if (!result.message.empty())
    {
        result.message += '\n';
    }
    result.message += message;
}

void update_error(double& maximum, const double actual, const double expected)
{
    maximum = std::max(maximum, std::abs(actual - expected));
}

void update_relative_error(double& maximum, const double actual, const double expected)
{
    maximum = std::max(maximum, std::abs(actual - expected) / std::max(1.0, std::abs(expected)));
}

double integer_power(const double base, const std::size_t exponent)
{
    double value = 1.0;
    for (std::size_t power = 0; power < exponent; ++power)
    {
        value *= base;
    }
    return value;
}

std::vector<double> quadrature_moments_1d(const std::vector<double>& nodes, const std::vector<double>& weights, const std::size_t count)
{
    std::vector<double> moments(count, 0.0);
    for (std::size_t order = 0; order < count; ++order)
    {
        for (std::size_t node = 0; node < nodes.size(); ++node)
        {
            moments[order] += weights[node] * integer_power(nodes[node], order);
        }
    }
    return moments;
}

qmomkit::MomentSourceUpdateOptions source_options()
{
    qmomkit::MomentSourceUpdateOptions options;
    options.ssp_rk3_options.estimate_error = true;
    return options;
}

const CaseDefinition* find_case(const std::string& name)
{
    static const std::vector<CaseDefinition> cases{
        {"order_dependent_relaxation", false, {0.2, 0.08, 1}, {0.2, 0.04, 6}},
        {"common_rate_relaxation", true, {0.1, 0.06, 1}, {0.1, 0.06, 5}},
        {"linear_coupled", false, {-0.2, 0.08, 1}, {-0.2, 0.06, 8}},
    };
    for (const CaseDefinition& candidate : cases)
    {
        if (candidate.name == name)
        {
            return &candidate;
        }
    }
    return nullptr;
}

std::vector<double> target_order_dependent_moments(const std::size_t count)
{
    return quadrature_moments_1d({-0.4, 0.7, 1.5}, {0.25, 0.35, 0.40}, count);
}

std::vector<double> exact_order_dependent(const std::vector<double>& initial, const double elapsed)
{
    const std::vector<double> target = target_order_dependent_moments(initial.size());
    std::vector<double> exact(initial.size(), 0.0);
    for (std::size_t order = 0; order < initial.size(); ++order)
    {
        const double lambda = order == 0 ? 0.0 : 0.15 * static_cast<double>(order);
        exact[order] = target[order] + (initial[order] - target[order]) * std::exp(-lambda * elapsed);
    }
    return exact;
}

std::vector<double> exact_linear_coupled(const std::vector<double>& initial, const double elapsed)
{
    constexpr double a = -0.2;
    constexpr double b = 0.15;
    constexpr double c = 0.1;
    constexpr double d = 0.6;
    const double exp_a = std::exp(a * elapsed);
    const double exp_c = std::exp(c * elapsed);
    const double exp_d = std::exp(d * elapsed);
    return {initial[0], initial[1] * exp_a, initial[2] * exp_c + b * initial[1] * (exp_a - exp_c) / (a - c), initial[3] * exp_d};
}

std::vector<double> exact_qmom_moments(const std::string& name, const std::vector<double>& initial, const double elapsed)
{
    if (name == "order_dependent_relaxation")
    {
        return exact_order_dependent(initial, elapsed);
    }
    if (name == "linear_coupled")
    {
        return exact_linear_coupled(initial, elapsed);
    }
    return {};
}

std::vector<double> qmom_source(const std::string& name, const std::vector<double>& moments)
{
    if (name == "order_dependent_relaxation")
    {
        const std::vector<double> target = target_order_dependent_moments(moments.size());
        std::vector<double> derivative(moments.size(), 0.0);
        for (std::size_t order = 1; order < moments.size(); ++order)
        {
            derivative[order] = -0.15 * static_cast<double>(order) * (moments[order] - target[order]);
        }
        return derivative;
    }
    if (name == "linear_coupled")
    {
        return {0.0, -0.2 * moments[1], 0.15 * moments[1] + 0.1 * moments[2], 0.6 * moments[3]};
    }
    return {};
}

struct TerminalNode
{
    double x = 0.0;
    double y = 0.0;
    double weight = 0.0;
};

std::vector<TerminalNode> target_common_rate_nodes()
{
    return {{-0.5, -1.0, 0.25}, {-0.5, 0.8, 0.20}, {1.2, -1.0, 0.20}, {1.2, 0.8, 0.35}};
}

qmomkit::RawMomentTensor moments_from_terminal_nodes(const std::vector<TerminalNode>& nodes, const std::vector<std::size_t>& extents)
{
    qmomkit::RawMomentTensor moments{extents, std::vector<double>(qmomkit::product_extents(extents), 0.0)};
    for (std::size_t flat = 0; flat < moments.values.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, extents);
        for (const TerminalNode& node : nodes)
        {
            moments.values[flat] += node.weight * integer_power(node.x, index[0]) * integer_power(node.y, index[1]);
        }
    }
    return moments;
}

qmomkit::RawMomentTensor target_common_rate_tensor(const std::vector<std::size_t>& extents)
{
    return moments_from_terminal_nodes(target_common_rate_nodes(), extents);
}

bool is_zeroth_index(const qmomkit::MultiIndex& index)
{
    return std::all_of(index.begin(), index.end(), [](const std::size_t value) { return value == 0; });
}

std::vector<double> exact_common_rate(const qmomkit::CqmomMomentData& initial, const double elapsed)
{
    constexpr double lambda = 0.45;
    const qmomkit::CqmomMomentData target =
        qmomkit::extract_cqmom_moments(target_common_rate_tensor(initial.layout.hierarchy.moment_extents()), initial.layout);
    std::vector<double> exact(initial.values.size(), 0.0);
    const double factor = std::exp(-lambda * elapsed);
    const std::vector<qmomkit::MultiIndex> indices = initial.layout.required_indices.to_multi_indices();
    for (std::size_t entry = 0; entry < initial.values.size(); ++entry)
    {
        exact[entry] =
            is_zeroth_index(indices[entry]) ? initial.values[entry] : target.values[entry] + (initial.values[entry] - target.values[entry]) * factor;
    }
    return exact;
}

std::vector<double> common_rate_source(const qmomkit::CqmomMomentData& initial_shape, const std::vector<double>& moments)
{
    constexpr double lambda = 0.45;
    const qmomkit::CqmomMomentData target =
        qmomkit::extract_cqmom_moments(target_common_rate_tensor(initial_shape.layout.hierarchy.moment_extents()), initial_shape.layout);
    std::vector<double> derivative(moments.size(), 0.0);
    const std::vector<qmomkit::MultiIndex> indices = initial_shape.layout.required_indices.to_multi_indices();
    for (std::size_t entry = 0; entry < moments.size(); ++entry)
    {
        if (!is_zeroth_index(indices[entry]))
        {
            derivative[entry] = -lambda * (moments[entry] - target.values[entry]);
        }
    }
    return derivative;
}

void compare_vectors(const std::vector<double>& actual, const std::vector<double>& expected, ErrorMetrics& errors)
{
    for (std::size_t index = 0; index < expected.size(); ++index)
    {
        update_error(errors.moment_abs, actual[index], expected[index]);
        update_relative_error(errors.moment_rel, actual[index], expected[index]);
    }
}

void compare_reconstruction(const std::vector<double>& reconstructed, const std::vector<double>& expected, ErrorMetrics& errors)
{
    for (std::size_t index = 0; index < expected.size(); ++index)
    {
        update_error(errors.reconstruction_abs, reconstructed[index], expected[index]);
        update_relative_error(errors.reconstruction_rel, reconstructed[index], expected[index]);
    }
}

void verify_source_diagnostics(const qmomkit::QmomMomentSourceUpdateResult& update, const double expected_dt, VerificationResult& result)
{
    if (!update.ok())
    {
        append_failure(result, "Moment-source update failed: " + update.message);
        return;
    }
    if (update.diagnostics.ssp_rk3_status != qmomkit::SspRk3Status::success || !update.diagnostics.accepted)
    {
        append_failure(result, "SSP-RK3 diagnostics did not report a successful accepted source step.");
    }
    if (update.diagnostics.rejected || update.diagnostics.retries != 0)
    {
        append_failure(result, "Unexpected SSP-RK3 retry occurred in an analytical source verification case.");
    }
    if (std::abs(update.diagnostics.time_step_used - expected_dt) > reconstruction_tolerance)
    {
        append_failure(result, "SSP-RK3 diagnostics reported the wrong accepted timestep.");
    }
}

std::vector<double> advance_qmom_plan(const CaseDefinition& definition, const std::vector<double>& initial, const TimePlan& plan,
                                      VerificationResult& result)
{
    std::vector<double> current = initial;
    double time = plan.t0;
    for (std::size_t step = 0; step < plan.steps; ++step)
    {
        const qmomkit::QmomMomentSourceUpdateResult update = qmomkit::advance_qmom_moments_with_source(
            current, time, plan.dt, [&definition](double, const qmomkit::QmomSourceContext& context)
            { return qmomkit::QmomMomentRate{qmom_source(definition.name, context.values())}; }, source_options());
        verify_source_diagnostics(update, plan.dt, result);
        if (!update.ok())
        {
            return current;
        }
        current = update.moments;
        time += plan.dt;
        update_error(result.errors.conservation_abs, current.front(), initial.front());
    }
    return current;
}

std::vector<double> advance_cqmom_plan(const qmomkit::CqmomMomentData& initial_moments, const TimePlan& plan, VerificationResult& result)
{
    std::vector<double> current = initial_moments.values;
    double time = plan.t0;
    for (std::size_t step = 0; step < plan.steps; ++step)
    {
        const qmomkit::QmomMomentSourceUpdateResult update = qmomkit::advance_qmom_moments_with_source(
            current, time, plan.dt, [&initial_moments](double, const qmomkit::QmomSourceContext& context)
            { return qmomkit::QmomMomentRate{common_rate_source(initial_moments, context.values())}; }, source_options());
        verify_source_diagnostics(update, plan.dt, result);
        if (!update.ok())
        {
            return current;
        }
        current = update.moments;
        time += plan.dt;
        update_error(result.errors.conservation_abs, current.front(), initial_moments.values.front());
    }
    return current;
}

void verify_qmom_reconstruction(const qmomkit::QmomCase& configured_qmom_case, const std::vector<double>& moments, VerificationResult& result)
{
    qmomkit::QmomCase final_qmom_case = configured_qmom_case;
    final_qmom_case.moments = qmomkit::MomentSet(moments);
    const qmomkit::QmomResult qmom = qmomkit::compute_qmom(final_qmom_case);
    if (!qmom.ok())
    {
        append_failure(result, "Final QMOM inversion failed: " + qmom.message);
        return;
    }
    const qmomkit::Quadrature1D quadrature{qmom.abscissas, qmom.weights, qmomkit::QuadratureStatus::success, qmom.message};
    compare_reconstruction(qmomkit::reconstruct_raw_moments(quadrature, moments.size()), moments, result.errors);
}

VerificationResult verify_qmom_case(const CaseDefinition& definition, const std::filesystem::path& configuration_path)
{
    VerificationResult result;
    const auto configured = qmomkit::config::load_qmom(configuration_path);
    if (!configured.ok())
    {
        append_failure(result, "Failed to load QMOM configuration.");
        return result;
    }

    const std::vector<double> initial = configured.value().moments.values();
    const qmomkit::QmomResult initial_qmom = qmomkit::compute_qmom(configured.value());
    if (!initial_qmom.ok())
    {
        append_failure(result, "Initial QMOM inversion failed: " + initial_qmom.message);
        return result;
    }

    for (const TimePlan& plan : {definition.one_step, definition.multi_step})
    {
        const std::vector<double> evolved = advance_qmom_plan(definition, initial, plan, result);
        const std::vector<double> exact = exact_qmom_moments(definition.name, initial, plan.dt * static_cast<double>(plan.steps));
        compare_vectors(evolved, exact, result.errors);
        update_error(result.errors.conservation_abs, evolved.front(), initial.front());
        verify_qmom_reconstruction(configured.value(), evolved, result);
        verify_qmom_reconstruction(configured.value(), exact, result);
    }
    return result;
}

void verify_cqmom_reconstruction(const qmomkit::CqmomCase& configured, const std::vector<double>& moments, VerificationResult& result)
{
    qmomkit::CqmomCase final_case = configured;
    final_case.moments.values = moments;
    const qmomkit::CqmomResult cqmom = qmomkit::compute_cqmom(final_case);
    if (!cqmom.ok())
    {
        append_failure(result, "Final CQMOM inversion failed: " + cqmom.message);
        return;
    }
    if (cqmom.quadrature.hierarchy().conditioning_order() != configured.moments.layout.hierarchy.conditioning_order())
    {
        append_failure(result, "CQMOM conditioning order was not preserved.");
    }
    if (cqmom.quadrature.adaptive_branch_counts_enabled() != configured.adaptive_branch_counts || !cqmom.quadrature.has_rectangular_storage())
    {
        append_failure(result, "CQMOM fixed/adaptive storage settings were not preserved.");
    }
    const qmomkit::CqmomMomentData reconstructed = qmomkit::reconstruct_cqmom_moments(cqmom.quadrature, configured.moments.layout);
    compare_reconstruction(reconstructed.values, moments, result.errors);
}

VerificationResult verify_cqmom_case(const CaseDefinition& definition, const std::filesystem::path& configuration_path)
{
    VerificationResult result;
    const auto configured = qmomkit::config::load_cqmom(configuration_path);
    if (!configured.ok())
    {
        append_failure(result, "Failed to load CQMOM configuration.");
        return result;
    }
    if (configured.value().moments.layout.hierarchy.conditioning_order() != std::vector<std::size_t>{1, 0} ||
        configured.value().qmom_options.adaptive || configured.value().adaptive_branch_counts)
    {
        append_failure(result, "CQMOM configuration did not preserve the expected non-identity fixed settings.");
    }
    const qmomkit::CqmomResult initial_cqmom = qmomkit::compute_cqmom(configured.value());
    if (!initial_cqmom.ok())
    {
        append_failure(result, "Initial CQMOM inversion failed: " + initial_cqmom.message);
        return result;
    }

    for (const TimePlan& plan : {definition.one_step, definition.multi_step})
    {
        const std::vector<double> evolved = advance_cqmom_plan(configured.value().moments, plan, result);
        const std::vector<double> exact = exact_common_rate(configured.value().moments, plan.dt * static_cast<double>(plan.steps));
        compare_vectors(evolved, exact, result.errors);
        update_error(result.errors.conservation_abs, evolved.front(), configured.value().moments.values.front());
        verify_cqmom_reconstruction(configured.value(), evolved, result);
        verify_cqmom_reconstruction(configured.value(), exact, result);
    }
    return result;
}

double maximum_abs_error(const std::vector<double>& actual, const std::vector<double>& expected)
{
    double maximum = 0.0;
    for (std::size_t index = 0; index < actual.size(); ++index)
    {
        update_error(maximum, actual[index], expected[index]);
    }
    return maximum;
}

void verify_linear_coupled_convergence(const std::filesystem::path& configuration_path, VerificationResult& result)
{
    const auto configured = qmomkit::config::load_qmom(configuration_path);
    if (!configured.ok())
    {
        append_failure(result, "Failed to load QMOM configuration for convergence check.");
        return;
    }
    const std::vector<double> initial = configured.value().moments.values();
    const std::vector<TimePlan> refinement_plans{{0.0, 0.12, 4}, {0.0, 0.06, 8}, {0.0, 0.03, 16}, {0.0, 0.015, 32}};
    const std::vector<double> exact = exact_linear_coupled(initial, 0.48);
    const CaseDefinition* const definition = find_case("linear_coupled");
    for (const TimePlan& plan : refinement_plans)
    {
        result.refinement_errors.push_back(maximum_abs_error(advance_qmom_plan(*definition, initial, plan, result), exact));
    }
    result.coarse_error = result.refinement_errors.front();
    result.refined_error = result.refinement_errors[1];

    for (std::size_t level = 1; level < result.refinement_errors.size(); ++level)
    {
        const double coarse_error = result.refinement_errors[level - 1];
        const double refined_error = result.refinement_errors[level];
        if (!(refined_error < coarse_error))
        {
            append_failure(result, "Refined SSP-RK3 source step did not reduce the linear-coupled final moment error.");
            continue;
        }
        result.observed_orders.push_back(std::log(coarse_error / refined_error) / std::log(2.0));
    }

    if (result.observed_orders.size() == refinement_plans.size() - 1 &&
        (result.observed_orders[result.observed_orders.size() - 2] < 2.8 || result.observed_orders.back() < 2.8))
    {
        append_failure(result, "The two finest linear-coupled SSP-RK3 observed orders must both be at least 2.8.");
    }
}

} // namespace

int main(const int argc, char** argv)
{
    if (argc != 3)
    {
        std::cerr << "Usage: moment_source_analytical <configuration.qmom> <case_identifier>\n";
        return 2;
    }
    const std::filesystem::path configuration_path = argv[1];
    const std::string case_identifier = argv[2];
    const CaseDefinition* const definition = find_case(case_identifier);
    if (definition == nullptr)
    {
        std::cerr << "Unknown moment-source analytical case '" << case_identifier << "'.\n";
        return 2;
    }

    VerificationResult result =
        definition->cqmom ? verify_cqmom_case(*definition, configuration_path) : verify_qmom_case(*definition, configuration_path);
    if (case_identifier == "linear_coupled")
    {
        verify_linear_coupled_convergence(configuration_path, result);
    }
    const bool passed = result.passed && result.errors.moment_abs <= moment_tolerance && result.errors.moment_rel <= moment_tolerance &&
                        result.errors.reconstruction_abs <= reconstruction_tolerance &&
                        result.errors.reconstruction_rel <= reconstruction_tolerance && result.errors.conservation_abs <= reconstruction_tolerance;

    // Silent on success. These numbers exist to diagnose a failure, and CTest
    // surfaces stderr through --output-on-failure.
    if (!(passed))
    {
        std::cerr << "QMOMKit moment-source analytical verification\n";
        std::cerr << "Configuration: " << configuration_path << '\n';
        std::cerr << "Moment tolerance: " << std::scientific << moment_tolerance << ", reconstruction tolerance: " << reconstruction_tolerance
                  << "\n\n";
        std::cerr << std::left << std::setw(32) << "case" << std::right << std::setw(14) << "moment abs" << std::setw(14) << "moment rel"
                  << std::setw(14) << "recon abs" << std::setw(14) << "recon rel" << std::setw(14) << "m0 abs" << std::setw(14) << "coarse err"
                  << std::setw(14) << "refined err" << std::setw(9) << "status" << '\n';
        std::cerr << std::string(125, '-') << '\n';
        std::cerr << std::left << std::setw(32) << definition->name << std::right << std::scientific << std::setprecision(4) << std::setw(14)
                  << result.errors.moment_abs << std::setw(14) << result.errors.moment_rel << std::setw(14) << result.errors.reconstruction_abs
                  << std::setw(14) << result.errors.reconstruction_rel << std::setw(14) << result.errors.conservation_abs << std::setw(14)
                  << result.coarse_error << std::setw(14) << result.refined_error << std::setw(9) << (passed ? "PASS" : "FAIL") << '\n';

        if (case_identifier == "linear_coupled")
        {
            std::cerr << "SSP-RK3 linear-coupled refinement errors:";
            for (const double error : result.refinement_errors)
            {
                std::cerr << ' ' << error;
            }
            std::cerr << '\n';
            if (result.observed_orders.size() >= 2)
            {
                std::cerr << "SSP-RK3 linear-coupled finest observed orders: " << result.observed_orders[result.observed_orders.size() - 2] << ' '
                          << result.observed_orders.back() << '\n';
            }
        }

        if (!result.passed)
        {
            std::cerr << result.message << '\n';
        }
    }

    return passed ? 0 : 1;
}
