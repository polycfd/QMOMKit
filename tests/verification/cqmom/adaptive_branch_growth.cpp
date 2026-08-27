#include "../verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"

#include <cmath>
#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

namespace
{

constexpr double tolerance = 5.0e-11;

qmomkit::RawMomentTensor raw_moments_from_terminals(const std::vector<double>& weights, const std::vector<std::vector<double>>& abscissas)
{
    const std::vector<std::size_t> extents{4, 4};
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat = 0; flat < moments.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, extents);
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            double contribution = weights[node];
            for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
            {
                contribution *= std::pow(abscissas[dimension][node], static_cast<int>(index[dimension]));
            }
            moments[flat] += contribution;
        }
    }
    return {extents, moments};
}

qmomkit::CqmomCase adaptive_case_from_moments(const qmomkit::RawMomentTensor& moments)
{
    std::ostringstream text;
    text.precision(17);
    text << "[phase_space]\n"
         << "names = x, y\n"
         << "node_counts = 2, 2\n"
         << "conditioning_order = 0, 1\n"
         << "adaptive_branch_counts = true\n\n"
         << "[qmom]\n"
         << "adaptive = true\n\n"
         << "[moments]\n"
         << "extents = 4, 4\n";
    for (std::size_t flat = 0; flat < moments.values.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, moments.extents);
        text << "M[" << index[0] << ',' << index[1] << "] = " << moments.values[flat] << '\n';
    }

    const auto parsed = qmomkit::config::parse_cqmom(text.str(), "adaptive-branch-growth");
    QMOMKIT_VERIFY_CHECK(parsed.ok());
    return parsed.value();
}

bool has_active_count(const std::vector<qmomkit::CqmomActiveBranchCount>& counts, const std::size_t level,
                      const qmomkit::ConditioningPath& parent_path, const std::size_t expected)
{
    for (const qmomkit::CqmomActiveBranchCount& count : counts)
    {
        if (count.conditioning_level == level && count.parent_path == parent_path && count.active_node_count == expected)
        {
            return true;
        }
    }
    return false;
}

bool has_growth_transition(const std::vector<qmomkit::CqmomNodeCountTransition>& transitions)
{
    for (const qmomkit::CqmomNodeCountTransition& transition : transitions)
    {
        if (transition.stage == qmomkit::PhaseSpaceEvolutionStage::post_source_inversion && transition.conditioning_level == 1 &&
            transition.parent_path == qmomkit::ConditioningPath{{1}} && transition.active_node_count_before == 1 &&
            transition.active_node_count_after == 2 && transition.kind == qmomkit::CqmomNodeCountTransitionKind::grew)
        {
            return true;
        }
    }
    return false;
}

void verify_adaptive_branch_reinjection()
{
    // The initial rule has one conditional y node for each x node.  The
    // prescribed final rule keeps the x marginal and splits only x = 1.
    const qmomkit::RawMomentTensor initial = raw_moments_from_terminals({0.2, 0.8}, {{-1.0, 1.0}, {0.0, 0.125}});
    const qmomkit::RawMomentTensor expected = raw_moments_from_terminals({0.2, 0.3, 0.5}, {{-1.0, 1.0, 1.0}, {0.0, -0.5, 0.5}});
    const qmomkit::CqmomState initial_state = adaptive_case_from_moments(initial).initial_state();

    QMOMKIT_VERIFY_CHECK(has_active_count(initial_state.active_branch_counts(), 0, {}, 2));
    QMOMKIT_VERIFY_CHECK(has_active_count(initial_state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 1));
    QMOMKIT_VERIFY_CHECK(has_active_count(initial_state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 1));

    const qmomkit::CqmomMomentSource source = [initial, expected](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t i = 0; i < rate.size(); ++i)
        {
            const std::size_t flat = qmomkit::flatten_multi_index(context.index(i), expected.extents);
            rate[i] = expected.values[flat] - initial.values[flat];
        }
        return qmomkit::CqmomMomentRate{rate};
    };

    const qmomkit::CqmomPhaseSpaceEvolutionResult result = qmomkit::advance_cqmom_moment_source(initial_state, 0.0, 1.0, source);
    QMOMKIT_VERIFY_CHECK(result.ok());
    QMOMKIT_VERIFY_CHECK(result.quadrature.has_value());
    QMOMKIT_VERIFY_CHECK(result.quadrature->total_nodes() == 3);
    const std::vector<qmomkit::ConditioningPath> expected_paths{qmomkit::ConditioningPath{{0, 0}}, qmomkit::ConditioningPath{{1, 0}},
                                                                qmomkit::ConditioningPath{{1, 1}}};
    QMOMKIT_VERIFY_CHECK(result.quadrature->terminal_paths() == expected_paths);
    QMOMKIT_VERIFY_CHECK(has_active_count(result.state.active_branch_counts(), 0, {}, 2));
    QMOMKIT_VERIFY_CHECK(has_active_count(result.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 1));
    QMOMKIT_VERIFY_CHECK(has_active_count(result.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 2));
    QMOMKIT_VERIFY_CHECK(has_growth_transition(result.diagnostics.cqmom_node_count_transitions));

    const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(*result.quadrature, expected.extents);
    QMOMKIT_VERIFY_CHECK(reconstructed.values.size() == expected.values.size());
    for (std::size_t flat = 0; flat < expected.values.size(); ++flat)
    {
        QMOMKIT_VERIFY_NEAR(reconstructed.values[flat], expected.values[flat], tolerance);
    }

    const std::vector<double>& weights = result.quadrature->weights();
    const std::vector<std::vector<double>>& abscissas = result.quadrature->abscissas();
    QMOMKIT_VERIFY_NEAR(weights[0], 0.2, tolerance);
    QMOMKIT_VERIFY_NEAR(weights[1], 0.3, tolerance);
    QMOMKIT_VERIFY_NEAR(weights[2], 0.5, tolerance);
    QMOMKIT_VERIFY_NEAR(abscissas[0][0], -1.0, tolerance);
    QMOMKIT_VERIFY_NEAR(abscissas[0][1], 1.0, tolerance);
    QMOMKIT_VERIFY_NEAR(abscissas[0][2], 1.0, tolerance);
    QMOMKIT_VERIFY_NEAR(abscissas[1][0], 0.0, tolerance);
    QMOMKIT_VERIFY_NEAR(abscissas[1][1], -0.5, tolerance);
    QMOMKIT_VERIFY_NEAR(abscissas[1][2], 0.5, tolerance);
}

} // namespace

int main()
{
    verify_adaptive_branch_reinjection();
    return 0;
}
