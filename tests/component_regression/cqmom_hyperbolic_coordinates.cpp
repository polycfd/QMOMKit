#include "../verification/verification_helpers.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace
{

using qmomkit::QmomFamily;

qmomkit::PhaseSpaceHierarchy hierarchy_of(std::vector<std::size_t> node_counts, std::vector<std::size_t> conditioning_order,
                                          const std::vector<QmomFamily>& families)
{
    qmomkit::HierarchyFamilies hierarchy_families;
    for (const QmomFamily family : families)
    {
        hierarchy_families.coordinates.push_back({family});
    }
    return {std::move(node_counts), std::move(conditioning_order), hierarchy_families};
}

// The raw moments of weighted points over the moment extents of the
// hierarchy. points[q] holds the coordinates of point q in phase-space order.
qmomkit::RawMomentTensor moments_of(const qmomkit::PhaseSpaceHierarchy& hierarchy, const std::vector<double>& weights,
                                    const std::vector<std::vector<double>>& points)
{
    const std::vector<std::size_t> extents = hierarchy.moment_extents();
    std::vector<double> values(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat = 0; flat < values.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, extents);
        for (std::size_t point = 0; point < weights.size(); ++point)
        {
            double contribution = weights[point];
            for (std::size_t dimension = 0; dimension < extents.size(); ++dimension)
            {
                contribution *= std::pow(points[point][dimension], static_cast<int>(index[dimension]));
            }
            values[flat] += contribution;
        }
    }
    return {extents, values};
}

qmomkit::CqmomCase case_of(const qmomkit::PhaseSpaceHierarchy& hierarchy, const qmomkit::RawMomentTensor& moments, const bool adaptive)
{
    qmomkit::CqmomCase cqmom_case{qmomkit::extract_cqmom_moments(moments, qmomkit::create_cqmom_moment_layout(hierarchy))};
    cqmom_case.qmom_options.adaptive = adaptive;
    return cqmom_case;
}

// The constant rate that carries the moments of the state onto the target
// over a source step of 1.
qmomkit::CqmomMomentSource source_from_to(const qmomkit::CqmomState& state, const qmomkit::RawMomentTensor& target)
{
    return [state, target](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t i = 0; i < rate.size(); ++i)
        {
            const qmomkit::MultiIndexView index = context.index(i);
            rate[i] = target.values[qmomkit::flatten_multi_index(index, target.extents)] - state.moment(index).value();
        }
        return qmomkit::CqmomMomentRate{rate};
    };
}

std::size_t active_count(const qmomkit::CqmomState& state, const std::size_t level, const qmomkit::ConditioningPath& parent_path)
{
    for (const qmomkit::CqmomActiveBranchCount& count : state.active_branch_counts())
    {
        if (count.conditioning_level == level && count.parent_path == parent_path)
        {
            return count.active_node_count;
        }
    }
    QMOMKIT_VERIFY_CHECK(false);
    return 0;
}

void check_identical(const qmomkit::CqmomResult& actual, const qmomkit::CqmomResult& expected)
{
    QMOMKIT_VERIFY_CHECK(actual.status == expected.status);
    QMOMKIT_VERIFY_CHECK(actual.message == expected.message);
    QMOMKIT_VERIFY_CHECK(actual.quadrature.weights() == expected.quadrature.weights());
    QMOMKIT_VERIFY_CHECK(actual.quadrature.abscissas() == expected.quadrature.abscissas());
    QMOMKIT_VERIFY_CHECK(actual.quadrature.terminal_paths() == expected.quadrature.terminal_paths());
    QMOMKIT_VERIFY_CHECK(actual.quadrature.has_rectangular_storage() == expected.quadrature.has_rectangular_storage());
}

void check_identical(const qmomkit::CqmomState& actual, const qmomkit::CqmomState& expected)
{
    QMOMKIT_VERIFY_CHECK(actual.moment_values() == expected.moment_values());
    QMOMKIT_VERIFY_CHECK(actual.quadrature().weights() == expected.quadrature().weights());
    QMOMKIT_VERIFY_CHECK(actual.quadrature().abscissas() == expected.quadrature().abscissas());
    QMOMKIT_VERIFY_CHECK(actual.active_branch_counts().size() == expected.active_branch_counts().size());
    for (std::size_t branch = 0; branch < expected.active_branch_counts().size(); ++branch)
    {
        const qmomkit::CqmomActiveBranchCount& lhs = actual.active_branch_counts()[branch];
        const qmomkit::CqmomActiveBranchCount& rhs = expected.active_branch_counts()[branch];
        QMOMKIT_VERIFY_CHECK(lhs.conditioning_level == rhs.conditioning_level && lhs.parent_path == rhs.parent_path &&
                             lhs.configured_node_count == rhs.configured_node_count && lhs.active_node_count == rhs.active_node_count);
    }
}

// The hierarchy has three hyperbolic coordinates with node count 2 and
// conditioning order {2, 0, 1}. The parent is phase-space dimension 2, the
// middle coordinate, which is both a child and a parent, is dimension 0, and
// the last coordinate is dimension 1. The points are a product of three values
// of each of dimensions 2 and 0 and four values of dimension 1, so the moments
// of every coordinate are interior. The values of dimension 1 are multiplied
// by spread and shifted by a quarter of dimension 2.
qmomkit::PhaseSpaceHierarchy three_coordinate_hierarchy()
{
    return hierarchy_of({2, 2, 2}, {2, 0, 1}, {QmomFamily::hyperbolic, QmomFamily::hyperbolic, QmomFamily::hyperbolic});
}

qmomkit::RawMomentTensor three_coordinate_moments(const qmomkit::PhaseSpaceHierarchy& hierarchy, const double spread)
{
    std::vector<double> weights;
    std::vector<std::vector<double>> points;
    for (const double z : {-1.0, 0.5, 2.0})
    {
        for (const double x : {-0.5, 0.25, 1.0})
        {
            for (const double y : {-1.0, -0.25, 0.5, 1.5})
            {
                weights.push_back(1.0 / 36.0);
                points.push_back({x, spread * y + 0.25 * z, z});
            }
        }
    }
    return moments_of(hierarchy, weights, points);
}

void test_branch_grows_one_node_per_source_step_to_the_hyperbolic_rule()
{
    // The coordinate x is classical with node count 1, and y is hyperbolic
    // with node count 2. The state starts from one point, so the branch of y
    // holds one node. Each source step carries the moments onto those of four
    // points with the same x, which are interior. The branch grows by one node
    // per step, from one node to the classical rule of two nodes and then to
    // the three nodes of HyQMOM, and it never jumps from one node to three. The
    // four points are not symmetric about their mean, so the classical rule of
    // two nodes differs from the two nodes of HyQMOM with node count 1.
    const qmomkit::PhaseSpaceHierarchy hierarchy = hierarchy_of({1, 2}, {0, 1}, {QmomFamily::classical, QmomFamily::hyperbolic});
    const qmomkit::RawMomentTensor four_points = moments_of(hierarchy, {0.25, 0.25, 0.25, 0.25}, {{0.5, -1.0}, {0.5, 0.0}, {0.5, 0.5}, {0.5, 2.0}});
    qmomkit::CqmomState state = qmomkit::CqmomPreparedCase{case_of(hierarchy, moments_of(hierarchy, {1.0}, {{0.5, 0.25}}), true)}.initial_state();
    QMOMKIT_VERIFY_CHECK(active_count(state, 1, qmomkit::ConditioningPath{{0}}) == 1);

    // The marginal moments of y are its conditional moments at the one node
    // of x, whose weight is 1.
    const qmomkit::QmomResult classical_two_nodes =
        qmomkit::compute_qmom(qmomkit::MomentSet({four_points.values[0], four_points.values[1], four_points.values[2], four_points.values[3]}), 2);
    QMOMKIT_VERIFY_CHECK(classical_two_nodes.ok() && classical_two_nodes.weights.size() == 2);
    for (std::size_t step = 1; step <= 4; ++step)
    {
        const qmomkit::CqmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_cqmom_moment_source(state, static_cast<double>(step - 1), 1.0, source_from_to(state, four_points));
        QMOMKIT_VERIFY_CHECK(result.ok());
        const std::size_t expected = std::min<std::size_t>(1 + step, 3);
        QMOMKIT_VERIFY_CHECK(active_count(result.state, 1, qmomkit::ConditioningPath{{0}}) == expected);
        QMOMKIT_VERIFY_CHECK(result.state.quadrature().total_nodes() == expected);
        if (expected == 2)
        {
            for (std::size_t node = 0; node < 2; ++node)
            {
                QMOMKIT_VERIFY_NEAR(result.state.quadrature().weights()[node], classical_two_nodes.weights[node], 1.0e-12);
                QMOMKIT_VERIFY_NEAR(result.state.quadrature().abscissas()[1][node], classical_two_nodes.abscissas[node], 1.0e-12);
            }
        }
        state = result.state;
    }
    const qmomkit::CqmomResult hyperbolic = qmomkit::compute_cqmom(case_of(hierarchy, four_points, true));
    QMOMKIT_VERIFY_CHECK(hyperbolic.ok() && hyperbolic.quadrature.total_nodes() == 3);
    for (std::size_t node = 0; node < 3; ++node)
    {
        QMOMKIT_VERIFY_NEAR(state.quadrature().weights()[node], hyperbolic.quadrature.weights()[node], 1.0e-12);
        QMOMKIT_VERIFY_NEAR(state.quadrature().abscissas()[1][node], hyperbolic.quadrature.abscissas()[1][node], 1.0e-12);
    }
}

void test_open_and_prepared_cases_agree_for_hyperbolic_coordinates()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy = three_coordinate_hierarchy();
    for (const bool adaptive : {false, true})
    {
        const qmomkit::CqmomCase cqmom_case = case_of(hierarchy, three_coordinate_moments(hierarchy, 1.0), adaptive);
        const qmomkit::CqmomPreparedCase prepared_case{cqmom_case};
        const qmomkit::CqmomResult open = qmomkit::compute_cqmom(cqmom_case);
        QMOMKIT_VERIFY_CHECK(open.ok());
        QMOMKIT_VERIFY_CHECK(open.quadrature.total_nodes() == 27);
        check_identical(qmomkit::compute_cqmom(prepared_case), open);
        check_identical(prepared_case.initial_state(), cqmom_case.initial_state());
    }
}

void test_moment_replacement_for_hyperbolic_coordinates()
{
    // Replacing the moments of a case gives the case built from the new
    // moments, through the tensor and the moment data overloads.
    const qmomkit::PhaseSpaceHierarchy hierarchy = three_coordinate_hierarchy();
    const qmomkit::RawMomentTensor replacement = three_coordinate_moments(hierarchy, 0.5);
    for (const bool adaptive : {false, true})
    {
        const qmomkit::CqmomPreparedCase base{case_of(hierarchy, three_coordinate_moments(hierarchy, 1.0), adaptive)};
        const qmomkit::CqmomCase expected_case = case_of(hierarchy, replacement, adaptive);
        const qmomkit::CqmomPreparedCase expected{expected_case};
        const qmomkit::CqmomResult expected_result = qmomkit::compute_cqmom(expected);
        QMOMKIT_VERIFY_CHECK(expected_result.ok());

        const auto from_tensor = qmomkit::config::with_initial_moments(base, replacement);
        const auto from_data = qmomkit::config::with_initial_moments(base, expected_case.moments);
        QMOMKIT_VERIFY_CHECK(from_tensor.ok() && from_data.ok());
        for (const qmomkit::CqmomPreparedCase* replaced : {&from_tensor.value(), &from_data.value()})
        {
            check_identical(qmomkit::compute_cqmom(*replaced), expected_result);
            check_identical(replaced->initial_state(), expected.initial_state());
        }
    }
}

void test_fixed_hyperbolic_failure_propagates_as_qmom_failed()
{
    // The coordinate x is classical at two points, and y is hyperbolic with
    // node count 2. The moments of four points of y at each x are interior. When the branch at
    // x = 1.5 holds only two points, its conditional moments lie on the
    // boundary, and fixed mode fails in y in compute_cqmom, state_candidate and
    // the post-source reconstruction of a source step.
    const qmomkit::PhaseSpaceHierarchy hierarchy = hierarchy_of({2, 2}, {0, 1}, {QmomFamily::classical, QmomFamily::hyperbolic});
    const std::vector<double> ys{-1.0, 0.0, 0.5, 2.0};
    std::vector<double> weights;
    std::vector<std::vector<double>> points;
    for (const double x : {-0.5, 1.5})
    {
        for (const double y : ys)
        {
            weights.push_back(0.125);
            points.push_back({x, y});
        }
    }
    const qmomkit::RawMomentTensor interior = moments_of(hierarchy, weights, points);
    const qmomkit::RawMomentTensor boundary = moments_of(hierarchy, {0.125, 0.125, 0.125, 0.125, 0.125, 0.375},
                                                         {{-0.5, -1.0}, {-0.5, 0.0}, {-0.5, 0.5}, {-0.5, 2.0}, {1.5, -1.0}, {1.5, 2.0}});

    const qmomkit::CqmomResult failed = qmomkit::compute_cqmom(case_of(hierarchy, boundary, false));
    QMOMKIT_VERIFY_CHECK(failed.status == qmomkit::CqmomStatus::qmom_failed);
    QMOMKIT_VERIFY_CHECK(failed.failed_phase_space_dimension == 1);

    const qmomkit::CqmomCase interior_case = case_of(hierarchy, interior, false);
    const qmomkit::CqmomPreparedCase prepared_case{interior_case};
    const qmomkit::CqmomStateCandidate candidate = prepared_case.state_candidate(case_of(hierarchy, boundary, false).moments);
    QMOMKIT_VERIFY_CHECK(candidate.status == qmomkit::StateCandidateStatus::reconstruction_failed);
    QMOMKIT_VERIFY_CHECK(candidate.reconstruction_status == qmomkit::CqmomStatus::qmom_failed);
    QMOMKIT_VERIFY_CHECK(candidate.message == failed.message);

    const qmomkit::CqmomState state = prepared_case.initial_state();
    QMOMKIT_VERIFY_CHECK(state.quadrature().total_nodes() == 6);
    const qmomkit::CqmomPhaseSpaceEvolutionResult evolved = qmomkit::advance_cqmom_moment_source(state, 0.0, 1.0, source_from_to(state, boundary));
    QMOMKIT_VERIFY_CHECK(evolved.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed);
    QMOMKIT_VERIFY_CHECK(evolved.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::post_source_inversion);
    QMOMKIT_VERIFY_CHECK(evolved.diagnostics.post_source_cqmom_status == qmomkit::CqmomStatus::qmom_failed);
    QMOMKIT_VERIFY_CHECK(evolved.state.moment_values() == state.moment_values());
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"branch_grows_one_node_per_source_step_to_the_hyperbolic_rule", test_branch_grows_one_node_per_source_step_to_the_hyperbolic_rule},
            {"open_and_prepared_cases_agree_for_hyperbolic_coordinates", test_open_and_prepared_cases_agree_for_hyperbolic_coordinates},
            {"moment_replacement_for_hyperbolic_coordinates", test_moment_replacement_for_hyperbolic_coordinates},
            {"fixed_hyperbolic_failure_propagates_as_qmom_failed", test_fixed_hyperbolic_failure_propagates_as_qmom_failed},
        });
}
